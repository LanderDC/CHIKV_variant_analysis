library(tidyverse)

ctrl <- rgb(215, 199, 179, maxColorValue = 255)
sof80 <-  rgb(72, 151, 151, maxColorValue = 255)
mpv10 <-  rgb(208, 142, 166, maxColorValue = 255)
combo <-  rgb(108, 132, 193, maxColorValue = 255)
mpv50 <-  rgb(183, 79, 115, maxColorValue = 255)
mpv100 <-  rgb(140, 57, 86, maxColorValue = 255)

cols <- c("CTRL"=ctrl, "SOF80"=sof80, "MPV10"=mpv10, "MPV50"=mpv50, "MPV100"=mpv100, "COMBO"=combo)
levels <- c("CTRL", "SOF80", "MPV10", "MPV50", "MPV100", "COMBO")

names_df <- read_tsv("data/mapped_reads_clean.tsv")

names_df <- names_df |> 
  separate(long_name, into = c("virus", "AG", "type", "concentration"), remove=F) |> 
  select(short_name, long_name, type, concentration, mapped_reads, deduplicated_reads) |> 
  mutate(concentration=ifelse(concentration == "Right", NA, concentration))

variant_df <- read_tsv("data/CHIKV_Sam_depth100.tsv", col_names = T) |> 
  rename_with(~ gsub("\\.tsv", "", .))

variant_df <- variant_df |> 
  select(-REGION, -GFF_FEATURE)

mutation_df <- variant_df |> 
  select(POS, REF, ALT, starts_with("PASS_")) |> 
  pivot_longer(cols=starts_with("PASS"), names_to = "sample", values_to = "real_variant") |> 
  mutate(sample=gsub("PASS_", "", sample),
        ALT=case_when(startsWith(ALT, "+") ~ "insertion",
                      startsWith(ALT, "-") ~ "deletion",
                      T ~ ALT),
        mutation=case_when(ALT %in% c("insertion", "deletion") ~ ALT,
                            T ~ glue::glue("{REF} \u2192 {ALT}"))) |> 
  filter(real_variant == T) |> 
  left_join(names_df, by=join_by("sample" == "short_name")) |> 
  select(long_name, type, concentration, POS, mutation, deduplicated_reads) |> 
  rename(sample=long_name)

# Identify POS and mutation combinations present in CTRL
ctrl_combinations <- mutation_df %>%
  select(-deduplicated_reads) |> 
  filter(type == "CTRL") %>%
  distinct() |> 
  mutate(value=T) |> 
  pivot_wider(names_from="sample", values_from = "value") |>
  # filter to remove mutations from list if they are not in all controls
  filter(if_all(starts_with("CHIKV"), ~ . == TRUE)) |> 
  # filter to remove mutations from list if they are in >=80% of the controls
  #filter(
  #  rowSums(across(starts_with("CHIKV"), ~ . == TRUE), na.rm = TRUE) >= 
  #    0.8 * ncol(across(starts_with("CHIKV")))
  #) %>%
  select(POS, mutation) %>%
  distinct()

# Filter out rows that match the CTRL combinations
filtered_df <- mutation_df %>%
  anti_join(ctrl_combinations, by = c("POS", "mutation")) |> 
  group_by(mutation, sample, type, concentration, deduplicated_reads) |> 
  summarise(n=n(), .groups = "drop")

final_df <- filtered_df |> 
  pivot_wider(names_from = mutation, values_from = n) |> 
  mutate_if(is.numeric, ~ replace_na(., 0)) |> 
  pivot_longer(cols=c(-sample, -type, -concentration, -deduplicated_reads), names_to = "mutation", values_to = "n") |>   
  mutate(type=ifelse(!is.na(concentration), paste0(type, concentration), type),
         type=type |> fct_relevel(levels),
         sample=sample |> fct_relevel(levels(gtools::mixedsort(unique(final_df$sample)))),
         mutation=mutation |> fct_relevel(c("deletion", "insertion", setdiff(rev(unique(filtered_df$mutation)), c("insertion", "deletion")))),
         normalized_mutations=n/deduplicated_reads*1000000)

heatmap <- final_df |> 
  ggplot(aes(x=sample, y=mutation, fill=n))+
  geom_tile()+
  scale_fill_gradient(low = "white", high = "blue")+
    facet_grid(~ type, scales = "free_x", space='free') +
  labs(fill="Mutations")+
  theme_minimal()+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        #plot.margin = margin(l = 20, unit = "mm"),
        panel.grid.major = element_blank(),)
ggsave("figures/mutation_count_samples.pdf", dpi=300, width=7, height=4.5, device = cairo_pdf)

#mean_df <- final_df |> 
#  group_by(type, concentration, mutation) |> 
#  summarise(mean=mean(n)) |> 
#  mutate(type=ifelse(!is.na(concentration), paste(type, concentration), type))
#
#
#mean_df$type <- factor(mean_df$type, levels = (c("CTRL", "SOF 80", "MPV 10", "MPV 50", "MPV 100", "COMBO")))
#
#mean_df |> 
#  ggplot(aes(x=type, y=mutation, fill=mean))+
#  geom_tile()+
#  scale_fill_gradient(low = "white", high = "blue")+
#  theme_minimal()+
#  theme(axis.text.x = element_text(angle = 45, hjust=1),
#        panel.grid.major = element_blank())
#ggsave("figures/mutation_count_compounds.pdf", dpi=300, width=5, device = cairo_pdf)

# Boxplots
total_df <- final_df |> 
  #left_join(names_df, by=join_by("sample" == "long_name")) |> 
  group_by(sample, type, concentration, deduplicated_reads) |> 
  summarise(total_mutations=sum(n), .groups = "drop") |> 
  mutate(normalized_mutations=total_mutations/deduplicated_reads*1000000) |> 
  select(-concentration)

# Wilcoxon tests
wtest <- total_df |> 
  # Filter outlier
  filter(sample != "CHIKV_AG_SOF_80_Right_ankle_50") |> 
  #filter(type != "CTRL") |> 
  mutate(type=as.character(type)) |> 
  rstatix::wilcox_test(normalized_mutations ~ type, 
                       alternative = "less",
                       p.adjust.method = "BH",
                       comparisons = list(#c("CTRL", "SOF80"), c("SOF80", "COMBO"), c("SOF80", "MPV10"),
                       c("COMBO", "MPV10"),c("COMBO", "MPV50"), c("COMBO", "MPV100")))
                       #ref.group = "COMBO")

wtest |> 
  select(-.y.) |> 
  write_tsv("results/wilcoxon_test_combo_mpv.tsv")

wtest <- wtest |> 
  cbind(tibble(xmin=c(6, 6, 6), xmax=c(3, 4, 5), y.position=4))

wtest2 <- total_df |> 
  # Filter outlier
  filter(sample != "CHIKV_AG_SOF_80_Right_ankle_50") |> 
  #filter(type != "CTRL") |> 
  mutate(type=as.character(type)) |> 
  rstatix::wilcox_test(normalized_mutations ~ type, 
                       alternative = "less",
                       p.adjust.method = "BH",
                       comparisons = list(c("CTRL", "SOF80")))
                       #ref.group = "COMBO")

wtest2 <- wtest2 |> 
  cbind(tibble(xmin=c(1), xmax=c(2), y.position=4))

sp <- c(0, 0)

boxplot <- total_df |> 
  ggplot(aes(x=type, y=normalized_mutations))+
  geom_boxplot(aes(color=type, fill=type), outliers = F)+
  geom_point(aes(color=type), size=1)+
  #ggpubr::geom_pwc(
  #  method = "wilcox_test",
  #  method.args=list(alternative="less"),
  #  p.adjust.method = "BH", 
  #  ref.group = "COMBO",
  #  label = "p.adj.format", hide.ns = "p.adj", show.legend = F, tip.length = 0.01, label.size=2
  #)+
  ggpubr::stat_pvalue_manual(wtest, label="p.adj.signif", step.increase = .1, tip.length = 0.01, label.size=c(2, sp, 4, sp, 4, sp), hide.ns = F)+
  ggpubr::stat_pvalue_manual(wtest2, label="p.adj.signif", step.increase = .1, tip.length = 0.01, label.size=c(2), hide.ns = F)+
  labs(y="Normalized mutation count")+
  scale_color_manual(values=cols)+
  scale_fill_manual(values=alpha(cols, 0.1))+
  scale_y_log10()+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        legend.position = "none")
boxplot
ggsave("figures/boxplots_total_mutations.pdf", dpi=300, width=7, height=5, device = cairo_pdf)

# Barplot
transitions <- c("A \u2192 G", "G \u2192 A", "C \u2192 T", "T \u2192 C")

barplot <- final_df |> 
  filter(sample != "CHIKV_AG_SOF_80_Right_ankle_50") |> 
  group_by(mutation, type) |> 
  summarise(avg=mean(n), .groups="drop") |>
  left_join(
    final_df |> 
      filter(type == "CTRL") |> 
      group_by(mutation) |> 
      summarise(ctrl_avg = mean(n), .groups = "drop"),
    by = "mutation"
  ) |> 
  mutate(scaled_avg = avg / ctrl_avg, 
         mutation=mutation |> fct_relevel(c(transitions, setdiff(unique(final_df$mutation), c("deletion", transitions)), "deletion"))) |> 
  select(mutation, type, scaled_avg) |> 
  filter(type!="CTRL") |> 
  ggplot(aes(y = scaled_avg, x = mutation, fill = type)) +
  geom_col(position = "dodge")+
  geom_hline(yintercept = 1, linetype=3, color = "black") +
  geom_vline(xintercept = 4.5, linetype = "dashed")+
  geom_vline(xintercept = 12.5, linetype = "dashed")+
  annotate("text", x=2.5, y=6, label = "Transitions", size = 2)+
  annotate("text", x=8.5, y=6, label = "Transversions", size = 2)+
  annotate("text", x=13.5, y=6, label = "Indels", size = 2)+
  labs(y="Fold change \n (mean mutation count scaled to control group)")+
  scale_fill_manual(values=cols)+
  #scale_x_discrete(limits=rev)+
  scale_y_continuous(n.breaks = 6, limits=c(0,6), expand = expansion(add=c(0, 0.2)))+
  theme_classic()+
  theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank(),
        legend.title=element_blank())
ggsave("figures/scaled_barplot.pdf", dpi=300, device=cairo_pdf, width=7, height=5)
 
library(patchwork)
(heatmap+theme(legend.key.size = unit(3, "mm"),
               legend.text = element_text(size=5),
               legend.title = element_text(size=6),
               legend.position = "right",
               legend.margin = margin(l = -10),
               axis.text = element_text(size=6.5),
               text=element_text(size=9))) / 
  (free(barplot+theme(legend.position = "bottom",
                      legend.key.size = unit(2, "mm"),
                      legend.text = element_text(size=5),
                      legend.margin = margin(t = -10),
                      axis.title=element_text(size=7),
                      axis.text = element_text(size=6.5))) | 
           free(boxplot+theme(axis.title=element_text(size=7), 
                              axis.text = element_text(size=6.5))))+
                                plot_annotation(tag_levels = 'A') & 
                                theme(plot.tag = element_text(size = 10, face="bold"))
ggsave("figures/FigureX.pdf", dpi=300, width=7, height=6.5, device=cairo_pdf)

(((heatmap+theme(legend.key.size = unit(3, "mm"),
               legend.text = element_text(size=5),
               legend.title = element_text(size=6),
               legend.position = "right",
               legend.margin = margin(l = -10),
               legend.key.spacing = unit(0.5, "mm"),
               panel.spacing = unit(.75, "mm"),
               axis.text = element_text(size=6.5),
               text=element_text(size=9)) |
  free(boxplot+theme(axis.title=element_text(size=7), 
                     axis.text.x = element_text(size=6.5, angle=45, hjust = 1, vjust=1),
                     axis.text.y = element_text(size=6.5))))+
  plot_layout(widths=c(1, .6))) /
  free(barplot+theme(legend.position = "bottom",
                      legend.key.size = unit(3, "mm"),
                      legend.text = element_text(size=6, margin = margin(l=0.5, unit = "mm")),
                      legend.margin = margin(t = -10),
                      axis.title = element_text(size=7),
                      axis.text = element_text(size=6.5),
                      axis.text.x = element_text(angle=45, hjust = 1, vjust=1))))+
  plot_layout(heights = c(.6, 1))+
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 10, face="bold"))
ggsave("figures/FigureX-b.pdf", dpi=300, width=7, height=6.5, device=cairo_pdf)

barplot_df <- final_df |> 
  filter(sample != "CHIKV_AG_SOF_80_Right_ankle_50") |> 
  group_by(mutation, type) |> 
  mutate(mutation_sum=sum(n)) |> 
  ungroup() |> 
  group_by(type) |> 
  mutate(total_sum=sum(n), nr_samples=length(unique(sample))) |> 
  ungroup() |> 
  mutate(relative=(mutation_sum/total_sum*100),
         mutation=mutation |> fct_relevel(c(transitions, setdiff(unique(final_df$mutation), c("deletion", transitions)), "deletion"))) |> 
  select(type, mutation, relative) |> 
  distinct()

barplot_df |> 
  #left_join(
  #  barplot_df |> 
  #    filter(type == "CTRL"),
  #  by = "mutation"
  #) |> 
  #filter(type.x != "CTRL") |> 
  ggplot(aes(x=mutation, y=relative, fill=type))+
  geom_col(position = "dodge")+
  geom_hline(yintercept = 1, linetype=3, color = "black") +
  geom_vline(xintercept = 4.5, linetype = "dashed")+
  geom_vline(xintercept = 12.5, linetype = "dashed")+
  annotate("text", x=2.5, y=20, label = "Transitions", size = 3)+
  annotate("text", x=8.5, y=20, label = "Transversions", size = 3)+
  annotate("text", x=13.5, y=20, label = "Indels", size = 3)+
  labs(y="Relative number of mutations (%)")+
  scale_fill_manual(values=cols)+
  #scale_x_discrete(limits=rev)+
  scale_y_continuous(n.breaks = 5, expand = expansion(add=c(0, 0.2)))+
  theme_classic()+
  theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank(),
        legend.title=element_blank())
ggsave("figures/relative_barplot.pdf", dpi=300, device=cairo_pdf, width=7, height=5)

barplot_df |> 
  filter(type=="SOF80") |> 
  group_by(type) |> 
  distinct() |> 
  summarise(sum(relative))

# Non-synonymous
# ==> no difference
ns_df <- variant_df |> 
  select(POS, REF, ALT, REF_AA, ALT_AA, starts_with("PASS_")) |> 
  pivot_longer(cols=starts_with("PASS"), names_to = "sample", values_to = "real_variant") |> 
  mutate(sample=gsub("PASS_", "", sample),
        ALT=case_when(startsWith(ALT, "+") ~ "insertion",
                      startsWith(ALT, "-") ~ "deletion",
                      T ~ ALT),
        mutation=case_when(ALT %in% c("insertion", "deletion") ~ ALT,
                            T ~ glue::glue("{REF} \u2192 {ALT}"))) |> 
  #Filter on real variants and  non-synonymous mutations
  filter(real_variant == T, REF_AA != ALT_AA) |> 
  left_join(names_df, by=join_by("sample" == "short_name")) |> 
  select(long_name, type, concentration, POS, mutation, deduplicated_reads) |> 
  rename(sample=long_name)

ns_filtered <- ns_df %>%
  anti_join(ctrl_combinations, by = c("POS", "mutation")) |> 
  group_by(mutation, sample, type, concentration, deduplicated_reads) |> 
  summarise(n=n(), .groups = "drop")

final_ns <- ns_filtered |> 
  pivot_wider(names_from = mutation, values_from = n) |> 
  mutate_if(is.numeric, ~ replace_na(., 0)) |> 
  pivot_longer(cols=c(-sample, -type, -concentration, -deduplicated_reads), names_to = "mutation", values_to = "n") |>   
  mutate(type=ifelse(!is.na(concentration), paste0(type, concentration), type),
         type=type |> fct_relevel(levels),
         sample=sample |> fct_relevel(levels(gtools::mixedsort(unique(final_df$sample)))),
         mutation=mutation |> fct_relevel(c("deletion", "insertion", setdiff(rev(unique(ns_filtered$mutation)), c("insertion", "deletion")))),
         normalized_mutations=n/deduplicated_reads*1000000)

final_ns |> 
  ggplot(aes(x=sample, y=mutation, fill=n))+
  geom_tile()+
  scale_fill_gradient(low = "white", high = "blue")+
    facet_grid(~ type, scales = "free_x", space='free') +
  labs(fill="Mutations")+
  theme_minimal()+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        #plot.margin = margin(l = 20, unit = "mm"),
        panel.grid.major = element_blank(),)
