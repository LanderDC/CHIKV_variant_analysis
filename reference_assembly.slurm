#!/bin/bash -l
#SBATCH --job-name="ref_assembly"
#SBATCH --nodes="1"
#SBATCH -M genius
#SBATCH -c 36
#SBATCH --time="1-00:00:00"
#SBATCH --account="lp_jm_virome_group"

#while read line; do
#sbatch --cluster=wice -o $line.log --export=ALL,study=CHIKV_Sam,sample=$line viper.slurm
#done < names.txt

threads=36

#Move to wherever your files are
cd ~/staging/LDC/CHIKV_Sam

#Activate conda environment
conda activate viper

#Make sure programs in conda viper environment are found first in PATH
PATH="$CONDA_PREFIX/bin:$PATH"

trimmomatic PE \
		-threads "$threads" \
		$sample.R1.fastq.gz $sample.R2.fastq.gz \
		"$sample".trimmed.R1.fastq.gz "$sample".R1.unpaired.fastq.gz \
		"$sample".trimmed.R2.fastq.gz "$sample".R2.unpaired.fastq.gz \
		HEADCROP 10 \ # Not performed
		ILLUMINACLIP:"$CONDA_PREFIX/share/trimmomatic/adapters/TruSeq3-PE.fa":2:30:7:1:true \
		LEADING:15 \
		TRAILING:15 \
		SLIDINGWINDOW:4:20 \
		MINLEN:50

cat "$sample".R1.unpaired.fastq.gz "$sample".R2.unpaired.fastq.gz > "$sample".trimmed.unpaired.fastq.gz
rm "$sample".R1.unpaired.fastq.gz
rm "$sample".R2.unpaired.fastq.gz

reference_assembly.sh -o $sample -1 $VSC_SCRATCH/$study/"$sample".trimmed.R1.fastq.gz -2 $VSC_SCRATCH/$study/"$sample".trimmed.R2.fastq.gz -u $VSC_SCRATCH/$study/"$sample".trimmed.unpaired.fastq.gz -i BNI-CHIKV_899_FJ959103.fasta -t $threads

cd "$sample"

# Remove PCR duplicates
samtools sort -n -@ 5 "$sample".sorted.bam | samtools fixmate -m -@ 5 - - | samtools sort -@ 5 - | samtools markdup -r -s -@ 5 - "$sample".deduplicated.bam

# Consensus calling
samtools consensus -f fasta -A --min-depth 5 -o $sample.fasta -a -C 10 --threads $threads $sample.deduplicated.bam
sed -i "s/>.*/&_$sample/" $sample.fasta

rm BNI*
rm *.sam
