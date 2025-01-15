#!/bin/bash

# Function to display help message
display_help() {
    echo "Usage: ./script.sh -o output_dir -1 read1.fastq -2 read2.fastq -u unpaired.fastq -i reference.fasta -t threads"
    echo "Options:"
    echo "  -o <output_dir>: Output directory"
    echo "  -1 <read1.fastq>: Paired-end read 1"
    echo "  -2 <read2.fastq>: Paired-end read 2"
    echo "  -u <unpaired.fastq>: Unpaired reads"
    echo "  -i <reference.fasta>: Reference genome in FASTA format"
    echo "  -t <threads>: Number of threads"
}

# Set a flag to track if any known options are provided
options_given=0

while getopts o:1:2:u:i:t: option; do
    options_given=1  # Update the flag to indicate that at least one known option is provided
    case "${option}" in
        o) output=${OPTARG} ;;
        1) read1=${OPTARG} ;;
        2) read2=${OPTARG} ;;
        u) unpaired=${OPTARG} ;;
        i) input=${OPTARG} ;;
        t) threads=${OPTARG} ;;
        *) display_help
           exit 1 ;;
    esac
done

# Check if no known options are provided
if [ $options_given -eq 0 ]; then
    display_help
    exit 1
fi

mkdir -p $output
cp $input $output
cd $output

bowtie2-build $input $input

echo
bowtie2 -x $input -1 $read1 -2 $read2 -p $threads -S $output.paired.sam
bowtie2 -x $input -U $unpaired -p $threads -S $output.unpaired.sam

samtools view -F 4 -bS $output.paired.sam -@ "$threads" | samtools sort - > $output.paired_sorted.bam
samtools view -F 4 -bS $output.unpaired.sam -@ "$threads" | samtools sort - > $output.unpaired_sorted.bam

samtools merge -f $output.bam $output.paired_sorted.bam $output.unpaired_sorted.bam

samtools sort $output.bam > $output.sorted.bam

rm $output.paired_sorted.bam
rm $output.unpaired_sorted.bam
rm $output.bam

samtools consensus -f fasta -A --min-depth 5 -o $output.fasta -a -C 10 --threads "$threads" $output.sorted.bam
sed -i "s/>.*/&_$output/" $output.fasta

rm ${input}*
rm *.sam
