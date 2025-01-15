#!/bin/bash -l
#SBATCH --job-name="variant_calling"
#SBATCH --nodes="1"
#SBATCH -M genius
#SBATCH -c 5
#SBATCH --time="12:00:00"
#SBATCH --account="lp_jm_virome_group"

#Move to wherever your files are
cd ~/staging/LDC/CHIKV_Sam

#Activate conda environment
conda activate viper

#Make sure programs in conda viper environment are found first in PATH
PATH="$CONDA_PREFIX/bin:$PATH"

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/lustre1/project/stg_00029/LDC/software/miniforge3-genius/lib

> ivar_files_depth100.txt

mkdir depth100
cd depth100

while read line; do

samtools mpileup -aa -A -d 0 -B -Q 0 --reference ../BNI-CHIKV_899_FJ959103.fasta ../$line/$line.deduplicated.bam | ivar variants -p $line -m 100 -r ../BNI-CHIKV_899_FJ959103.fasta -g ../FJ959103.gff

echo "$line.tsv" >> ../ivar_files_depth100.txt

done < ../names.txt

ivar filtervariants -p CHIKV_Sam_depth100 -f ../ivar_files_depth100.txt -t 0
