#!/bin/bash -l
#SBATCH -D /home/twinkle1/projects/seed_color_differential_expression_analysis
#SBATCH -t 2:00:00
#SBATCH -J STAR
#SBATCH -o logs/STAR/dea/mappingLog-%j.txt
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=8gb
#SBATCH --array 0-25

module load star/2.7.8a

# there are a total of 26 samples, array goes from 0-25


### MAIN

# create array of read fastq files (R1 only):
SOURCE_DIR=/scratch/twinkle1/rnaSeq_STAR_input
FILES=("$SOURCE_DIR"/*trimmed_1P.fq.gz)

INSTRING=/scratch/twinkle1/rnaSeq_STAR_input/
OUTSTRING=/scratch/twinkle1/rnaSeq_STAR_output/

# only for testing puposes:
#echo "${FILES[0]}"
#echo "${FILES[0]/R1.fastq.gz/R2.fastq.gz}"
#echo "${FILES[0]/$INSTRING/$OUTSTRING}"


# run STAR after genome index creation
mkdir -p /scratch/twinkle1/rnaSeq_STAR_output

STAR --runThreadN 8 \
	--runMode alignReads \
	--readFilesCommand zcat \
	--outSAMtype BAM SortedByCoordinate \
	--genomeDir /scratch/twinkle1/STAR_index_dea \
	--outFileNamePrefix "${FILES["${SLURM_ARRAY_TASK_ID}"]/$INSTRING/$OUTSTRING}"_ \
	--readFilesIn "${FILES["${SLURM_ARRAY_TASK_ID}"]}" "${FILES["${SLURM_ARRAY_TASK_ID}"]/trimmed_1P.fq.gz/trimmed_2P.fq.gz}"

