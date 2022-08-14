#!/bin/bash -l
#SBATCH -D /home/twinkle1/projects/seed_color_differential_expression_analysis/
#SBATCH -t 4:00:00
#SBATCH -J trimmomatic
#SBATCH -o logs/trimmomatic/mappingLog-%j.txt
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --mem=8gb
#SBATCH --array 0-25
#SBATCH --mail-user=twinkle1@smail.uni-koeln.de
#SBATCH --mail-type=ALL

module load trimmomatic/0.39

# there are a total of 26 samples, array goes from 0-25


### PREPARE INPUT

# create a new directory with all input files, name them according to the sample set

mkdir -p /scratch/twinkle1/rnaSeq_STAR_input

# make array of all read file from february
FEB=(/scratch/twinkle1/rnaSeq_seeds/Feb/AM*)

# make array of all read files from december
DEC=(/scratch/twinkle1/rnaSeq_seeds/Dec13/AM*)

# rename the files in order to prevent any issues
for index in "${!FEB[@]}";
do
	mv "${FEB[$index]}" "${FEB[$index]/AM/FEB}"
done

for index in "${!DEC[@]}";
do
       mv "${DEC[$index]}" "${DEC[$index]/AM/DEC}"
done

cp /scratch/twinkle1/rnaSeq_seeds/Feb/* /scratch/twinkle1/rnaSeq_STAR_input/
cp /scratch/twinkle1/rnaSeq_seeds/Dec13/* /scratch/twinkle1/rnaSeq_STAR_input/


### MAIN

# create array of read fastq files (R1 only):
SOURCE_DIR=/scratch/twinkle1/rnaSeq_STAR_input
FILES=("$SOURCE_DIR"/*R1.fastq.gz)

# run trimmomatic, use 6 threads, taking advantage of the baseout function to name output files
# use the sequencing adapters send by the sequencing center as custom fasta file
# forward and reverse read adapters are indicated in the fasta file by the /1 and /2 suffixes

java -jar $TRIMMOMATIC/trimmomatic.jar PE \
	-threads 6 \
	"${FILES["${SLURM_ARRAY_TASK_ID}"]}" \
	"${FILES["${SLURM_ARRAY_TASK_ID}"]/R1.fastq.gz/R2.fastq.gz}" \
	-baseout "${FILES["${SLURM_ARRAY_TASK_ID}"]/R1.fastq.gz/trimmed.fq.gz}" \
	ILLUMINACLIP:data/sequencing_adapters/custom_adapters.fa:2:30:10
