#!/bin/bash -l
#SBATCH -D /home/twinkle1/projects/seed_color_differential_expression_analysis
#SBATCH -t 40:00:00
#SBATCH -J QC
#SBATCH -o logs/trimmomatic_qc/mappingLog-%j.txt
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --mem=8gb

# Use featurecounts to create the countmatrix from STAR mappings, also perform different quality control measures

# load necessary modules

source $CONDA_PREFIX/etc/profile.d/conda.sh
conda activate featurecounts

module load fastqc/0.11.9

GTFIN=/home/twinkle1/projects/A_hypochondriacus_reannotation/polished_genome_annotation/annotation/Ahypochondriacus_2.2_polished_corrected.gtf

# run featurecounts
# feature is CDS instead of exon, count on meta features (gene_id), use 10 threads, count paired-end fragments, only if both reads mapped
FCOUT=/scratch/twinkle1/featurecounts/

mkdir $FCOUT
featureCounts -a $GTFIN \
	-t CDS \
	-T 10 \
	-p \
	-B \
	-g gene_id \
	-o "$FCOUT"counts.txt \
	/scratch/twinkle1/rnaSeq_STAR_output/*1P.fq.gz_Aligned.sortedByCoord.out.bam


# Quality control
# Initialize the output directory
QCOUT=/scratch/twinkle1/rnaSeq_qc/

mkdir -p $QCOUT
# run quality control fastqc
fastqc -t 10 -o $QCOUT /scratch/twinkle1/rnaSeq_STAR_input/*P.fq.gz


module load samtools/1.13


# run quality control qualimap on the generated bam files
# run rnaseq mode for each file, qualimap takes as input a bam sorted by name
# -p = strand specific protocol, -pe = paired-end sequencing data, -s = file is sorted by name

# prepare input gtf file
GTFQM=/scratch/twinkle1/temp.gtf
sed 's/CDS/exon/' $GTFIN > $GTFQM

# run qualimap
for file in /scratch/twinkle1/rnaSeq_STAR_output/*out.bam
do
	OUTFILE=$(echo "${file/\/scratch\/twinkle1\/rnaSeq_STAR_output\//}")
	OUTFILEEDIT=$(echo "${OUTFILE/_trimmed_1P.fq.gz_Aligned.sortedByCoord.out.bam/}")
	samtools sort -n -T /scratch/twinkle1/ -@ 8 $file -O bam > "${file/trimmed_1P.fq.gz_Aligned.sortedByCoord.out.bam/name_sorted.bam}"
	rm -r "$QCOUT""$OUTFILEEDIT"
	mkdir "$QCOUT""$OUTFILEEDIT"
	qualimap rnaseq -bam "${file/trimmed_1P.fq.gz_Aligned.sortedByCoord.out.bam/name_sorted.bam}" \
		-outdir "$QCOUT""$OUTFILEEDIT" \
		-outfile $OUTFILEEDIT \
		-gtf $GTFQM \
		-p strand-specific-reverse \
		-pe \
		-s \
		--java-mem-size=4G
done

rm $GTFQM

# run multiqc to combine the results from fastqc and qualimap into a single report
MULTIQCOUT=/scratch/twinkle1/multiqc_seed
mkdir -p $MULTIQCOUT

multiqc -o $MULTIQCOUT $QCOUT

# copy files to different partition
mkdir -p data/featurecounts/
cp "$FCOUT"* data/featurecounts/


