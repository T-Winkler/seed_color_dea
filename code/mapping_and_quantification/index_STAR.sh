#!/bin/bash -l
#SBATCH -D /home/twinkle1/projects/seed_color_differential_expression_analysis/
#SBATCH -t 1:00:00
#SBATCH -J STAR
#SBATCH -o logs/STAR/dea/mappingLog-%j.txt
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=16gb
#SBATCH --job-name="index_STAR"

# IMPORTANT: it is necessary to start the script from the cheops1 partition because of the specified STAR version

module load star/2.7.8a

# Index the reference genome
# only run once per reference genome

# 8 threads, genome generation mode
# sjdbOverhang and genomeSAindexNbases settings specific for the amaranth reference assembly v2.1
# more specific settings: use the polished, softmasked reference assembly, also use the reorded version
# as SJDB file, use the new genome annotation

mkdir -p /scratch/twinkle1/STAR_index_dea/

STAR --runThreadN 8 \
	--runMode genomeGenerate \
	--genomeDir /scratch/twinkle1/STAR_index_dea \
	--sjdbOverhang 99 \
	--genomeSAindexNbases 13 \
	--sjdbGTFfeatureExon CDS \
	--genomeFastaFiles /home/twinkle1/projects/A_hypochondriacus_reannotation/polished_genome_annotation/assembly/Ahypochondriacus_2.2_polished.softmasked.fasta \
	--sjdbGTFfile /home/twinkle1/projects/A_hypochondriacus_reannotation/polished_genome_annotation/annotation/Ahypochondriacus_2.2_polished_corrected.gtf


