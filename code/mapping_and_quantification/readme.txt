## RNAseq read mapping and quantification

All code used for processing, mapping, quality control and quantification of RNAseq reads from developing A. hypochondriacus seeds. 

### Script order:

- code/mapping_and_quantification/rename_samples.sh
renaming of samples, calls the two helper scripts call_rename_Dec.sh and call_rename_Feb9.sh

- code/mapping_and_quantification/adapter_trimming.sh
trimming of adapter sequences using trimmomatic

- code/mapping_and_quantification/index_STAR.sh
indexing of the polished reference genome to prepare read mapping using STAR

- code/mapping_and_quantification/run_STAR.sh
map RNA sequencing reads to the polished reference genome using STAR

- code/mapping_and_quantification/QC_and_counting.sh
perform quality control and quantify the number of read pairs for each gene in genome annotation v2.2
