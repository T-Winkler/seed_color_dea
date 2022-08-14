#!/bin/bash


# automatically edit the names of the shortread seed fastq files to match the sample names
# read in the information from the sample names file and use it to replace the sample names
# remove unneccessary parts of the name as well as minus character symbols 
# this is the improved version of the script, that can copy the data from the raw data folder


#### SETUP

# input sample names file
input=/projects/ag-stetter/markus/amaranth_raw_data/rnaSeq_seeds_Dec2021/mstetter_MAS03_November13/Sample_Names.tab

# raw data
RAW=/projects/ag-stetter/markus/amaranth_raw_data/rnaSeq_seeds_Dec2021/mstetter_MAS03_November13/

# read in from the sample names file
ID=($(awk '{print $1}' $input))
LABEL=($(awk '{print $2}' $input))


# ID[0] is the header
#echo ${ID[1]}
#echo ${LABEL[9]}

# /scratch/ directory with the rnaSeq reads
INDIR=/scratch/twinkle1/rnaSeq_seeds/Dec13/
mkdir -p "$INDIR"

# copy all files from the raw data directory
cp "$RAW"A0* "$INDIR"

#### MAIN

# did not find a way to automatically set the length of the array-1 as the number of iterations
for ((i=1; i<=9; i++))
do
	#ls "$INDIR"*"${ID[i]}"* | sed "s/A0068.*L001/${LABEL[i]}/"		
	# rename all files according to the sample list
	rename "${ID[i]}" "${LABEL[i]}" "$INDIR"/A006* 
done


# remove unncecessary part of the name:
for file in "$INDIR"A0068*
do
	mv "$file" ${file//A006850165_/} # remove unnecessary parts of the name
done


# remove minus characters
for file in "$INDIR"*fastq.gz
do
	mv "$file" ${file//-/_} # remove minus character with underscore
done


# remove unnecessary part of the name:
for file in "$INDIR"*fastq.gz
do
	mv "$file" ${file//_S*_R/_R} # remove unnecessary parts of the name
done

for file in "$INDIR"*fastq.gz
do
  	mv "$file" ${file//_001/} # remove unnecessary parts of the name
done

# added this part afterwards:
# add the resequencing files to the renamed files:
cat /projects/ag-stetter/markus/amaranth_raw_data/rnaSeq_seeds_Dec2021/additonal_reads_MAS03_November13/A006200198_159068_S67_L002_R1_001.fastq.gz >> "$INDIR"AM_00364_4_R1.fastq.gz
cat /projects/ag-stetter/markus/amaranth_raw_data/rnaSeq_seeds_Dec2021/additonal_reads_MAS03_November13/A006200198_159068_S67_L002_R2_001.fastq.gz >> "$INDIR"AM_00364_4_R2.fastq.gz
