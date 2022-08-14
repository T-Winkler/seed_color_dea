#!/bin/bash


# automatically edit the names of the shortread seed fastq files to match the sample names
# read in the information from the sample names file and use it to replace the sample names
# remove unneccessary parts of the name as well as minus character symbols 
# this is the improved version of the script, that can copy the data from the raw data folder



##### copy all files from the first february sample to the new scratch directory

# raw data
RAW=/projects/ag-stetter/markus/amaranth_raw_data/rnaSeq_seeds_Dec2021/mstetter_MAS04_Februar4/

# /scratch/ directory with the rnaSeq reads
INDIR=/scratch/twinkle1/rnaSeq_seeds/Feb/
mkdir -p "$INDIR"

# copy all files from the raw data directory
cp "$RAW"A0* "$INDIR"



##### concatenate all file from the second february sample with the reads from the first sample

# raw data
RAWNINE=/projects/ag-stetter/markus/amaranth_raw_data/rnaSeq_seeds_Dec2021/mstetter_MAS04_Februar9/

# make array of all read file from feb9
FEBNINE=("$RAWNINE"A00*)

# make array of all read files from feb4
FEBFOUR=("$INDIR"A00*)

for ((i=0; i<=33; i++))
do
	cat "${FEBNINE[$i]}" >> "${FEBFOUR[$i]}"
done


##### continue by renaming sample according to the Feb4 sample name list
# sample list file
input=/projects/ag-stetter/markus/amaranth_raw_data/rnaSeq_seeds_Dec2021/mstetter_MAS04_Februar4/Sample_Names.tab

# read in from the sample names file
ID=($(awk '{print $1}' $input))
LABEL=($(awk '{print $2}' $input))



#### MAIN

# replace sample ID with label
for ((i=1; i<=17; i++))
do
	#ls "$INDIR"*"${ID[i]}"* | sed "s/A0068.*L001/${LABEL[i]}/"		
	# rename all files according to the sample list
	rename "${ID[i]}" "${LABEL[i]}" "$INDIR"/A006* 
done


# remove unncecessary part of the name:
for file in "$INDIR"A0068*
do
	mv "$file" ${file//A006850177_/} # remove unnecessary parts of the name
done

# remove minus characters
for file in "$INDIR"*fastq.gz
do
	mv "$file" ${file//-/_} # replace minus character with underscore
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

