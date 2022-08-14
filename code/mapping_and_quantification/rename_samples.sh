#!/bin/bash

# Rename the samples returned from sequencing.
# Copy the samples from the raw data directory to the /scratch/ partition, combine resequencing files with the original files and prepare adapter removal
# call the two helper scripts in order:

# first process reads from December
bash code/mapping_and_quantification/call_rename_Dec.sh

# process reads from February
bash code/mapping_and_quantification/call_rename_Feb9.sh
