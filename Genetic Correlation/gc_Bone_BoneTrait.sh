#!/bin/bash

# Activate the LDSC environment
cd /LDSC/ldsc
source activate ldsc

# Set the path to the LDSC software
LDSC_PATH=/LDSC/ldsc

# Set the paths to the two folders containing munge.sumstats.gz files
FOLDER1=/Data/sumstats/Bone
FOLDER2=/Data/sumstats/BoneTrait

# Create a new folder to store the output results
OUTPUT_FOLDER=/Data/Rg/Bone_BoneTrait
mkdir -p "$OUTPUT_FOLDER"

for file1 in "$FOLDER1"/*_munge.sumstats.gz; do
    for file2 in "$FOLDER2"/*_munge.sumstats.gz; do
    filename1=$(basename "$file1" _munge.sumstats.gz)
    filename2=$(basename "$file2" _munge.sumstats.gz)
    "$LDSC_PATH"/ldsc.py \
    --rg "$file1","$file2" \
    --ref-ld-chr "$LDSC_PATH"/eur_w_ld_chr/ \
    --w-ld-chr "$LDSC_PATH"/eur_w_ld_chr/ \
    --out "$OUTPUT_FOLDER"/"$filename1"_"$filename2"
    done
done
    
