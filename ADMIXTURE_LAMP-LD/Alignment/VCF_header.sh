#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=20

#extract out the header of the file to figure out individuals in annotated file
cat AMPPD_GBA.vcf | grep "^#[^#]" > AMPPD_header.txt
