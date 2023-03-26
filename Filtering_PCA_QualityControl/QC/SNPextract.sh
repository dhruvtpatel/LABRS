#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=20

vcftoolsold=/projects/b1049/genetics_programs/vcftools_old/cpp/vcftools

$vcftoolsold --gzvcf ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.RAWID.vcf4.2.vcf.gz --snps ppmi_gwas.ID.txt --recode --recode-INFO-all --out ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.RAWID.filtered
