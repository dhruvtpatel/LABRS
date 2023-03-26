#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=5

module load vcftools/0.1.17

vcftools --gzvcf ppmi_gwas_preQC_chr22_hwe.vcf.gz --max-missing 0.9 --recode --recode-INFO-all --stdout | gzip -c > ppmi_gwas_preQC_finished_chr22.vcf.gz

