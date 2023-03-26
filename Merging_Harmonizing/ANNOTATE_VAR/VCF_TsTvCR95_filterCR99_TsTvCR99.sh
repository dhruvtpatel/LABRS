#!/bin/bash
#MSUB -A b1042
#MSUB -q genomics
#MSUB -l walltime=24:00:00
#MSUB -l nodes=1:ppn=6 

cd $PBS_O_WORKDIR

# VCF Ts/Tv summmary

/projects/b1049/genetics_programs/vcftools/cpp/vcftools --vcf /projects/b1049/bernabe/QC_NoPass_GC20_DP8_CR10_HWE_PDGC.V3.vcf.recode.vcf.recode.vcf --TsTv-summary --out /projects/b1049/bernabe/TsTv_filtered_CR95 &

/projects/b1049/genetics_programs/vcftools/cpp/vcftools --vcf /projects/b1049/bernabe/QC_NoPass_GC20_DP8_CR10_HWE_PDGC.vcf.recode.vcf --max-missing 0.01 --recode --recode-INFO-all --out /projects/b1049/bernabe/QC_NoPass_GC20_DP8_CR99_HWE_PDGC.V4 &

/projects/b1049/genetics_programs/vcftools/cpp/vcftools --vcf /projects/b1049/bernabe/QC_NoPass_GC20_DP8_CR99_HWE_PDGC.V4.recode.vcf --TsTv-summary --out /projects/b1049/bernabe/TsTv_filtered_CR99
