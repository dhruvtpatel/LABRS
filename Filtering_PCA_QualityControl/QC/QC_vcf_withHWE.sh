#!/bin/bash
#MSUB -A b1042
#MSUB -q genomics
#MSUB -l walltime=24:00:00
#MSUB -l nodes=1:ppn=6 

cd $PBS_O_WORKDIR

# VCF quality controls

/projects/b1049/genetics_programs/vcftools/cpp/vcftools --vcf /projects/b1049/bernabe/step2.PDGC.vcf --remove-filtered-all --minGQ 20 --minDP 8 --max-missing 0.1 --hwe 0.000001 --recode --recode-INFO-all --out /projects/b1049/bernabe/QC_NoPass_GC20_DP8_CR10_HWE_PDGC.vcf
