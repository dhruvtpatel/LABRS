#!/bin/bash
#MSUB -A b1042
#MSUB -q genomics
#MSUB -l walltime=24:00:00
#MSUB -l nodes=1:ppn=6 

cd $PBS_O_WORKDIR

/projects/b1049/genetics_programs/vcftools/cpp/vcftools --vcf /projects/b1049/bernabe/step2.PDGC.vcf --minDP 8 --minGQ 20 --hwe 0.000001 --remove-filtered-all --recode --recode-INFO-all --out /projects/b1049/bernabe/PDGC_QC_noCR
