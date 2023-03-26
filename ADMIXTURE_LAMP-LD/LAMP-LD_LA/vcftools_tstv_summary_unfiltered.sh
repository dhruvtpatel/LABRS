#!/bin/bash
#MSUB -A b1042
#MSUB -q genomics
#MSUB -l walltime=24:00:00
#MSUB -l nodes=1:ppn=6 

cd $PBS_O_WORKDIR

# Ts/Tv summary for the exome-seq data with vcftools

/projects/b1049/genetics_programs/vcftools/cpp/vcftools --vcf /projects/b1049/bernabe/step2.PDGC.vcf --TsTv-summary --out /projects/b1049/bernabe/TsTv_unfiltered
