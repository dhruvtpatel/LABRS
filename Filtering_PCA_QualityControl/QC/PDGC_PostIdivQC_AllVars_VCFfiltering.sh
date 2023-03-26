#!/bin/bash
#MSUB -A b1042
#MSUB -q genomics
#MSUB -l walltime=24:00:00
#MSUB -l nodes=1:ppn=6 

cd $PBS_O_WORKDIR

# Filter CR90 VCF with good-quality individuals

/projects/b1049/genetics_programs/vcftools/cpp/vcftools --vcf /projects/b1049/bernabe/PDGC_QC_CR90.rawID.vcf --keep /projects/b1049/bernabe/PDGC_PstQC_indivList.txt --recode --recode-INFO-all --out /projects/b1049/bernabe/PDGC_PostIdivQC_AllVars
