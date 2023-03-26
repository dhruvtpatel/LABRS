#!/bin/bash
#MSUB -A b1042
#MSUB -q genomics
#MSUB -l walltime=24:00:00
#MSUB -l nodes=1:ppn=6 

cd $PBS_O_WORKDIR

# Converting VCF to Plink

/projects/b1049/genetics_programs/plink --vcf /projects/b1049/bernabe/PDGC_QC_CR90.rawID.vcf --double-id --biallelic-only strict list --vcf-half-call missing --make-bed --out /projects/b1049/bernabe/PDGC_QC_CR90.rawID
