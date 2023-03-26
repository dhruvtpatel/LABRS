#!/bin/bash
#MSUB -A b1042
#MSUB -q genomics
#MSUB -l walltime=24:00:00
#MSUB -l nodes=1:ppn=6 

cd $PBS_O_WORKDIR

/projects/b1049/genetics_programs/plink --bfile /projects/b1049/bernabe/PDGC_QC_CR90.rawID --update-sex /projects/b1049/bernabe/PDGC_SEX_INFO.txt --pheno /projects/b1049/bernabe/PDGC_PHENO_INFO.txt --mind 0.1 --make-bed --out /projects/b1049/bernabe/PDGC_QC_CR90.rawID_wSEX_wPHENO_MIND1

