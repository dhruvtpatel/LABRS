#!/bin/bash
#MSUB -A b1042
#MSUB -q genomics
#MSUB -l walltime=24:00:00
#MSUB -l nodes=1:ppn=6 

cd $PBS_O_WORKDIR

/projects/b1049/genetics_programs/plink --bfile /projects/b1049/bernabe/PDGC_QC_CR90.rawID_wSEX_wPHENO_MIND05 --maf 0.01 --indep-pairwise 50 5 0.2 --out /projects/b1049/bernabe/PDGC_QC_CR90.rawID_wSEX_wPHENO_MIND05_MAF01

/projects/b1049/genetics_programs/plink --bfile /projects/b1049/bernabe/PDGC_QC_CR90.rawID_wSEX_wPHENO_MIND05 --remove /projects/b1049/bernabe/PDGC_high_HetRateIndiv.txt --extract /projects/b1049/bernabe/PDGC_QC_CR90.rawID_wSEX_wPHENO_MIND05_MAF01.prune.in --make-bed --out /projects/b1049/bernabe/PDGC_QC_CR90.rawID_wSEX_wPHENO_MIND05_MAF01_NoHighHe_pruned
