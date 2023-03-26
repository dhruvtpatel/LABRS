#!/bin/bash
#MSUB -A b1042
#MSUB -q genomics
#MSUB -l walltime=24:00:00
#MSUB -l nodes=1:ppn=6 

cd $PBS_O_WORKDIR

/projects/b1049/genetics_programs/plink --bfile /projects/b1049/bernabe/PDGC_QC_CR90.rawID_wSEX_wPHENO_MIND05 --remove /projects/b1049/bernabe/PDGC_high_HetRateIndiv.txt --make-bed --out /projects/b1049/bernabe/PDGC_QC_CR90.rawID_wSEX_wPHENO_MIND05_NoHighHe

/projects/b1049/genetics_programs/plink --bfile /projects/b1049/bernabe/PDGC_QC_CR90.rawID_wSEX_wPHENO_MIND05_NoHighHe --remove /projects/b1049/bernabe/fail-IBD-QC.txt --make-bed --out /projects/b1049/bernabe/PDGC_QC_CR90.rawID_wSEX_wPHENO_MIND05_NoHighHe_NoIBD_AllVars
