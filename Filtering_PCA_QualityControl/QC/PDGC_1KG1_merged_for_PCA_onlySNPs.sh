#!/bin/bash
#MSUB -A b1042
#MSUB -q genomics
#MSUB -l walltime=24:00:00
#MSUB -l nodes=1:ppn=6 

cd $PBS_O_WORKDIR


/projects/b1049/genetics_programs/plink --vcf /projects/b1049/bernabe/PDGC_PostIdivQC.chr1-22.1kg3.filtered.merged_OnlySNPs.recode.vcf --double-id --biallelic-only strict list --vcf-half-call missing --make-bed --out /projects/b1049/bernabe/PDGC_PostIdivQC.chr1-22.1kg3.filtered.merged_OnlySNPs

/projects/b1049/genetics_programs/plink --bfile /projects/b1049/bernabe/PDGC_PostIdivQC.chr1-22.1kg3.filtered.merged_OnlySNPs --maf 0.01 --allow-no-sex --indep-pairwise 50 5 0.2 --out /projects/b1049/bernabe/PDGC_PostIdivQC.chr1-22.1kg3.filtered.merged_OnlySNPs

/projects/b1049/genetics_programs/plink --bfile /projects/b1049/bernabe/PDGC_PostIdivQC.chr1-22.1kg3.filtered.merged_OnlySNPs --maf 0.01 --allow-no-sex --extract /projects/b1049/bernabe/PDGC_PostIdivQC.chr1-22.1kg3.filtered.merged_OnlySNPs.prune.in --make-bed --out /projects/b1049/bernabe/PDGC_PostIdivQC.chr1-22.1kg3.filtered.merged_OnlySNPs_MAF01_pruned
