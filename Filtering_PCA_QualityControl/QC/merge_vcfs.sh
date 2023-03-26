#!/bin/bash
#MSUB -A b1042
#MSUB -q genomics
#MSUB -l walltime=24:00:00
#MSUB -l nodes=1:ppn=6 

cd $PBS_O_WORKDIR

/projects/b1049/genetics_programs/bcftools/bcftools merge -m id /projects/b1049/bernabe/PDGC_PostIdivQC_Shared_1kg3Vars.recode.vcf.gz /projects/b1049/bernabe/1kg3/chr1-22.1kg3.MultiSplit.IndelLeftNorm.RawID.filtered.vcf.gz -o /projects/b1049/bernabe/PDGC_PostIdivQC.chr1-22.1kg3.filtered.merged.vcf
