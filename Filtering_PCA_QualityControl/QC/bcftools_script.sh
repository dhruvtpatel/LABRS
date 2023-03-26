#!/bin/bash
#MSUB -A b1042
#MSUB -q genomics
#MSUB -l walltime=24:00:00
#MSUB -l nodes=1:ppn=6 

cd $PBS_O_WORKDIR

# left normalizartion of indels with BCFtools

/projects/b1049/genetics_programs/bcftools/bcftools norm -m-both -o /projects/b1049/bernabe/step1.PDGC.vcf /projects/b1049/kronos/exome/VCF/pd_jan_20_2015_FILTERED.vcf

/projects/b1049/genetics_programs/bcftools/bcftools norm -f /projects/b1049/genetics_refs/fasta/hg19.fa -o /projects/b1049/bernabe/step2.PDGC.vcf /projects/b1049/bernabe/step1.PDGC.vcf
