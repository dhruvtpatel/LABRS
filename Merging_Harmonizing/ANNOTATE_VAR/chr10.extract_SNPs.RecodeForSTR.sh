#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=16000

/projects/b1049/genetics_programs/plink_Jun_2019/plink --bfile chromosome_data/STR_analysis/chr10.UK.MF.US.GER.FRcas_strand.lifted_remake_RawID.QC_10072019 --recode vcf bgz --out chromosome_data/STR_analysis/chr10.UK.MF.US.GER.FRcas_strand.lifted_remake_RawID.QC_10072019.recoded --extract refPanel_rawID.SNPs.bim --real-ref-alleles --a2-allele refPanel_rawID.bim 5 2 
