#!/bin/bash

#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 2:00:00
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=100000

plink --bfile UK.MF.US.GER.FRcas_strand.lifted_remake_RawID.QC_10072019 --bmerge chr1-22.1KG3_phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.RAWID.filtered_2020.bed chr1-22.1KG3_phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.RAWID.filtered_2020.bim chr1-22.1KG3_phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.RAWID.filtered_2020.fam --make-bed --aec --allow-no-sex --out sample_remake_QC_merge_1kg3_test
