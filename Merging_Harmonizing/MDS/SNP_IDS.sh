#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --time=05:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=1G
#SBATCH --job-name=SNP_IDS
#SBATCH --output=Step5_3.1_output.log

#a regular bash comment
cd /projects/b1049/SajivH
module load vcftools/0.1.17
for i in {1..20}
do
vcftools --gzvcf ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.RAWID.vcf.gz --snps GWAS_Training_Step5.ID.txt --recode --recode-INFO-all --out ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.RAWID.filtered
done


