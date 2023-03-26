#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=5G
#SBATCH --job-name=makeVCF-check
#SBATCH --error=makeVCF.error
#SBATCH --output=makeVCF

module load plink/1.9
cd /projects/b1049/SajivH

for chnum in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23};
  do
  plink --bfile GWAS_Raw_Sorted_AllIndivQC-updated-chr$chnum --recode vcf --chr $chnum --out GWAS_Raw_Sorted_AllIndivQC$chnum 
done
