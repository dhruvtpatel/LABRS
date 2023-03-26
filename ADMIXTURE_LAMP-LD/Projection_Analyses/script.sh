#!/bin/sh

#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 4:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=28
#SBATCH --mail-user=dpatel4@imsa.edu
#SBATCH --mail-type=END
#SBATCH --mem=20000

plink2 --bfile ft_hwe --king-cutoff 0.177 --make-bed --out ft_ped

