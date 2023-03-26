#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=40G
#SBATCH --job-name=HRC-check
#SBATCH --error=HRC.error
#SBATCH --output=HRC

cd /projects/b1049/SajivH
module load perl/5.26

perl HRC-1000G-check-bim.pl -b GWAS_Raw_Sorted_AllIndivQC.bim -f GWAS_Raw_Sorted_AllIndivQC.frq -r HRC.r1-1.GRCh37.wgs.mac5.sites.tab -h

