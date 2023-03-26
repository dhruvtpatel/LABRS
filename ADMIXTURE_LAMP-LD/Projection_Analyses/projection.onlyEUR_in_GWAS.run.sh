#!/bin/bash
#MSUB -A b1042
#MSUB -q genomics
#MSUB -l walltime=48:00:00
#MSUB -l nodes=1:ppn=22

cd $PBS_O_WORKDIR

/projects/b1049/genetics_programs/admixture_linux-1.3.0/admixture -j20 -P projection.onlyEUR_in_GWAS.bed 4
