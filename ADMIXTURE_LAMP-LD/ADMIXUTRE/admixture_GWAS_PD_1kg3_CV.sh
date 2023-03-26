#!/bin/bash
#MSUB -A b1042
#MSUB -q genomics
#MSUB -l walltime=24:00:00
#MSUB -l nodes=1:ppn=4

cd $PBS_O_WORKDIR

for K in 3 4 5 6 7 8 9 10
do

/projects/b1049/genetics_programs/admixture_linux-1.3.0/admixture --cv UK.MF.US.GER.FRcas_strand.lifted_1KG3_merged_RawID_AllPheno_prunned.bed $K | tee log${K}.out

done
