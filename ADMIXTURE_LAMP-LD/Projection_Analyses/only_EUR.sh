#!/bin/bash
#MSUB -A b1042
#MSUB -q genomics
#MSUB -l walltime=48:00:00
#MSUB -l nodes=1:ppn=22

cd $PBS_O_WORKDIR

/projects/b1049/genetics_programs/admixture_linux-1.3.0/admixture -j20 UK.MF.US.GER.FRcas_strand.lifted_1KG3_merged_RawID_AllPheno_prunned_ONLY_EUR.bed 5
