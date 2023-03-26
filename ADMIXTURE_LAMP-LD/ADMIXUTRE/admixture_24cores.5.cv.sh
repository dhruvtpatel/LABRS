
#!/bin/bash
#MSUB -A b1042
#MSUB -q genomics
#MSUB -l walltime=48:00:00
#MSUB -l nodes=1:ppn=24

cd $PBS_O_WORKDIR

/projects/b1049/genetics_programs/admixture_linux-1.3.0/admixture --cv UK.MF.US.GER.FRcas_strand.lifted_1KG3_merged_RawID_AllPheno_prunned.bed 5 -j24 | tee log5.out
