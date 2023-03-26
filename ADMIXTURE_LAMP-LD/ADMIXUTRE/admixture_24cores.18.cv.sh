
#!/bin/bash
#MSUB -A b1042
#MSUB -q genomicslong
#MSUB -l walltime=240:00:00
#MSUB -l nodes=1:ppn=4

cd $PBS_O_WORKDIR

/projects/b1049/genetics_programs/admixture_linux-1.3.0/admixture --cv UK.MF.US.GER.FRcas_strand.lifted_1KG3_merged_RawID_AllPheno_prunned.bed 18 -j4 | tee log18.out
