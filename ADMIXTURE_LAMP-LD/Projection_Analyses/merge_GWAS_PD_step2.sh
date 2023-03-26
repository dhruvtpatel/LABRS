#!/bin/bash
#MSUB -A b1042
#MSUB -q genomics
#MSUB -l walltime=24:00:00
#MSUB -l nodes=1:ppn=2
    
cd $PBS_O_WORKDIR

perl /projects/b1049/bernabe/GWAS_GSD/bimAnnotationUpdate.pl UK.MF.US.GER.FRcas+strand.lifted_toMerge.bim < refPanel_rawID.bim > UK.MF.US.GER.FRcas+strand.lifted_toMerge.bim.info
