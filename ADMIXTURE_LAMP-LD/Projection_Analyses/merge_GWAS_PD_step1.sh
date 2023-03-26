#!/bin/bash
#MSUB -A b1042
#MSUB -q genomics
#MSUB -l walltime=24:00:00
#MSUB -l nodes=1:ppn=2
    
cd $PBS_O_WORKDIR

perl /projects/b1049/bernabe/GWAS_GSD/bimAnnotationUpdate.pl /projects/b1049/kronos/gwas/UK.MF.US.GER.FRcas+strand.lifted.bim < /projects/b1049/bernabe/GWAS_GSD/refPanel.bim > UK.MF.US.GER.FRcas+strand.lifted.bim.info
