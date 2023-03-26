#!/bin/bash
#MSUB -A b1042
#MSUB -q genomics
#MSUB -l walltime=24:00:00
#MSUB -l nodes=1:ppn=24 

cd $PBS_O_WORKDIR

module load perl

perl /projects/b1049/genetics_programs/annovar_2017/annovar/table_annovar.pl /projects/b1049/bernabe/PDGC_VarsQC_FullIndivQc.rawID_NoMaf.avinput /projects/b1049/genetics_programs/annovar_2017/annovar/humandb/ -protocol dbnsfp33a -operation f -build hg19 -nastring NA -thread 24
