#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=20

#awk attempt 1
awk '$12=="GBA" {print}' /projects/b1049/jing/annovar/withfre.nospace.short.nozero.polymorphic.chr1-22.AMP.hg38_multianno.txt > AMP.hg38_GBA_mutations.txt

#awk attempt 2
