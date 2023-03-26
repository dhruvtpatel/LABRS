#!/bin/bash
#MSUB -A b1042
#MSUB -q genomics
#MSUB -l walltime=24:00:00
#MSUB -l nodes=1:ppn=6 

cd $PBS_O_WORKDIR

module load perl


perl /projects/b1049/genetics_programs/annovar_2017/annovar/annotate_variation.pl /projects/b1049/bernabe/PDGC_VarsQC_FullIndivQc.rawID_NoMaf.avinput.hg19_avsnp150_filtered.hg19_gnomad_genome_filtered.hg19_ALL.sites.2015_08_filtered.hg19_exac03_filtered.hg19_esp6500siv2_all_filtered /projects/b1049/genetics_programs/annovar_2017/annovar/humandb/ -filter -build hg19 -dbtype cg69
