
#!/bin/bash
#MSUB -A b1042
#MSUB -q genomics
#MSUB -l walltime=24:00:00
#MSUB -l nodes=1:ppn=1

cd $PBS_O_WORKDIR

/projects/b1049/genetics_programs/tabix-0.2.6/bgzip /projects/b1049/bernabe/VCF_per_chromosome/PDGC_VarsQC_FullIndivQc.rawID_NoMaf_autosomes.NoMonomorfic.chr16.recode.vcf
