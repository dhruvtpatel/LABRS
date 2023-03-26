
#!/bin/bash
#MSUB -A b1042
#MSUB -q genomics
#MSUB -l walltime=24:00:00
#MSUB -l nodes=1:ppn=1

/projects/b1049/genetics_programs/vcftools/cpp/vcftools --vcf /projects/b1049/bernabe/PDGC_VarsQC_FullIndivQc.rawID_NoMaf_autosomes.NoMonomorfic.recode.vcf --chr 4 --recode --recode-INFO-all --out /projects/b1049/bernabe/VCF_per_chromosome/PDGC_VarsQC_FullIndivQc.rawID_NoMaf_autosomes.NoMonomorfic.chr4
