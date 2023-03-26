#!/bin/bash
#MSUB -A b1042
#MSUB -l naccesspolicy=singlenode
#MSUB -l walltime=48:00:00 
#MSUB -q genomicsburst
cd $PBS_O_WORKDIR

#/projects/b1049/genetics_programs/bcftools/bcftools merge -m both /projects/b1049/Niccolo_NGS/genomes/bams/SS4009013/9013_ExpHunter_12032018.vcf  /projects/b1049/Niccolo_NGS/genomes/bams/SS4009014/9014_ExpHunter_12032018.vcf  /projects/b1049/Niccolo_NGS/genomes/bams/SS4009015/9015_ExpHunter_12032018.vcf   /projects/b1049/Niccolo_NGS/genomes/bams/SS4009016/9016_ExpHunter_12032018.vcf  /projects/b1049/Niccolo_NGS/genomes/bams/SS4009017/9017_ExpHunter_12032018.vcf  /projects/b1049/Niccolo_NGS/genomes/bams/SS4009020/9020_ExpHunter_12032018.vcf  /projects/b1049/Niccolo_NGS/genomes/bams/SS4009021/9021_ExpHunter_12032018.vcf  /projects/b1049/Niccolo_NGS/genomes/bams/SS4009022/9022_ExpHunter_12032018.vcf  /projects/b1049/Niccolo_NGS/genomes/bams/SS4009023/9023_ExpHunter_12032018.vcf  /projects/b1049/Niccolo_NGS/genomes/bams/SS4009030/9030_ExpHunter_12032018.vcf -o /projects/b1049/Niccolo_NGS/genomes/bams/nm_gemomes_xlinkedregion_exphunter_19032018.vcf

/projects/b1049/genetics_programs/vcftools/bin/vcf-sort -c NM_SS400_EXOMEonly_02102017_EXOME_SNP.recal.snps.indel_VQSR_hardFiltered.vcf_alleleReduction.vcf.recode.vcf > NM_SS400_EXOMEonly_sorted.vcf

#/projects/b1049/genetics_programs/vcftools/bin/vcftools --vcf /projects/b1042/LubbeLab/temp/VIPs_genome_13032018.vcf --plink --remove-filtered-all --remove-filtered-geno LowDepth --remove-filtered-geno LowGQ --remove-filtered-geno LowGQmean --min-alleles 2 --max-alleles 2 --minQ 30 --minDP 5 --out /projects/b1042/LubbeLab/temp/VIPs_genome_13032018_qc_15032018

#common vars only, ld prune, pi het and ibd