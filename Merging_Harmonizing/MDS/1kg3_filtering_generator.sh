for i in {1..22}
do
cat <<EOF > 1kg3.filtering.chr${i}.sh
#!/bin/bash
#MSUB -A b1042
#MSUB -q genomics
#MSUB -l walltime=24:00:00
#MSUB -l nodes=1:ppn=6 

cd \$PBS_O_WORKDIR

/projects/b1049/genetics_programs/vcftools/cpp/vcftools --gzvcf /projects/b1049/bernabe/1kg3/chr${i}.MultiSplit.IndelLeftNorm.RawID.vcf.gz --snps /projects/b1049/bernabe/PDGC_PostIdivQC_AllVars.variants.txt --recode --recode-INFO-all --out /projects/b1049/bernabe/1kg3/chr${i}.MultiSplit.IndelLeftNorm.RawID.filtered

EOF
done
