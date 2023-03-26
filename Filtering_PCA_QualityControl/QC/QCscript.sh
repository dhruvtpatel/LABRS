#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=20

vcftoolsold=/projects/b1049/genetics_programs/vcftools_old/cpp/vcftools
plink=/projects/b1049/genetics_programs/plink_feb2018/plink

#PreQC for each chromosome (ppmi)
#$vcftoolsold --gzvcf /projects/b1049/ppmi/WGS_2019/ppmi.feb.2018.chr1.vqsr.vcf.gz --remove-filtered-all --minDP 5.0 --minGQ 20.0 --recode --recode-INFO-all --stdout | gzip -c > ppmi_gwas_preQC_chr1.vcf.gz

#PreQC all data
#$vcftoolsold --vcf IPDGC/pd_jan_20_2015_FILTERED.vcf --remove-filtered-all --minDP 5.0 --minGQ 20.0 --recode --recode-INFO-all --stdout | gzip -c > IPDGC_preQC.vcf.gz

#Hardy Weinberg equilibrium QC
#$vcftoolsold --gzvcf IPDGC_preQC.vcf.gz --hwe 0.000001 --recode --recode-INFO-all --stdout | gzip -c > IPDGC_preQC_hwe.vcf.gz

#Max-missing QC
#$vcftoolsold --gzvcf IPDGC_preQC_hwe.vcf.gz --max-missing 0.8 --recode --recode-INFO-all --stdout | gzip -c > IPDGC_preQC_finished.vcf.gz

#Combine VCFs (ppmi)
#export PERL5LIB=/projects/b1049/genetics_programs/vcftools_old/perl

#/projects/b1049/genetics_programs/vcftools_old/perl/vcf-concat ppmi_gwas_preQC_finished_chr1.vcf.gz ppmi_gwas_preQC_finished_chr2.vcf.gz ppmi_gwas_preQC_finished_chr3.vcf.gz ppmi_gwas_preQC_finished_chr4.vcf.gz ppmi_gwas_preQC_finished_chr5.vcf.gz ppmi_gwas_preQC_finished_chr6.vcf.gz ppmi_gwas_preQC_finished_chr7.vcf.gz ppmi_gwas_preQC_finished_chr8.vcf.gz ppmi_gwas_preQC_finished_chr9.vcf.gz ppmi_gwas_preQC_finished_chr10.vcf.gz ppmi_gwas_preQC_finished_chr11.vcf.gz ppmi_gwas_preQC_finished_chr12.vcf.gz ppmi_gwas_preQC_finished_chr13.vcf.gz ppmi_gwas_preQC_finished_chr14.vcf.gz ppmi_gwas_preQC_finished_chr15.vcf.gz ppmi_gwas_preQC_finished_chr16.vcf.gz ppmi_gwas_preQC_finished_chr17.vcf.gz ppmi_gwas_preQC_finished_chr18.vcf.gz ppmi_gwas_preQC_finished_chr19.vcf.gz ppmi_gwas_preQC_finished_chr20.vcf.gz ppmi_gwas_preQC_finished_chr21.vcf.gz ppmi_gwas_preQC_finished_chr22.vcf.gz ppmi_gwas_preQC_finished_chrX.vcf.gz | gzip -c > ppmi_gwas_preQC_finished.vcf.gz

#Change the SNP IDs (Run in shell)
#zcat IPDGC_preQC_finished.vcf.gz | grep '^#' > IPDGC_preQC_finished_header.vcf
#zcat IPDGC_preQC_finished.vcf.gz | grep -v '^#' | awk 'BEGIN {OFS="\t"}{$3=$1"_"$2"_"$4"_"$5} {print}' | cat IPDGC_preQC_finished_header.vcf - | gzip -c > IPDGC_preQC_finished_RawID.vcf.gz
#rm IPDGC_preQC_finished_header.vcf

#From here, swap ppmi_gwas_preQC_finished.vcf.gz with the RawID file "ppmi_gwas_preQC_finished_RawID.vcf.gz"

#TsTv for whole genome
#$vcftoolsold --gzvcf IPDGC_preQC_finished.vcf.gz --TsTv-summary

#Transfering to plink format
#$plink --vcf IPDGC_preQC_finished_RawID.vcf.gz --double-id --biallelic-only strict list --vcf-half-call missing --make-bed --out IPDGC_WES

#Update sex (Do missingness and heterozygosity first)
#$plink -noweb -bfile IPDGC_WES --update-sex IPDGC_SEX_INFO.txt --make-bed --out IPDGC_WES_sex

#Extract pseudoautosomal region (ppmi)
#plink --bfile IPDGC_WES_sex --split-x 2781479 155701383 'no-fail' --make-bed --out IPDGC_WES_sex_par

#Sex check
#$plink --bfile IPDGC_WES_sex_par --set-hh-missing --mind 0.1 --maf 0.05 --check-sex 0.5 0.5 --out IPDGC_WES_sexcheck --noweb
#grep PROBLEM ppmi_gwas_sexcheck.sexcheck | cut -f3 > fail-sexcheckqc.txt

#Heterozygosity and missingness
#$plink --bfile IPDGC_WES --missing --out IPDGC_WES_missing
#$plink --bfile IPDGC_WES --het --out IPDGC_WES_het

#Create graph with heterozygosity rate and rate of missing SNPs per individual (shell)
#module load R/3.6.0

#Install geneplotter through Bioconductor
#open R by typing "R"
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")

#BiocManager::install("geneplotter")

#Run script for graph (shell)
#Rscript imiss-vs-het.brent.Rscript

#Export file into computer to view and use thresholds to get rid of outliers using command prompt
#pscp bml7208@quest.it.northwestern.edu:/projects/b1049/brent/QC/IPDGC_WES_imiss-vs-het.pdf C:\Users\sjl2148\Desktop

#Create file with Family and Individual IDs of outliers (Go back to sex check)
#cat > fail-misshetqc.txt
#Family ID Individual ID
#^C

#Duplicated or related individuals
#LD prune data
#$plink --bfile IPDGC_WES --maf 0.01 --indep-pairwise 50 5 0.2 --out IPDGC_WES_LDprune

#Generate pairwise IBS data based on LD-pruned data
#$plink --bfile IPDGC_WES --extract IPDGC_WES_LDprune.prune.in --genome --out IPDGC_WES_dup_ind

#Rename file in order to run perl script
#cp IPDGC_WES_dup_ind.genome IPDGC_WES_missing.genome

#Identify individuals with IBD > 0.185 and put them into a file
#perl run-IBD-QC.pl IPDGC_WES_missing

#PCA
#Download 1k genomes project files
#wget (link address for each file on http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/)

#Create bim file of 1k genomes project sites
#zcat ../1000gen/ALL.wgs.shapeit2_integrated_v1a.GRCh38.20181129.sites.vcf.gz | grep -v "#" | awk '{print $1 "\t" $1"_"$2"_"$4"_"$5 "\t" "0" "\t" $2 "\t" $4 "\t" $5}' > 1kgGRCh38.bim

#Using bimAnnotationUpdate.pl, pull common variants between your own  data and 1000 genomes data
#perl ../../scripts/perl_scripts/bimAnnotationUpdate.pl IPDGC_WES.bim < /projects/b1049/morgan/1kg3p3.bim > IPDGC_WES.bim.info

#Examine SNPs carefully
#cut -f 7,8 IPDGC_WES.bim.info | sort | uniq -c

#Keep good SNPs
awk '$8>0{print $2}' IPDGC_WES.bim.info > IPDGC_WES.snps 

#$plink --bfile ppmi_gwas --extract ppmi_gwas.snps --make-bed --out ppmi_gwas.goodsnps

#Save file just to have a backup
#cp ppmi_gwas.goodsnps.bim ppmi_gwas.goodsnps.bim.orig

#Verify results
#perl ../scripts/perl_scripts/bimAnnotationUpdate.pl ppmi_gwas.goodsnps.bim.orig < 1kgGRCh38.bim > ppmi_gwas.goodsnps.bim.orig.info

#Generate new bim file with corrected SNP ID and allele coding/orientation
#awk '{if($8==2){a=$14;$14=$15;$15=a};print $10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15}' \ppmi_gwas.goodsnps.bim.orig.info > ppmi_gwas.goodsnps.bim

#Change SNP IDs for 1kg vcfs
#for i in {1..22}
#do
#zcat ../1000gen/ALL.chr${i}.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz | grep '^#' > ALL.chr${i}.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.header

#zcat ../1000gen/ALL.chr${i}.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz | grep -v '^#' | awk 'BEGIN {OFS="\t"}{$3=$1"_"$2"_"$4"_"$5} {print}' | cat ALL.chr${i}.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.header - | gzip -c >  ALL.chr${i}.shapeit2_integrated_v1a.GRCh38.20181129.phased.RAWID.vcf.gz

#rm ALL.chr${i}.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.header
#done

#Take the files and change the file version string from 4.3 to 4.2
#for i in {1..22}
#do
#zcat ALL.chr${i}.shapeit2_integrated_v1a.GRCh38.20181129.phased.RAWID.vcf.gz | sed 's/##fileformat=VCFv4.3/##fileformat=VCFv4.2/' | gzip -c > ALL.chr${i}.shapeit2_integrated_v1a.GRCh38.20181129.phased.RAWID.vcf4.2.vcf.gz
#done

#Take SNP IDs from ppmi and extract from the 1kg files
#for i in {1..22}
#do
#$vcftoolsold --gzvcf ALL.chr${i}.shapeit2_integrated_v1a.GRCh38.20181129.phased.RAWID.vcf4.2.vcf.gz --snps ppmi_gwas.ID.txt --recode --recode-INFO-all --out ALL.chr${i}.shapeit2_integrated_v1a.GRCh38.20181129.phased.RAWID.filtered
#done 

#Concatenate all chromosomes from the resulting VCFs
#for i in {1..22}
#do
#echo ALL.chr${i}.shapeit2_integrated_v1a.GRCh38.20181129.phased.RAWID.filtered.recode.vcf >> vcf_input.txt
#done

#Take the text file and convert to VCF with VCF-concat
#export PATH=$PATH:/projects/b1049/genetics_programs/vcftools_old/perl/
#export PERL5LIB=/projects/b1049/genetics_programs/vcftools_old/perl/
#vcf-concat `cat vcf_input.txt` | gzip -c > chr1-22.shapeit2_integrated_v1a.GRCh38.20181129.phased.RAWID.filtered.recode.vcf.gz

#Take VCF and convert to plink
#$plink --vcf chr1-22.shapeit2_integrated_v1a.GRCh38.20181129.phased.RAWID.filtered.recode.vcf.gz --double-id --biallelic-only strict list --vcf-half-call missing --make-bed --out chr1-22.shapeit2_integrated_v1a.GRCh38.20181129.phased.RAWID.filtered

#Merge 1kg file with ppmi
#$plink --bfile ppmi_gwas.goodsnps --bmerge chr1-22.shapeit2_integrated_v1a.GRCh38.20181129.phased.RAWID.filtered --make-bed --out file.1kgGRCh38.merged

#Prune data, but this time with 1kg data
#$plink --bfile file.1kgGRCh38.merged --maf 0.01 --indep-pairwise 50 5 0.2 --out ppmi_gwas.1kgGRCh38.merged
#$plink --bfile file.1kgGRCh38.merged --extract ppmi_gwas.1kgGRCh38.merged.prune.in --make-bed --out ppmi_gwas.1kgGRCh38.merged.pruned

#Create pca-populations.txt
#cat > pca-populations.txt
#3
#4
#5
#6
#7

#Change the superpopulation information in the 6th column of ppmi_gwas.1kgGRCh38.merged.fam

#Create .pedind and .pedsnp files
#cp ppmi_gwas.1kgGRCh38.merged.pruned.populations.fam ppmi_gwas.1kgGRCh38.merged.pruned.pedind
#cp ppmi_gwas.1kgGRCh38.merged.pruned.bim ppmi_gwas.1kgGRCh38.merged.pruned.pedsnp

#Run PCA through EIGENSOFT software (10eigenvectors with 10 cores)
#export PATH=$PATH:/projects/b1049/genetics_programs/EIG-6.1.4/bin/
#smartpca.perl -i ppmi_gwas.1kgGRCh38.merged.pruned.bed -a ppmi_gwas.1kgGRCh38.merged.pruned.pedsnp -b ppmi_gwas.1kgGRCh38.merged.pruned.pedind -o ppmi_gwas.1kgGRCh38.merged.pruned.pca -p ppmi_gwas.1kgGRCh38.merged.pruned.plot -e ppmi_gwas.1kgGRCh38.merged.pruned.eval -l ppmi_gwas.1kgGRCh38.merged.pruned.log -k 10 -t 10 -w pca-populations.txt

#If PCA fails with IDs that are too long, do --remove-indels in ppmi vcf file and repeat steps

#Create a file with all individuals failing QC. No duplicates are allowed
#cat > fail_qc_inds.txt
#Family ID, Individual ID
#^C

#Remove QC-failed individuals (There cannot be repeated individuals)
#$plink --bfile ppmi_gwas --remove fail_qc_inds.txt --maf 0.01 --geno 0.05 --hwe 1E-6 --mind 0.05 --make-bed --out clean_ppmi_gwas

#Add in the superpopulation information in the final fam file and you're finished
