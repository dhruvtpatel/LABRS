#!/bin/bash
#December 2016 - modified to work on b1042
#MSUB -A b1042
#MSUB -l naccesspolicy=singlenode
#MSUB -l walltime=336:00:00 
#MSUB -q genomicsburst
cd $PBS_O_WORKDIR
##first need to load the correct version of java (v1.8)
module load java

#Program to (i) genotype gVCF files in batches, and (ii) filter resultant combined VCF

#Step_by_step:
#1)Combine gVCF files into a joint gVCF file, 
#2)Genotype joint gVCF file
#3)VQSR on batch (set variant filter level)

#STEP 1 - you need to define your directories:
iDirectory=/projects/b1049/genomes/ON_5179/GVCFs #folder where our gVCFs are located #put all gVCFs in one folder
oDirectory=/projects/b1049/genomes/ON_5179/VCFs #folder where our VCF files will end up
temp=/projects/b1049/temp

##STEP 2 - set your filter level here #I generally use 99.9
VQSRFILTERLEVEL="99.9"

#STEP 3 - provide a list of your samples: (make sure you run with 10 more samples otherwise will cause inbreeding coeffecient error)
#works best with 50-100 samples! go over and you may run into computation problems #CAN I GET THIS TO WORK ON MORE????

mygVCFs="PD178
SS6005179
SJ_2001
FC08701
BM_2000
PD23100" 

#STEP 4 - give your analysis batch a name:
BATCH="ON_LL_BM_DN_exome_10092017" 

#These are our resources
GATK=/projects/b1049/genetics_programs/gatk/GenomeAnalysisTK.jar
GENOMEREF=/projects/b1049/genetics_refs/fasta/human_g1k_v37.fasta		#/projects/b1049/genetics_refs/fasta/hg19.fa
INTERVALS="/projects/b1049/genetics_refs/INTERVALS.bed"
ANNOVAR=/projects/b1049/genetics_programs/annovar/table_annovar.pl
ANNOdbs="refGene,genomicSuperDups,esp6500siv2_ea,dbnsfp30a,cg69,exac03,1000g2015aug_all,snp132,clinvar_20160302"
HARDfilterArgs=" --genotypeFilterExpression \"DP < 5\" --genotypeFilterName \"LowDepth\" --genotypeFilterExpression \"GQ < 20.0 && GQ > 0.0\" --genotypeFilterName \"LowGQ\" --filterExpression \"GQ_MEAN < 35.0 && GQ_MEAN > 0.0\" --filterName \"LowGQmean\""
VQSRargsSNP=" -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR --maxGaussians 5 "
VQSRresourcesSNP=" -resource:hapmap,known=false,training=true,truth=true,prior=15.0 /projects/b1049/genetics_refs/hapmap_3.3.hg19.sites.vcf -resource:omni,known=false,training=true,truth=true,prior=12.0 /projects/b1049/genetics_refs/1000G_omni2.5.hg19.sites_chr.vcf -resource:1000G,known=false,training=true,truth=false,prior=10.0 /projects/b1049/genetics_refs/1000G_phase1.snps.high_confidence.hg19.sites_chr.vcf -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /projects/b1049/genetics_refs/dbsnp_138.hg19_chr.vcf"
VQSRargsINDEL=" -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 --maxGaussians 4 -an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum "
VQSRresourcesINDEL="-resource:mills,known=false,training=true,truth=true,prior=12.0 /projects/b1049/genetics_refs/Mills_and_1000G_gold_standard.indels.hg19_modified_chr.vcf -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /projects/b1049/genetics_refs/dbsnp_138.hg19_chr.vcf"
chIDs="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT"
VCFtools=/projects/b1049/genetics_programs/vcftools/bin/vcftools

#############################################################

#STEP 5 - let's remove any old input lists:
#rm $temp/${BATCH}_gVCF_inputs.txt

#STEP 6 - let's make a list of our .gVCF file inputs and paths:
#for nID in $mygVCFs; do
#	echo "--variant $iDirectory/${nID}/${nID}.raw.snps.indels.EXOME.g.vcf " >> $temp/${BATCH}_gVCF_inputs.txt
#done

#STEP 7 - CombineGVCFs [CombineGVCFS] #include -L if you want to focus on exome only 
java -jar $GATK -R $GENOMEREF -T CombineGVCFs  `cat $temp/${BATCH}_gVCF_inputs.txt` -o $temp/${BATCH}.g.vcf -L $INTERVALS 

#STEP 8 - Joint Genotyping [GenotypeGVCFS] #note we really need multithreading here??? -nct 24 does not work!!
#for chr in $chIDs; do
#	java -jar $GATK -R $GENOMEREF -L $INTERVALS -T GenotypeGVCFs --variant $temp/${BATCH}.g.vcf -o $temp/${BATCH}_EXOME.vcf #--includeNonVariantSites
#done

#STEP 9 - now we have to join all of our split .vcf files back together.
#java -cp $GATK org.broadinstitute.gatk.tools.CatVariants -R $GENOMEREF -V $temp/${BATCH}_1.vcf -V $temp/${BATCH}_2.vcf -V $temp/${BATCH}_3.vcf -V $temp/${BATCH}_4.vcf -V $temp/${BATCH}_5.vcf -V $temp/${BATCH}_6.vcf -V $temp/${BATCH}_7.vcf -V $temp/${BATCH}_8.vcf -V $temp/${BATCH}_9.vcf -V $temp/${BATCH}_10.vcf -V $temp/${BATCH}_11.vcf -V $temp/${BATCH}_12.vcf -V $temp/${BATCH}_13.vcf -V $temp/${BATCH}_14.vcf -V $temp/${BATCH}_15.vcf -V $temp/${BATCH}_16.vcf -V $temp/${BATCH}_17.vcf -V $temp/${BATCH}_18.vcf -V $temp/${BATCH}_19.vcf -V $temp/${BATCH}_20.vcf -V $temp/${BATCH}_21.vcf -V $temp/${BATCH}_22.vcf -V $temp/${BATCH}_X.vcf -V $temp/${BATCH}_Y.vcf -V $temp/${BATCH}_MT.vcf -out $temp/${BATCH}.vcf 

java -jar $GATK -R $GENOMEREF -T SelectVariants -V $temp/${BATCH}_EXOME.vcf -selectType SNP -o $temp/${BATCH}_snps.vcf -L $INTERVALS

#STEP 10 - we must now do some filtering/recalibration of the calls #first, hard filter DP and GQ #remember to include -L from here on if you want to focus on exome only
java -jar $GATK -R $GENOMEREF -T VariantFiltration -V $temp/${BATCH}_snps.vcf --genotypeFilterExpression "DP < 5" --genotypeFilterName "LowDepth" --genotypeFilterExpression "GQ < 20.0 && GQ > 0.0" --genotypeFilterName "LowGQ" --genotypeFilterExpression "FS > 60.0" --genotypeFilterName "SB_snp" --genotypeFilterExpression "QD < 2.0" --genotypeFilterName "LowQD" --genotypeFilterExpression "MQ < 40.0" --genotypeFilterName "LowMQ" --genotypeFilterExpression "MQRankSum < -12.5" --genotypeFilterName "LowMQRankSum" --genotypeFilterExpression "ReadPosRankSum < -8.0" --genotypeFilterName "LowReadPosRankSum" -o $temp/${BATCH}_snps_filtered.vcf -L $INTERVALS #--filterExpression "GQ_MEAN < 35.0 && GQ_MEAN > 0.0" --filterName "LowGQmean"

#$JAVA -jar $GATK -R $GENOMEREF -L $INTERVALS -T VariantFiltration --genotypeFilterExpression "DP < 8" --genotypeFilterName "LowDepth" -V $temp/${BATCH}.vcf  -o $OUT/${BATCH}_HF1.vcf
#$JAVA -jar $GATK -R $GENOMEREF -L $INTERVALS -T VariantFiltration --genotypeFilterExpression "GQ < 20.0 && GQ > 0.0" --genotypeFilterName "LowGQ"  -V $OUT/${BATCH}_HF1.vcf   -o $OUT/${BATCH}_HF2.vcf
#$JAVA -jar $GATK -R $GENOMEREF -L $INTERVALS -T VariantFiltration --filterExpression "GQ_MEAN < 35.0 && GQ_MEAN > 0.0" --filterName "LowGQmean" -V $OUT/${BATCH}_HF2.vcf -o $OUT/${BATCH}_HF3.vcf

java -jar $GATK -R $GENOMEREF -T SelectVariants -V $temp/${BATCH}_EXOME.vcf -selectType INDEL -o $temp/${BATCH}_indels.vcf -L $INTERVALS

java -jar $GATK -R $GENOMEREF -T VariantFiltration -V $temp/${BATCH}_EXOME.vcf --genotypeFilterExpression "FS > 200.0" --genotypeFilterName "SB_indel" --genotypeFilterExpression "QD < 2.0" --genotypeFilterName "LowQD" --genotypeFilterExpression "ReadPosRankSum < -20.0" --genotypeFilterName "LowReadPosRankSum" -o $temp/${BATCH}_indels_filtered.vcf -L $INTERVALS

java -cp $GATK org.broadinstitute.gatk.tools.CatVariants -R $GENOMEREF -V $temp/${BATCH}_snps_filtered.vcf -V $temp/${BATCH}_indels_filtered.vcf -out $temp/${BATCH}_hardFiltered.vcf 

#STEP 11 - let's flag filtered variants, rather than setting them to No Call
java -jar $GATK -R $GENOMEREF -T SelectVariants --variant $temp/${BATCH}_snps_filtered.vcf --variant $temp/${BATCH}_indels_filtered.vcf -o $temp/${BATCH}_hardFiltered.vcf -L $INTERVALS

####VQSR runs best on many files - if using a few, just use the hard filtering step

#STEP 12 - VARIANT QUALITY SCORE RECALIBRATION - check tranch filter levels
#VQSR Step 1 [VariantRecalibrator] - SNPs
java -jar $GATK -R $GENOMEREF -T VariantRecalibrator -input $temp/${BATCH}_hardFiltered.vcf -mode SNP -recalFile $temp/${BATCH}_SNP.recal -tranchesFile $temp/${BATCH}_SNP.tranches -rscriptFile $temp/${BATCH}.plots.R $VQSRargsSNP $VQSRresourcesSNP -L $INTERVALS 
#VQSR Step 2 [ApplyRecalibration] - SNPs
java -jar $GATK -R $GENOMEREF -T ApplyRecalibration -mode SNP --ts_filter_level $VQSRFILTERLEVEL -input $temp/${BATCH}_hardFiltered.vcf -tranchesFile $temp/${BATCH}_SNP.tranches -recalFile $temp/${BATCH}_SNP.recal -o $temp/${BATCH}_hardFiltered_SNP.recal.snps.vcf -L $INTERVALS
#VQSR Step 3 [VariantRecalibrator] - INDELs
java -jar $GATK -R $GENOMEREF -T VariantRecalibrator -input $temp/${BATCH}_hardFiltered.vcf -mode INDEL -recalFile $temp/${BATCH}_INDEL.recal -tranchesFile $temp/${BATCH}_INDEL.tranches -rscriptFile $temp/${BATCH}_INDEL.plots.R $VQSRargsINDEL $VQSRresourcesINDEL -L $INTERVALS
#VQSR Step 4 [ApplyRecalibration] - INDELs
java -jar $GATK -R $GENOMEREF -T ApplyRecalibration -mode INDEL --ts_filter_level $VQSRFILTERLEVEL -input $temp/${BATCH}_hardFiltered_SNP.recal.snps.vcf -tranchesFile $temp/${BATCH}_INDEL.tranches -recalFile $temp/${BATCH}_INDEL.recal -o $temp/${BATCH}_hardFiltered_SNP.recal.snps.indel.vcf -L $INTERVALS

#STEP 13 - let's flag variants with too many alleles #I typically use 8 as a cutoff?
$VCFtools --vcf $temp/${BATCH}_hardFiltered_SNP.recal.snps.indel.vcf --max-alleles 8 --recode --recode-INFO-all --out $temp/${BATCH}_hardFiltered_sorted.vcf_alleleReduction.vcf

#STEP 14 - let's flag variant that fail our VQSR filtering, and produce our final VCF file
java -jar $GATK -R $GENOMEREF -T SelectVariants -V $temp/${BATCH}_hardFiltered_sorted.vcf_alleleReduction.vcf.recode.vcf -o $temp/${BATCH}_hardFiltered_sorted.vcf_alleleReduction_VQSR.vcf -L $INTERVALS

#STEP 15 - let's remove all temporary files
#rm $temp/${BATCH}_*

exit

rm $OUT/${BATCH}_INDEL.recal
rm $OUT/${BATCH}_SNP.recal
rm $OUT/${BATCH}_SNP.tranches
rm $OUT/${BATCH}_INDEL.plots.R
rm $OUT/${BATCH}_INDEL.recal
rm $OUT/${BATCH}_INDEL.recal.idx
rm $OUT/${BATCH}_INDEL.tranches
rm $OUT/${BATCH}.plots.R
rm $OUT/${BATCH}.recal.snps.indel_alleleReduction.vcf.log
rm $OUT/${BATCH}.recal.snps.indel_alleleReduction.vcf.recode.vcf
rm $OUT/${BATCH}.recal.snps.indel_alleleReduction.vcf.recode.vcf.idx
rm $OUT/${BATCH}.recal.snps.indel.vcf
rm $OUT/${BATCH}.recal.snps.indel.vcf.idx
rm $OUT/${BATCH}.recal.snps.indel.vcf.vcfidx
rm $OUT/${BATCH}_SNP.recal
rm $OUT/${BATCH}_SNP.recal.idx
rm $OUT/${BATCH}_SNP.tranches
rm $OUT/${BATCH}_SNP.tranches.pdf
rm $OUT/${BATCH}_HF4.recal.snps.indel_alleleReduction.vcf.recode.vcf
rm $OUT/${BATCH}_HF4.recal.snps.indel_alleleReduction.vcf 
rm $OUT/${BATCH}_HF4_SNP.recal.snps.indel.vcf
rm $OUT/${BATCH}_HF4_SNP.recal.snps.vcf
rm $OUT/${BATCH}_HF4.vcf
rm $OUT/${BATCH}_HF3.vcf
rm $OUT/${BATCH}_HF2.vcf
rm $OUT/${BATCH}_HF1.vcf 
rm $OUT/${BATCH}_HF1.vcf.idx
rm $OUT/${BATCH}_HF2.vcf.idx
rm $OUT/${BATCH}_HF3.vcf.idx
rm $OUT/${BATCH}_HF4.recal.snps.indel_alleleReduction.vcf.recode.vcf.idx
rm $OUT/${BATCH}_HF4.vcf.idx
rm $OUT/${BATCH}_HF4_SNP.recal.snps.indel.vcf.idx
rm $OUT/${BATCH}_HF4_SNP.recal.snps.vcf.idx




#housekeeping

rm $temp/${BATCH}_gVCF_inputs.txt

rm $temp/$BATCH.g.vcf
rm $temp/$BATCH.g.vcf.idx
rm $temp/${BATCH}_SNP.recal.snps.vcf
rm $temp/${BATCH}_SNP.recal.snps.vcf.idx

#remove the chromosome .vcf files
#place '#' in front of the files you wish to keep

rm $temp/${BATCH}_1.vcf
rm $temp/${BATCH}_2.vcf
rm $temp/${BATCH}_3.vcf
rm $temp/${BATCH}_4.vcf
rm $temp/${BATCH}_5.vcf
rm $temp/${BATCH}_6.vcf
rm $temp/${BATCH}_7.vcf
rm $temp/${BATCH}_8.vcf
rm $temp/${BATCH}_9.vcf
rm $temp/${BATCH}_10.vcf
rm $temp/${BATCH}_11.vcf
rm $temp/${BATCH}_12.vcf
rm $temp/${BATCH}_13.vcf
rm $temp/${BATCH}_14.vcf
rm $temp/${BATCH}_15.vcf
rm $temp/${BATCH}_16.vcf
rm $temp/${BATCH}_17.vcf
rm $temp/${BATCH}_18.vcf
rm $temp/${BATCH}_19.vcf
rm $temp/${BATCH}_20.vcf
rm $temp/${BATCH}_21.vcf
rm $temp/${BATCH}_22.vcf
rm $temp/${BATCH}_X.vcf
rm $temp/${BATCH}_Y.vcf

rm $temp/${BATCH}_1.vcf.idx
rm $temp/${BATCH}_2.vcf.idx
rm $temp/${BATCH}_3.vcf.idx
rm $temp/${BATCH}_4.vcf.idx
rm $temp/${BATCH}_5.vcf.idx
rm $temp/${BATCH}_6.vcf.idx
rm $temp/${BATCH}_7.vcf.idx
rm $temp/${BATCH}_8.vcf.idx
rm $temp/${BATCH}_9.vcf.idx
rm $temp/${BATCH}_10.vcf.idx
rm $temp/${BATCH}_11.vcf.idx
rm $temp/${BATCH}_12.vcf.idx
rm $temp/${BATCH}_13.vcf.idx
rm $temp/${BATCH}_14.vcf.idx
rm $temp/${BATCH}_15.vcf.idx
rm $temp/${BATCH}_16.vcf.idx
rm $temp/${BATCH}_17.vcf.idx
rm $temp/${BATCH}_18.vcf.idx
rm $temp/${BATCH}_19.vcf.idx
rm $temp/${BATCH}_20.vcf.idx
rm $temp/${BATCH}_21.vcf.idx
rm $temp/${BATCH}_22.vcf.idx
rm $temp/${BATCH}_X.vcf.idx
rm $temp/${BATCH}_Y.vcf.idx

rm $temp/${BATCH}_INDEL.recal
rm $temp/${BATCH}_INDEL.tranches

rm $temp/${BATCH}_SNP.recal
rm $temp/${BATCH}_SNP.tranches


rm $temp/${BATCH}_INDEL.plots.R
rm $temp/${BATCH}_INDEL.recal
rm $temp/${BATCH}_INDEL.recal.idx
rm $temp/${BATCH}_INDEL.tranches
rm $temp/${BATCH}.plots.R
rm $temp/${BATCH}.recal.snps.indel_alleleReduction.vcf.log
rm $temp/${BATCH}.recal.snps.indel_alleleReduction.vcf.recode.vcf
rm $temp/${BATCH}.recal.snps.indel_alleleReduction.vcf.recode.vcf.idx
#rm $temp/${BATCH}.recal.snps.indel_alleleReduction_VQSRfiltered.vcf # keep maybe ?
#rm $temp/${BATCH}.recal.snps.indel_alleleReduction_VQSRfiltered.vcf.idx # # keep maybe ?
rm $temp/${BATCH}.recal.snps.indel.vcf
rm $temp/${BATCH}.recal.snps.indel.vcf.idx
rm $temp/${BATCH}.recal.snps.indel.vcf.vcfidx
rm $temp/${BATCH}_SNP.recal
rm $temp/${BATCH}_SNP.recal.idx
rm $temp/${BATCH}_SNP.tranches
rm $temp/${BATCH}_SNP.tranches.pdf

