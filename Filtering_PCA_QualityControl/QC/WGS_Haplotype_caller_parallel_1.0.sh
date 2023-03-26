#!/bin/bash

#pipeline to align WGS data from HiSeq FastQ files
#November 2015 - modified to run on Quest
#HaplotypeCaller genomewide
#Input multiple .bam files

#if using -L ${chr} change to -o $oFolder/${myIDs}/${myIDs}.raw.snps.indels.{chr}.g.vcf
#note options! input dir, output dir and sample ID! pres 
#sh WGS_Haplotype_caller_pipeline.sh -i /data/kronos/apittman/Exome_Data_Analysis/Whole_Genome_Sequencing_Pipeline/Unaligned -o /data/kronos/apittman/Exome_Data_Analysis/Whole_Genome_Sequencing_Pipeline/WGS_aligned -l P57-04

iFolder=/projects/b1049/Niccolo_NGS/genomes
oFolder=/projects/b1042/LubbeLab/Haplo_Parallel_test

myIDs=SS4009022

########################################################################################

novoalign="/projects/b1049/genetics_programs/novocraft/novoalign"
novoalignArguments="--rOQ --hdrhd 3 -H -k -o SoftClip -t 320 -F ILM1.8" # -c 4 specifies multi-threading but it's disabled in free version of novoalign
bwa="/projects/b1049/genetics_programs/bwa-0.7.12/bwa"
indexedgenome="/projects/b1049/genetics_refs/novoalign/human_g1k_v37.fasta.k15.s2.novoindex"
samtools="/projects/b1049/genetics_programs/samtools/samtools"
picard="/projects/b1049/genetics_programs/picard/picard.jar"
gatk="/projects/b1049/genetics_programs/gatk/GenomeAnalysisTK.jar"
genomeFASTA="/projects/b1049/genetics_refs/fasta/human_g1k_v37.fasta"
genomeFASTAI="/projects/b1049/genetics_refs/fasta/human_g1k_v37.fastai"
knownINDELS="/projects/b1049/genetics_refs/Mills_and_1000G_gold_standard.indels.hg19_modified.vcf"
samtoolsMEM="5000000000" 
picardArguments="TMP_DIR=/projects/b1042/LubbeLab/tempHap ASSUME_SORTED=TRUE REMOVE_DUPLICATES=FALSE VALIDATION_STRINGENCY=LENIENT"

INTERVALS="/projects/b1049/genetics_refs/INTERVALS.bed" #add an interval file if you want (-L) for GATK !!!

########################################################################################

chIDs="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT"

########################################################################################

#while getopts i:o:l: option
#do
#        case "${option}"
#        in
#                i) iFolder=${OPTARG};;
#                o) oFolder=${OPTARG};;
#				l) myIDs=${OPTARG};;
#
# esac

#done

if [ ! -e $oFolder ]; then mkdir $oFolder; echo "making output directory"; fi

########################################################################################

#First we need to make a list of inputs: 

if [  -e ${oFolder}/${myIDs}/${myIDs}_bam.inputs.txt ]; then rm ${oFolder}/${myIDs}/${myIDs}_bam.inputs.txt; fi

output=${oFolder}/${myIDs}/${myIDs}
# changed oFolder to iFolder in below line
BAMFILES=${iFolder}/${myIDs}/${myIDs}_*_realigned.bam

for file in $BAMFILES;do
	echo $BAMFILES
	echo "-I $file \\" >> ${oFolder}/${myIDs}/${myIDs}_bam.inputs.txt
done

########################################################################################

#Haplotype Caller [HaplotypeCaller] # lets generate a gVCF file for each chromosome of each sample individually:


for chr in $chIDs; do
        echo "#!/bin/bash                       
        #MSUB -A b1042
        #MSUB -l naccesspolicy=singlenode
        #MSUB -l walltime=78:00:00 
        #MSUB -l nodes=1:ppn=2
        #MSUB -q genomicslong
        cd $PBS_O_WORKDIR
        ##first need to load the correct version of java (v1.8)
        module load java" > ${oFolder}/${myIDs}/${myIDs}_${chr}_HaplotypeCaller_Job.sh

	echo "java -jar $gatk -R $genomeFASTA -L ${chr} -T HaplotypeCaller --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 `cat ${oFolder}/${myIDs}/${myIDs}_bam.inputs.txt`-o $oFolder/${myIDs}/${myIDs}.raw.snps.indels.${chr}.g.vcf" >> ${oFolder}/${myIDs}/${myIDs}_${chr}_HaplotypeCaller_Job.sh
	echo " " >> ${oFolder}/${myIDs}/${myIDs}_${chr}_HaplotypeCaller_Job.sh
	
#submission to b1042:

msub ${oFolder}/${myIDs}/${myIDs}_${chr}_HaplotypeCaller_Job.sh
sleep 10

done
########################################################################################
########################################################################################
exit

