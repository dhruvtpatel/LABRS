#!/bin/bash
#MSUB -A b1042
#MSUB -l naccesspolicy=singlenode
#MSUB -l walltime=48:00:00 
#MSUB -q genomics
#MSUB -o output_files/
#MSUB -e output_files/
#MSUB -l nodes=1:ppn=1

#pipeline to align WGS data from HiSeq FastQ files
#November 2015 - modified to run on b1049 (September 2016)
#September 2017 - modified to run on the compute nodes, launching jobs via qsub  - jon
#splits alignment into 8 jobs (by HiSeq lane)
#includes all .bam file processing
#end result includes multiple (small) .bam files per sample ready for variant calling

#Example usage:
#note options! input dir, output dir and sample ID! press <ENTER> 
#sh WGS_alignment_pipeline.sh -i /data/kronos/apittman/Exome_Data_Analysis/Whole_Genome_Sequencing_Pipeline/Unaligned -o /data/kronos/apittman/Exome_Data_Analysis/Whole_Genome_Sequencing_Pipeline/Unaligned -l sample_xxx

core_count=24
#######################################################################################

##NOTE: easiest option is to change the names of your FastQs to reflect that shown below

#######################################################################################

#Enter the path to your input and output folders
#Enter Sample ID to process --> should have folder with this ID containing relevant FastQs

# Input folder
iFolder=/projects/b1049/genomes/JS/RawData
# Output folder
oFolder=/projects/b1042/LubbeLab/testtemp/genomes/JS20
myIDs=PD26700

########################################################################################

bwa="/projects/b1049/genetics_programs/bwa-0.7.12/bwa"
indexedgenome="/projects/b1049/genetics_refs/novoalign/human_g1k_v37.fasta.k15.s2.novoindex"
samtools="/projects/b1049/genetics_programs/samtools/samtools"
picard="/projects/b1049/genetics_programs/picard/picard.jar"
gatk="/projects/b1049/genetics_programs/gatk/GenomeAnalysisTK.jar"
genomeFASTA="/projects/b1049/genetics_refs/fasta/human_g1k_v37.fasta"
genomeFASTAI="/projects/b1049/genetics_refs/fasta/human_g1k_v37.fastai"
knownINDELS="/projects/b1049/genetics_refs/Mills_and_1000G_gold_standard.indels.hg19_modified.vcf"
samtoolsMEM="5000000000" 
picardArguments="TMP_DIR=/projects/b1042/LubbeLab/testtemp ASSUME_SORTED=TRUE REMOVE_DUPLICATES=FALSE VALIDATION_STRINGENCY=LENIENT" 

########################################################################################

echo `date` >> $oFolder/output.log

echo "
	====================================
	Whole Genome Sequencing Alignment 
	====================================
	"	>> $oFolder/output.log



echo "
	===================================
	STEP ONE - Checking for FastQ files
	===================================
	" >> $oFolder/output.log
	
########################################################################################

#while getopts i:o:l: option
#do
#        case "${option}"
#        in
#               i) iFolder=${OPTARG};; #as defined by your command line
#               o) oFolder=${OPTARG};;
#				l) myIDs=${OPTARG};;

#  esac

#done

if [ ! -e $oFolder ]; then mkdir $oFolder; echo "Making output directory" >> $oFolder/output.log; fi #checking if output folder exists, and creating it if not found

########################################################################################

n=0
while (( $n <= 7 ))
	do
	n=$(( n+1 ))

##NOTE: easiest option is to change the names of your FastQs to reflect that shown below

if [ -e ${iFolder}/${myIDs}/*L${n}_1*gz ]; then	#-e ${iFolder}/${myIDs}/*L00${n}_R1*gz
	
	echo "Files exist" >> $oFolder/output.log
	echo "There are FASTQ files from lane ${n}" >> $oFolder/output.log

	for nID in $myIDs; do

		folder=${iFolder}/${nID}
		
		for file1 in `find $folder -name *L${n}_1*gz`; do	#for file1 in `find $folder -name *L00${n}_R1*gz`; do
		((nfiles=nfiles+2))
		file2=`echo $file1 | sed -e 's/_1/_2/g'`
		inputFiles=`echo "<(zcat "$file1") <(zcat "$file2")"` >> $oFolder/output.log
		
		echo $inputFiles >> $oFolder/output.log
		
		output=${oFolder}/${nID}/${nID}
		
		if [ ! -e ${oFolder}/${nID} ]; then mkdir ${oFolder}/${nID}; fi 



		
echo "
	===================================
	STEP TWO - Writing the scripts
	===================================
	" >> $oFolder/output.log

echo "
	done !
	" >> $oFolder/output.log

		echo "#!/bin/bash
			#PBS -A b1042
			#PBS -l naccesspolicy=singlenode
			#PBS -l walltime=48:00:00 
			#PBS -q genomics
                        #PBS -l nodes=1:ppn=$core_count
                        
			cd $PBS_O_WORKDIR
			##first need to load the correct version of java (v1.8)
			module load java
			"  > $oFolder/J_script_${nID}_${n}.sh	#do the number of nodes specified here (n=4) match what is defined in script output below as set out above? #change walltime to match length of time for sh to run
		
		echo "#These are the pipeline steps for sample ${nID} of lane ${n} of the HiSeq"  >> $oFolder/J_script_${nID}_${n}.sh
			
##now lets write to script each step in the pipeline:		
##"$bwa mem $inputFiles $genomeFASTA -R $'@RG\tID:${nID}\tSM:${nID}\tLB:${nID}\tPL:ILLUMINA' -t 4 -M > ${output}_${n}.sam" >> $oFolder/J_script_${nID}_${n}.sh
##"$novoalign -o SAM $'@RG\tID:${nID}\tSM:${nID}\tLB:${nID}\tPL:ILLUMINA' -f $inputFiles -d $indexedgenome $novoalignArguments > ${output}_${n}.sam" >> $oFolder/J_script_${nID}_${n}.sh
		
		echo " " >> $oFolder/J_script_${nID}_${n}.sh	#this step aligns reads in the FastQ to the reference genome to create a SAM file
		#Commented out by jon - getting different results with different core counts - 
		echo "$bwa mem -M -t $core_count -R $'@RG\tID:${nID}\tSM:${nID}\tLB:${nID}\tPL:ILLUMINA' $genomeFASTA $inputFiles > ${output}_${n}.sam" >> $oFolder/J_script_${nID}_${n}.sh
                #Trying this with novoalign instead -jon
                #echo "$novoalign -o SAM $'@RG\tID:${nID}\tSM:${nID}\tLB:${nID}\tPL:ILLUMINA' -f $inputFiles -d $indexedgenome $novoalignArguments > ${output}_${n}.sam" >> $oFolder/J_script_${nID}_${n}.sh #end novoalign insert -jon
		echo " " >> $oFolder/J_script_${nID}_${n}.sh	#this step creates a BAM file
		echo "$samtools view -bS -t $genomeFASTAI -@$core_count -o ${output}_${n}.bam ${output}_${n}.sam" >> $oFolder/J_script_${nID}_${n}.sh
		echo " " >> $oFolder/J_script_${nID}_${n}.sh	#this step sorts the BAM file
		echo "$samtools sort -m $samtoolsMEM -@$core_count ${output}_${n}.bam -o ${output}_${n}_sorted.bam" >> $oFolder/J_script_${nID}_${n}.sh
		echo " " >> $oFolder/J_script_${nID}_${n}.sh	#this creates an index file for the BAM file
		echo "$samtools index ${output}_${n}_sorted.bam" >> $oFolder/J_script_${nID}_${n}.sh
		echo " " >> $oFolder/J_script_${nID}_${n}.sh	#
		echo "java -Xmx10g -jar $picard MarkDuplicates $picardArguments INPUT=${output}_${n}_sorted.bam OUTPUT=${output}_${n}_sorted_unique.bam METRICS_FILE=${output}_${n}_picard_metrics.out" >> $oFolder/J_script_${nID}_${n}.sh
		echo " " >> $oFolder/J_script_${nID}_${n}.sh
		echo "$samtools index ${output}_${n}_sorted_unique.bam" >> $oFolder/J_script_${nID}_${n}.sh
		echo " " >> $oFolder/J_script_${nID}_${n}.sh
		echo java -jar $gatk -T RealignerTargetCreator -nt $core_count -R $genomeFASTA -o ${output}_${n}_sorted_unique.bam.list -I ${output}_${n}_sorted_unique.bam --known $knownINDELS >> $oFolder/J_script_${nID}_${n}.sh #ensure -nt equals the number of cores specified using MSUB
		echo " " >> $oFolder/J_script_${nID}_${n}.sh
		echo java -jar $gatk -T IndelRealigner -nt $core_count -R $genomeFASTA -targetIntervals ${output}_${n}_sorted_unique.bam.list -I ${output}_${n}_sorted_unique.bam -o ${output}_${n}_sorted_unique_realigned.bam --knownAlleles $knownINDELS >> $oFolder/J_script_${nID}_${n}.sh
		echo " " >> $oFolder/J_script_${nID}_${n}.sh
		echo "$samtools index ${output}_${n}_sorted_unique_realigned.bam" >> $oFolder/J_script_${nID}_${n}.sh
##removing intermediate files:	
	
		echo "
	
			###Remove Intermediate BAM files leaving only final processed_GATK_bam file:

			if [ -e ${output}_${n}_sorted_unique_realigned.bam ] 
				then
					rm ${output}_${n}.bam
					rm ${output}_${n}.bam.bai
					rm ${output}_${n}_sorted.bam
					rm ${output}_${n}_sorted.bam.bai
					rm ${output}_${n}.sam
					rm ${output}_${n}_sorted_unique.bam
					rm ${output}_${n}_sorted_unique.bam.bai
					rm ${output}_${n}_sorted_unique.bam.list
					rm ${output}_${n}_sorted_unique_realigned.bai
					echo 'old BAM and intermediate BAM.bai all removed'
				else
					echo '! something has gone wrong '
				fi

		

			" >> $oFolder/J_script_${nID}_${n}.sh
			
			
			echo "exit" >> $oFolder/J_script_${nID}_${n}.sh
	
		done
		
	done
		
	else echo "There are NO FASTQ files from lane ${n}" >> $oFolder/output.log
fi	

done


##now last (but not least) we will automatically submit all jobs to the SGE of b1049

	echo "
	=======================================
	STEP THREE - Submitting jobs to QUEST
	=======================================
	" >> $oFolder/output.log
	
n=0
while (( $n <= 7 ))
	do
	n=$(( n+1 ))

		for nID in $myIDs; do

			if [ -e ${iFolder}/${myIDs}/*L${n}_1*gz ]; then
				sleep 2
				qsub $oFolder/J_script_${nID}_${n}.sh
				
			else 
			echo "There is NO job from lane ${n}" >> $oFolder/output.log

			fi
		done	

done

exit

