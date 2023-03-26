#!/bin/bash
#SBATCH -A b1042
#SBATCH -N 1
#SBATCH -t 240:00:00
#SBATCH -p genomicslong
#SBATCH --ntasks-per-node=10
#SBATCH --mem=50000
cd $PBS_O_WORKDIR
#first enter path to and name of VCF file
FILELOCATION=/projects/b1049/brent/QC
FILE="brent_ppmi_chr6"
ANNOVAR=/projects/b1049/genetics_programs/annovar_2017/annovar
CONVERT=/projects/b1049/genetics_programs/annovar_2017/annovar/convert2annovar.pl
TABLE=/projects/b1049/genetics_programs/annovar_2017/annovar/table_annovar.pl
#$CONVERT -format vcf4 $FILELOCATION/$FILE -allsample -withfreq -includeinfo -outfile $FILELOCATION/$FILE.avinput
#for old vcfs i.e. not run on HaplotypeCaller, use the command line below instead
#$CONVERT -format vcf4old $FILELOCATION/$FILE -includeinfo > $FILELOCATION/$FILE.avinput
cd $ANNOVAR
$TABLE $FILELOCATION/$FILE.avinput $ANNOVAR/humandb/ -buildver hg19 -protocol refGene,genomicSuperDups,gnomad_genome,esp6500siv2_all,cadd13gt20,gerp++gt2,dbnsfp33a,cg69,exac03,1000g2015aug_all,avsnp150,clinvar_20170905 -remove -otherinfo -operation g,r,f,f,f,f,f,f,f,f,f,f -nastring . -out $FILELOCATION/brent_chr6_park2__annotated
#$TABLE $FILELOCATION/dk_annovar.txt $ANNOVAR/humandb/ -buildver hg19 -protocol refGene,genomicSuperDups,gnomad_genome,esp6500siv2_all,cadd13gt20,hg19_gerp++gt2,dbnsfp33a,cg69,exac03,1000g2015aug_all,avsnp150,clinvar_20170905 -remove -otherinfo -operation g,r,f,f,f,f,f,f,f,f,f,f -nastring .
cd $FILELOCATION