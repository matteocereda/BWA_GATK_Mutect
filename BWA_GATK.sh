#!/bin/sh

## Set the working Directorty to be the current one, i.e. where you are submitting the script
#$ -cwd

#$ -j y

## Job Name, can be anything##

## Set the SHELL type to be sh##
#$ -S /bin/sh


# -------------------------------------------------------------------------------------------------#

function count_raw_reads {
	printf "## Raw reads in fastq (not filtered)\n" >> $1
	IFS=$'\x0A'
	# R1
	names=( $(ls -dtr $2) )
	rawreads=0
	for file in ${names[@]}
	do
	printf "${file}\t" >> $1
	cnt=$(zcat $file | grep $4 | wc -l)
	printf "${cnt}\n" >> $1
	rawreads=$(($rawreads + $cnt))
	done;
	printf "# TOT RAW READS R1:\t ${rawreads}\n" >> $1
	printf "\n" >> $1
	# R2
	names=( $(ls -dtr $3) )
	rawreads=0
	for file in ${names[@]}
	do
	printf "${file}\t" >> $1
	cnt=$(zcat $file | grep $4 | wc -l)
	printf "${cnt}\n" >> $1
	rawreads=$(($rawreads + $cnt))
	done;
	printf "# TOT RAW READS R2:\t ${rawreads}\n" >> $1
	printf "\n" >> $1
}


function filer_raw_reads {
	local path=$1
	local ffiltered=$2
	local flog=$3
	local HISEQ=$4
	local R1_reads=( $(ls -dtr *"R1"*) )
	local R2_reads=( $(ls -dtr *"R2"*) )
	local old_dir=$(pwd)
	cd $path
	printf "## Filtered reads in fastq\n" >> $flog
	printf "$R1_reads\t" >> $flog
	zgrep -A 3 '^@.* [^:]*:N:[^:]*:' "$path/$R1_reads" | grep -v -- '^--$' | gzip -c > "$path/$ffiltered/filtered_$R1_reads" ;
	cnt=$(zgrep $HISEQ "$path/$ffiltered/filtered_$R1_reads" | wc -l)
	printf "${cnt}\n" >> $flog
	printf "# TOT FILTERED READS R1:\t ${cnt}\n" >> $flog
	printf "\n" >> $flog
	printf "$R2_reads\t" >> $flog
	zgrep -A 3 '^@.* [^:]*:N:[^:]*:' "$path/$R2_reads" | grep -v -- '^--$' | gzip -c > "$path/$ffiltered/filtered_$R2_reads" ;
	cnt=$(zgrep $HISEQ "$path/$ffiltered/filtered_$R2_reads" | wc -l)
	printf "${cnt}\n" >> $flog
	printf "# TOT FILTERED READS R2:\t ${cnt}\n" >> $flog
	printf "\n" >> $flog
	cd $old_dir
	echo $old_dir
}

# -------------------------------------------------------------------------------------------------#
  
# CONFIGURATION
# --------------

sample=$1 # Name of the sample in input
folder=$2

module load bwa/0.7.12 
module load samtools/0.1.18
module load bedtools
module load R
module load picard-tools/1.74
module load jdk/1.7.0_07 
module load perl

# where are the fastq files
path="/home/FC/ClonalExpansion/WholeExome/$folder/$sample"

# Where the results will be stored
root="/home/FC/ClonalExpansion/WholeExome/$folder/$sample"

# On target configurations:
ftarget='/home/FC/Capture_kits/Human_V5/S04380110/S04380110_Covered.bed'

# ref genome
genome_ref='/home/FC/DB/Genomes/Hs_v37/hg19.fa'

# samtools
samtls='samtools'
GATK="java -Xmx16g -jar /home/ceredam/Software/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar"
DBSNP138="/home/ceredam/CRC/Mutect/dbsnp_138.hg19.vcf"
PICARD="java -Xmx16g -jar /apps/picardtools/1.74"

# END configuration


# MAIN
# ------
# aln folder
align='mutect'

# filtered fastq folder
ffiltered='filtered-fastq'

OLD_IFS=$IFS
mkdir -p $root
cd $root

mkdir -p "$path/$ffiltered"
mkdir -p "$root/$align"


HISEQ='HISEQ'

current=$(date +"%Y%m%d")
flog="${root}/$align/${current}_$sample.log"

# We are assunming ONLY two fastq files!!!
R1_reads=( $(ls -dtr "$path/"*"R1"*) )
R2_reads=( $(ls -dtr "$path/"*"R2"*) )
filer_raw_reads $path $ffiltered $flog $HISEQ


# # we move in the aln directory
cd "$root/$align"

fout="$sample.sam"

# get the fastq filterd file paths
R1_reads=( $(ls -dtr "$path/$ffiltered/"*"R1"*) )
R2_reads=( $(ls -dtr "$path/$ffiltered/"*"R2"*) )

# echo ".. Allinment and stats (Step 02)"
# bwa index -a bwtsw $genome_ref 
bwa mem -p -M -t 6 -R "@RG\tID:$sample\tSM:$sample" -v 1 $genome_ref $R1_reads $R2_reads > $fout

$samtls view -S -F 4 -b -u -L $ftarget -o $fout.bam $fout
# java -jar $PICARD/SamFormatConverter I=$fout O=$fout.bam CREATE_INDEX=TRUE

$PICARD/AddOrReplaceReadGroups.jar I=$fout.bam O=$fout.RG.bam RGID=$sample RGLB=$sample RGPL="ILLUMINA" RGPU="unkn-0.0" RGSM=$sample VALIDATION_STRINGENCY=LENIENT 
$PICARD/SortSam.jar I=$fout.RG.bam O=$fout.sort.bam SO=coordinate
$PICARD/ReorderSam.jar I=$fout.sort.bam O=$fout.sort.reorder.bam R=${genome_ref} CREATE_INDEX=true

# Remove dublicates
$PICARD/MarkDuplicates.jar I=$fout.sort.reorder.bam O=$fout.dedupped.bam M=dedup_metrics.txt CREATE_INDEX=true	

#Run Picard metrics collection tools, interpret output

$PICARD/CollectGcBiasMetrics.jar R=$genome_ref I=$fout.dedupped.bam O=gcbias_metrics.txt CHART=gcbias_metrics.pdf S=gcbias_summ_metrics.txt
$PICARD/CollectAlignmentSummaryMetrics.jar R=$genome_ref I=$fout.dedupped.bam O=alignment_summ_metrics.txt 

# up to here ok


# Base Recalibration

#A. Build the recalibration model

$GATK -T BaseRecalibrator -R $genome_ref -I $fout.dedupped.bam  -knownSites $DBSNP138 -o recalibration.table 

#B. Apply the recalibration model

$GATK -T PrintReads -R $genome_ref -I $fout.dedupped.bam -BQSR recalibration.table -o $fout.dedupped.recal.bam

