#!/bin/sh

## Set the working Directorty to be the current one, i.e. where you are submitting the script
#$ -cwd

#$ -j y

## Job Name, can be anything##
#$ -N Mut.uh6_T

## Set the SHELL type to be sh##
#$ -S /bin/sh

#$ -M matteo.cereda@kcl.ac.uk
#$ -m e

module load jdk/1.7.0_07 

sample=$1
normal=$2
folder=$3

bam="/home/FC/ClonalExpansion/WholeExome/$folder/$sample/mutect/$sample.sam.dedupped.recal.bam"
Nbam="/home/FC/ClonalExpansion/WholeExome/$folder/$normal/mutect/$normal.sam.dedupped.recal.bam"

outfolder="/home/FC/ClonalExpansion/WholeExome/$folder/$sample/mutect"


MUTECT="java -Xmx16g -jar /home/ceredam/Software/mutect-1.1.7.jar --analysis_type MuTect  --enable_extended_output " #--useOriginalQualities  --fraction_contamination --fix_misencoded_quality_scores --max_alt_allele_in_normal_count 2
COSMIC="/home/FC/DB/COSMIC/74/Cosmic.hg19.vcf"
DBSNP138="/home/ceredam/CRC/Mutect/dbsnp_138.hg19.vcf"
genome_ref='/home/FC/DB/Genomes/Hs_v37/hg19.fa'

$MUTECT --reference_sequence ${genome_ref} --cosmic ${COSMIC} --dbsnp ${DBSNP138} --input_file:normal ${Nbam} --input_file:tumor ${bam} --out ${outfolder}/${sample}.${normal}.mutect.1.17 

