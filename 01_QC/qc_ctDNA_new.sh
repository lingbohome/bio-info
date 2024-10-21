#!/bin/bash

##### Description: Use for quality control process with input paired-end sequencing data of ctDNA.
##### Usage: sh qc_ctDNA.sh -f [R1 fastq file] -r [R2 fastq file] -p [Sample name] -o [Output directory]
##### Requirements:
#####     Fasta format file of adapter sequences was prepaired, and the environment variable of the file was set to $ADAPTER.
#####     Softwares  Trimmomatic-0.36 were installed, and the environment variables of executable program were set as follows:
#####     $TRIMMOMATIC: Trimmomatic-0.36
#####     openjdk version "1.8.0_121"
#####     OpenJDK Runtime Environment (Zulu 8.20.0.5-linux64) (build 1.8.0_121-b15)
#####     OpenJDK 64-Bit Server VM (Zulu 8.20.0.5-linux64) (build 25.121-b15, mixed mode)

while getopts f:r:p:o: opt
do
    case $opt in
        f)
        r1_fq=$OPTARG
        ;;
        r)
        r2_fq=$OPTARG
        ;;
        p)
        sampleID=$OPTARG
        ;;
        o)
        outDir=$OPTARG
        ;;
    esac
done
##软件路径
TRIMMOMATIC=trimmomatic-0.36.jar
echo "*** QC by FastQC-v0.11.5 and Trimmomatic-0.36 ***"
echo "Begin at:" `date`
$FASTQC -t 8 -o $outDir $r1_fq $r2_fq
unzip -d $outDir $outDir/${sampleID}_R1*_fastqc.zip
unzip -d $outDir $outDir/${sampleID}_R2*_fastqc.zip
cp $outDir/${sampleID}_R1*_fastqc/fastqc_data.txt $outDir/$sampleID.R1.fastqc_data.txt
cp $outDir/${sampleID}_R2*_fastqc/fastqc_data.txt $outDir/$sampleID.R2.fastqc_data.txt
rm -rf $outDir/${sampleID}_R1*_fastqc $outDir/${sampleID}_R2*_fastqc
java   -jar $TRIMMOMATIC PE -threads 6 -phred33 $r1_fq $r2_fq $outDir/$sampleID.R1_clean.fq.gz $outDir/$sampleID.R1_unpaired.fq.gz $outDir/$sampleID.R2_clean.fq.gz $outDir/$sampleID.R2_unpaired.fq.gz ILLUMINACLIP:$ADAPTER:2:30:10:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:72
echo "End at:" `date`
echo "*** Finished QC by FastQC-v0.11.5 and Trimmomatic-0.36 ***"