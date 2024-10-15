#!/bin/bash

##### Description: Use for alignment process with input paired-end clean data files of ctDNA acquired from quality control process.
##### Usage: sh aln_ctDNA.sh -f [R1 clean fastq file] -r [R2 clean fastq file] -p [Sample name] -t [Target region bed file] -o [Output directory]
##### Requirements: 
#####     Fasta file of reference genome sequences(hg19) was downloaded, and the environment variable of the file was set to $REF.
#####     Softwares of BWA-0.7.12-r1039, picard-2.20.0 and samtools-v0.1.19 were installed, and the environment variables of executable program were set as follows:
#####         $PICARD: picard-2.20.0
#####         $SAMTOOLS_0_1_19: samtools-v0.1.19

while getopts f:r:p:t:o: opt
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
        t)
        target=$OPTARG
        ;;
        o)
        outDir=$OPTARG
        ;;
    esac
done

##加载环境变量
source hg19_env.sh
######例子
####   sh /home/wangce/workdir/bin/ctDNA_v5/aln_ctDNA.sh  -f /data3/Projects/panel338_2021/240428_A00582_1255_AHWWNGDSX7/01_QC/202411_CL02080.R1_clean.fq.gz  -r /data3/Projects/panel338_2021/240428_A00582_1255_AHWWNGDSX7/01_QC/202411_CL02080.R2_clean.fq.gz  -t /home/wangce/workdir/database/humandb/panel/338genes_20210624.bed  -p 202411_CL02080  -o /data3/Projects/panel338_2021/240428_A00582_1255_AHWWNGDSX7/02_aln  >>/data3/Projects/panel338_2021/240428_A00582_1255_AHWWNGDSX7/scriptLog/02_aln_202411_CL02080.o  2>>/data3/Projects/panel338_2021/240428_A00582_1255_AHWWNGDSX7/scriptLog/02_aln_202411_CL02080.e
echo "Begin at:" `date`
echo "*** Alignment by BWA-0.7.12-r1039 ***"
bwa mem -R "@RG\tID:$sampleID\tSM:$sampleID\tPL:illumina\tDS:hg19\tCN:clgene" -t 8 -M $REF $r1_fq $r2_fq > $outDir/$sampleID.bam
echo "*** Finished alignment by BWA-0.7.12-r1039 ***"
echo "*** Bam sorting by samtools-v1.8 ***"
$SAMTOOLS_1_8 sort -m 4G -@ 8 -T $outDir/$sampleID.sort -o $outDir/$sampleID.sort.bam $outDir/$sampleID.bam
if [ -f "$outDir/$sampleID.sort.bam" ];then
	$SAMTOOLS_0_1_19 index $outDir/$sampleID.sort.bam
	rm $outDir/$sampleID.bam
fi
echo "*** Finished sorting bam by samtools-v1.8 ***"
echo "*** Target region filter by samtools-v0.1.19 ***"
$SAMTOOLS_0_1_19 view -b -@ 8 -q 30 -L $target -o $outDir/$sampleID.target.bam $outDir/$sampleID.sort.bam
$SAMTOOLS_0_1_19 index $outDir/$sampleID.target.bam
echo "*** Finished target region filter by samtools-v0.1.19 ***"
echo "End at:" `date`
