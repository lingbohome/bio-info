#!/bin/bash

##### Description: Use for post alignment process with input bam file of ctDNA acquied from alignment process.
##### Usage: sh aln_post_ctDNA.sh -b [Sorted bam file] -t [Target region bed file] -p [Sample name] -o [Output directory]
##### Requirements: 
#####     Databases required and the environment variables are as follows:
#####         $REF: Reference genome sequences(hg19, fasta format)
#####     Softwares of gencore-v0.14.0, picard-2.20.0 and samtools-v0.1.19 were installed, and the environment variables of executable program were set as follows:
#####         $GENCORE: gencore-v0.14.0
#####         $PICARD: picard-2.20.0
#####         $SAMTOOLS_1_8: samtools-v1.8

while getopts b:t:m:p:o: opt
do
    case $opt in
        b)
        bam=$OPTARG
        ;;
        t)
        target=$OPTARG
        ;;
        m)
        medicineDir=$OPTARG
        ;;
        p)
        sampleID=$OPTARG
        ;;
        o)
        outDir=$OPTARG
        ;;
    esac
done

source hg19_env.sh

###例子
### sh /home/wangce/workdir/bin/ctDNA_v5/snvIndel_calling_ctDNA_v2.sh  -t /data3/Projects/panel338_2021/240428_A00582_1255_AHWWNGDSX7/02_aln/202411_CL02080.rmdup.target.bam  -l /home/wangce/workdir/database/humandb/panel/338genes_20210624.bed  -c /home/wangce/workdir/database/humandb/ctDNA_report/pon/panel338_Mutect2_2021/vcf/40samples.pon.vcf.gz  -p 202411_CL02080  -o /data3/Projects/panel338_2021/240428_A00582_1255_AHWWNGDSX7/03_var  >>/data3/Projects/panel338_2021/240428_A00582_1255_AHWWNGDSX7/scriptLog/03_snvInDel_calling_202411_CL02080.o  2>>/data3/Projects/panel338_2021/240428_A00582_1255_AHWWNGDSX7/scriptLog/03_snvInDel_calling_202411_CL02080.e
echo "*** Post alignment process by gencore-v0.14.0, picard-2.20.0 and samtools-v0.1.19 ***"
echo "Begin at:" `date`
{
    ##### Generate consensus reads to reduce sequencing noises and remove duplications 
    $GENCORE -i $bam -o $outDir/$sampleID.rmdup.unsorted.bam -r $REF -a 0.7 -j $outDir/$sampleID.report.json -h $outDir/$sampleID.report.html
    ##### Sort bam
    $SAMTOOLS_1_8 sort -m 4G -@ 4 -T $outDir/$sampleID.rmdup -o $outDir/$sampleID.rmdup.bam $outDir/$sampleID.rmdup.unsorted.bam
    $SAMTOOLS_1_8 index $outDir/$sampleID.rmdup.bam
    ##### Remove bam
    rm $outDir/$sampleID.rmdup.unsorted.bam
    ##### Target bam
    $SAMTOOLS_1_8 view -b -@ 4 -L $target -o $outDir/$sampleID.rmdup.target.bam $outDir/$sampleID.rmdup.bam
    $SAMTOOLS_1_8 index $outDir/$sampleID.rmdup.target.bam
} &
{
    ##### Remove duplicate reads
    java -Xmx16g -jar $PICARD MarkDuplicates I=$bam O=$outDir/$sampleID.sv_rmdup.bam M=$outDir/$sampleID.sv_rmdup.metrix VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true
    $SAMTOOLS_1_8 index $outDir/$sampleID.sv_rmdup.bam
    ##### Target bam
    #$SAMTOOLS_1_8 view -b -@ 4 -L $medicineDir/fusion_gene.bed -o $outDir/$sampleID.sv_rmdup.sv_target.bam $outDir/$sampleID.sv_rmdup.bam
    #$SAMTOOLS_1_8 index $outDir/$sampleID.sv_rmdup.sv_target.bam
} &
wait
echo "End at:" `date`
echo "*** Finished post alignment process by gencore-v0.14.0, picard-2.20.0 and samtools-v1.8 ***"
