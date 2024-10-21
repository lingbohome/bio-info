#!/bin/bash
# Name: microsatellite_v2.sh

##### Description: Used for microsatellite status detection with single tumor tissue.
##### Usage: sh microsatellite_v2.sh -n [Constructed normal bam file for MSIsensor] -t [Tumor pickseq bam file for MSIsensor] -T [Tumor bam file] -b [Target region bed file] -r [Baseline file for MSIsensor-pro] -B [Bed file for MSIFinder] -c [Ctr file for MSIFinder] -s [Sample name] -o [Output directory]
##### Requirements: 
#####     MSIsensor
#####     MSIsensor-pro
#####     MSIFinder
#####     microsatellites.list (hg19) file for MSIsensor

while getopts n:t:T:b:r:B:c:s:o: opt
do
    case $opt in
        n)
        normal_bam=$OPTARG
        ;;
        t)
        tumor_pickseq_bam=$OPTARG
        ;;
        T)
        tumor_bam=$OPTARG
        ;;
        b)
        target=$OPTARG
        ;;
        r)
        baseline_for_msisensor_pro=$OPTARG
        ;;
        B)
        bed_for_msifinder=$OPTARG
        ;;
        c)
        ctr_for_msifinder=$OPTARG
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

echo "Begin at:" `date`

echo "*** Analysis by MSIsensor (selected sites) ***"
$MSISENSOR msi -d $MS_V2 -n $normal_bam -t $tumor_pickseq_bam -e $target -o $outDir/$sampleID
msisensor_score=`cat $outDir/$sampleID | tail -1 | cut -f 3`
msisensor_judge=`echo "$msisensor_score > 70" | bc`
echo -e "MSIsensor_Status" > $outDir/$sampleID.msisensor
if [ $msisensor_judge -eq 1 ];then
    echo -e "MSI" >> $outDir/$sampleID.msisensor
else
    echo -e "MSS" >> $outDir/$sampleID.msisensor
fi
paste -d "\t" $outDir/$sampleID $outDir/$sampleID.msisensor > $outDir/$sampleID.msisensor.tmp
rm -rf $outDir/$sampleID.msisensor $outDir/$sampleID $outDir/${sampleID}_germline $outDir/${sampleID}_somatic $outDir/${sampleID}_dis

echo "*** Analysis by MSIsensor-pro ***"
###加载动态库
export LD_LIBRARY_PATH=/data1/workdir/software/anaconda3/lib:$LD_LIBRARY_PATH
$MSISENSOR_PRO pro -d $baseline_for_msisensor_pro -t $tumor_bam -o $outDir/$sampleID
msisensor_pro_score=`cat $outDir/$sampleID | tail -1 | cut -f 3`
msisensor_pro_judge=`echo "$msisensor_pro_score > 20" | bc`
echo -e "MSIsensor-pro_Status" > $outDir/$sampleID.msisensor-pro
if [ $msisensor_pro_judge -eq 1 ];then
    echo -e "MSI" >> $outDir/$sampleID.msisensor-pro
else
    echo -e "MSS" >> $outDir/$sampleID.msisensor-pro
fi
paste -d "\t" $outDir/$sampleID $outDir/$sampleID.msisensor-pro > $outDir/$sampleID.msisensor-pro.tmp
rm -rf $outDir/$sampleID.msisensor-pro $outDir/$sampleID $outDir/${sampleID}_unstable $outDir/${sampleID}_all $outDir/${sampleID}_dis

echo "*** Analysis by MSIFinder ***"
python3 $MSIFINDER -fastahack $FASTAHACK -rf $REF -bam $tumor_bam -bed $bed_for_msifinder -ctr $ctr_for_msifinder -o $outDir/ -p $sampleID
mv $outDir/$sampleID.MSI.xls $outDir/$sampleID.msifinder.tmp
rm -rf $outDir/${sampleID}_MSIscore.xls

echo "*** Summary of results ***"
msisensor_result=`cat $outDir/$sampleID.msisensor.tmp | tail -1`
msisensor_pro_result=`cat $outDir/$sampleID.msisensor-pro.tmp | tail -1`
msifinder_score=`cat $outDir/$sampleID.msifinder.tmp |tail -1 | cut -f 2`
msifinder_judge=`cat $outDir/$sampleID.msifinder.tmp |tail -1 | cut -f 1`
rm -rf $outDir/$sampleID.msisensor.tmp $outDir/$sampleID.msisensor-pro.tmp $outDir/$sampleID.msifinder.tmp
counter=0
if [ $msisensor_judge -eq 1 ];then
    counter=$((counter+1))
fi
if [ $msisensor_pro_judge -eq 1 ];then
    counter=$((counter+1))
fi
if [ $msifinder_judge == "MSI-H" ];then
    counter=$((counter+1))
fi
echo "Final MSI results: $counter of 3"
if [ $counter -ge 2 ];then
    final_judge="MSI"
else
    final_judge="MSS"
fi
echo -e "#Sample\tJudgement\tTotal_Sites(MSIsensor)\tUnstable_Sites(MSIsensor)\tUnstable_Ratio(MSIsensor)\tStatus(MSIsensor)\tTotal_Sites(MSIsensor-pro)\tUnstable_Sites(MSIsensor-pro)\tUnstable_Ratio(MSIsensor-pro)\tStatus(MSIsensor-pro)\tProba(MSIFinder)\tStatus(MSIFinder)" > $outDir/$sampleID.microsatellite.status.xls
echo -e "$sampleID\t$final_judge ($counter of 3)\t$msisensor_result\t$msisensor_pro_result\t$msifinder_score\t$msifinder_judge" >> $outDir/$sampleID.microsatellite.status.xls

echo "End at:" `date`
