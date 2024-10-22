#import matplotlib
#matplotlib.use('Agg')
import pandas as pd
import pysam
import math
from collections import defaultdict
import argparse as apa
import traceback
import re
import sys
import os
import warnings
import configparser
from sklearn.ensemble import RandomForestClassifier
warnings.filterwarnings("ignore")
#sns.set(color_codes=True)


def get_left_right_base(bed, ref, specbp,fastahack):
    bed = pd.read_table(bed)
    left_base, right_base = [], []
    for r, s, e in zip(bed["chr"], bed["start"], bed["end"]):
        cmdLeft = "%s -r %s:%s-%s %s" % (fastahack, r, s - specbp, s - 1, ref)
        cmdRight = "%s -r %s:%s-%s %s" % (fastahack, r, e + 1, e + specbp, ref)
        left_base.append(str.upper(os.popen(cmdLeft).read().strip()))
        right_base.append(str.upper(os.popen(cmdRight).read().strip()))
    bed["left_base"] = left_base
    bed["right_base"] = right_base
    return bed


def msi_analysis(bed2, bam, dp, flank, mapq, sup_reads):
    Total_Reads = []
    Number_of_Peaks = []
    infos = []
    qcs = []
    sam = pysam.AlignmentFile(bam, "rb")
    for r, s, e, left, right, MS in zip(bed2["chr"], bed2["start"], bed2["end"], bed2["left_base"], bed2["right_base"]
            , bed2["MS"]):
        startpos, endpos = s - flank, e + flank
        TStat = defaultdict(int)
        moudle = r'%s((%s)*)%s' % (left, MS, right)
        for read in sam.fetch(r, startpos, endpos):
            if read.is_duplicate is not True and read.is_secondary is not True and read.mapping_quality >= mapq \
                    and read.is_qcfail is not True:
                seq = read.seq
                if re.search(moudle, seq):
                    msLength = len(re.search(moudle, seq).group(1))
                    TStat[msLength - (e - s + 1)] += 1
        Reads = sum(TStat.values())
        Peaks = 0
        info = ""
        for msLength, read in TStat.items():
            # if read / Reads >= cutoff and read >= sup_reads:
            if read >= sup_reads:
                Peaks += 1
                info += str(msLength) + ":" + str(read / Reads) + ":" + str(read) + " "
        Total_Reads.append(Reads)
        Number_of_Peaks.append(Peaks)
        infos.append(info.strip())
        if Reads < dp:
            qcs.append("fail")
        else:
            qcs.append("pass")
    bed2 = bed2[["MSID", "chr", "start", "end"]]
    bed2["qcs"] = qcs
    bed2["Total_Reads"] = Total_Reads
    bed2["Number_of_Peaks"] = Number_of_Peaks
    bed2["IndelLength:AlleleFraction:SupportingCalls"] = infos
    return bed2


def normalization(result, control, sup_reads):
    control = control[["MSID", "Average_Total_Reads"]]
    control.columns = ["MSID", "Baseline_Average_Total_Reads"]
    result = pd.merge(result, control, on="MSID")
    IndelLength_SupportingCalls = []
    peaks = []
    for reads1, reads2, info in zip(result["Total_Reads"],
                                    result["Baseline_Average_Total_Reads"],
                                    result["IndelLength:AlleleFraction:SupportingCalls"]):
        if reads2 == 0:
            ratio = 1
        else:
            ratio = float(reads1/reads2)
        _IndelLength_SupportingCalls = ""
        Peaks = 0
        for i in info.split(" "):
            if str(i) == '':
                pass
            else:
                _IndelLength_SupportingCalls += str(i.split(":")[0]) + ":" + str(int(i.split(":")[-1])/ratio) + " "
                if int(i.split(":")[-1])/ratio >= sup_reads:
                    Peaks += 1
        IndelLength_SupportingCalls.append(_IndelLength_SupportingCalls.strip())
        peaks.append(Peaks)
    result["Normalized_Number_of_Peaks"] = peaks
    result["Normalized_IndelLength:SupportingCalls"] = IndelLength_SupportingCalls
    return result


def parse_msi(result, control_file, pex, sup_reads):
    baseline = pd.read_table(control_file)
    result = normalization(result, baseline, sup_reads)
    result.to_csv(pex + "_MSIscore.xls", index=False, sep="\t")

    ##预测样本是否阳性或者阴性
    result = pd.merge(baseline, result, on="MSID", how="inner")
    y_test = result["Normalized_Number_of_Peaks"].as_matrix().reshape(1, -1)
    result.index = result["MSID"]
    del result["MSID"]

    positive_s = [sample for sample in result.columns if "positive" in sample]
    negative_s = [sample for sample in result.columns if "negative" in sample]
    result = result[positive_s + negative_s]
    x_train = result.T
    y_train = ["MSI-H"] * len(positive_s) + ["MSS"] * len(negative_s)
    finalPredict(x_train, y_train, y_test, pex)
    return result


def finalPredict(x_train, y_train, y_test, pex):
    clf = RandomForestClassifier(n_estimators=30, max_depth=50, random_state=0)
    clf = clf.fit(x_train, y_train)
    proba = clf.predict_proba(y_test)[0][0]
    print(proba)
    f = open(pex + ".MSI.xls", "w")
    f.write("MSI\tproba\n")
    if proba >= 0.6:
        f.write("MSI-H\t" + str(proba) + "\n")
    elif 0.4 <= proba < 0.6:
        f.write("MSI-L\t" + str(proba) + "\n")
    else:
        f.write("MSS\t" + str(proba) + "\n")
    f.close()


def main():
    parser = apa.ArgumentParser(prog="convert")
    parser.add_argument("-rf", "--ref_genome", required=True, type=str,
                        help="referance genome, such as ucsc.hg19.fasta")
    parser.add_argument('-bam', '--bam', required=True, type=str,
                        help='bam file')
    parser.add_argument('-bed', '--bed', required=True, type=str,
                        help='panel8 str bed file')
    parser.add_argument('-ctr', "--control", required=False, type=str,
                        help='control str result file')
    parser.add_argument('-mp', "--multiplier", default=2, required=False, type=float,
                        help='if peaks greater than Average_Number_Peaks + multiplier * Standard_Deviation, '
                             'this loci will be supposed to a unstable loci, default is 2')
    parser.add_argument("-ct", "--cutoff", default=0.01, required=False, type=float,
                        help="if a AlleleFraction greater than this parameter, "
                             "this Allele will be supposed to a peak, default is 0.01")
    parser.add_argument("-sr", "--sup_reads", default=3, required=False, type=int,
                        help="if support reads greater than this parameter, "
                             "this Allele will be supposed to a peak, default is 3")
    parser.add_argument("-trd", "--threshold", default=0.2, required=False, type=float,
                        help="if msing_score greater than this parameter, "
                             "this sample will be supposed to MSI, default is 0.2")
    parser.add_argument("-mapq", "--mapq", default=20, type=int,
                        help="mapping_quality for read, default is 20")
    parser.add_argument("-dp", "--DP", default=30, type=int,
                        help="depth for all base loci, default is 50")
    parser.add_argument("-flank", action="store", dest="flank", default=150,
                        help="the flank of ms loc")
    parser.add_argument("-fastahack", action="store", dest="fastahack",
                        help="the path of fastahack")
    parser.add_argument("-specbp", action="store", dest="specbp", default=5,
                        help="left or right number of base")
    parser.add_argument("-p", "--patient", type=str,
                        help="output file prefix")
    parser.add_argument("-o", "--outdir", type=str, default=os.getcwd(),help="output file dir, default is current dir")
    parser.add_argument("-mkc", "--make_control", default=False, type=bool, choices=(True, False), help="if you want to make new control file, please run create_baseline.py follow run this script, default is False")
    args = parser.parse_args()
    ###########database import
    os.chdir(args.outdir)
    bam = args.bam
    pex = args.patient
    make_control = args.make_control
    ref = args.ref_genome
    bed = args.bed
    sup_reads = args.sup_reads   ##if support reads greater than this parameter,this Allele will be supposed to a peak, default is 3
    mapq = args.mapq                    ##mapping_quality for read, default is 20
    dp = args.DP                         ###"depth for all base loci, default is 30
    flank = args.flank             ##the flank of ms loc
    specbp = args.specbp             ##left or right number of base
    fastahack = args.fastahack
    control_file = args.control
    bed2 = get_left_right_base(bed, ref, specbp, fastahack)
    result = msi_analysis(bed2, bam, dp, flank, mapq, sup_reads)
    if make_control:
        result.to_csv(pex + "_MSIscore.xls", index=False, sep="\t")
    if not make_control:
        parse_msi(result, control_file, pex, sup_reads)


if __name__ == "__main__":
    main()
