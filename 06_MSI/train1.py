import matplotlib
matplotlib.use('Agg')
import pandas as pd
import pysam
import math
from collections import defaultdict
import argparse as apa
import seaborn as sns
import matplotlib.pyplot as plt
import traceback
import re
import sys
import os
import warnings
warnings.filterwarnings("ignore")
sns.set(color_codes=True)


def get_left_right_base(fastahack, bed, ref, specbp):
#    fastahack = "/data1/workdir/wangce/software/fastahack/fastahack"
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


def msi_analysis(bed2, bam, cutoff, dp, flank, mapq, sup_reads):
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


def plot(pex, result):
    result = result.loc[result["qcs"] == "pass"]
    num = math.ceil(len(result) / 3)
    numi = 1
    for ms, info in zip(result["MSID"], result["IndelLength:AlleleFraction:SupportingCalls"]):
        x = []
        ax = plt.subplot(num, 3, numi)
        numi = numi + 1
        xticks = []
        for peak in info.split(" "):
            x += [int(peak.split(":")[0])] * int(peak.split(":")[2])
            xticks.append(int(peak.split(":")[0]))
        plt.title(ms, fontsize=4)
        sns.kdeplot(x, bw=.15, ax=ax)
        plt.xticks(fontsize=4)
        ax.set_yticklabels("", size=0)
        ax.tick_params(width=0.1, length=0.01, pad=0.01)
        plt.tight_layout(pad=0.05)
    plt.savefig(pex + "_MSI.pdf")


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


def parse_msi(result, qc_data, multiplier, threshold, pex, sup_reads):
    control = pd.read_table(qc_data)
    result = normalization(result, control, sup_reads)
    control["value"] = control["Average_Number_Peaks"] + multiplier * control["Standard_Deviation"]
    control = {msi_loci: value for msi_loci, value in zip(control["MSID"], control["value"])}
    msi_status = []
    for qc, peak, name in zip(result["qcs"], result["Normalized_Number_of_Peaks"], result["MSID"]):
        if qc == "fail" or name not in control.keys():
            msi_status.append("-")
        else:
            if peak >= control[name]:
                msi_status.append("unstable")
            else:
                msi_status.append("stable")
    result["msi_status"] = msi_status
    result.to_csv(pex + "_MSIscore.xls", index=False, sep="\t")
    result = result.loc[result["qcs"] == "pass"]
    jinbao = ["BAT25", "BAT26", "NR24", "NR21", "MONO27", "D5S346", "D2S123", "NR27"]
    jinbao_result = result.loc[result["MSID"].isin(jinbao)]
    f = open(pex + ".MSI.xls", "w")
    if len(result) == 0:
        msing_score = 0
    else:
        msing_score = len(result.loc[result["msi_status"] == "unstable"]) / len(result)

    if len(jinbao_result.loc[jinbao_result["msi_status"] == "unstable"]) >= 2 or msing_score >= threshold:
        f.write("MSI-H")
    elif len(jinbao_result.loc[jinbao_result["msi_status"] == "unstable"]) == 1 or 0.1 < msing_score < threshold:
        f.write("MSI-L")
    else:
        f.write("MSS")
    return result


def main(args):
    os.chdir(args.outdir)
    fastahack = args.fastahack
    ref = args.ref_genome
    bam = args.bam
    bed = args.bed
    pex = args.patient
    cutoff = args.cutoff
    sup_reads = args.sup_reads
    qc_data = args.qc_data
    multiplier = args.multiplier
    threshold = args.threshold
    dp = args.DP
    mapq = args.mapq
    flank = args.flank
    specbp = args.specbp
    bed2 = get_left_right_base(fastahack, bed, ref, specbp)
    result = msi_analysis(bed2, bam, cutoff, dp, flank, mapq, sup_reads)
    parse_msi(result, qc_data, multiplier, threshold, pex, sup_reads)


if __name__ == "__main__":
    usage = "\n\tMSI analysis"
    parser = apa.ArgumentParser(prog="convert")
    parser.add_argument("-fastahack", "--fastahack", required=True, type=str,
                        help="the path of fastahack")
    parser.add_argument("-rf", "--ref_genome", required=True, type=str,
                        help="referance genome, such as ucsc.hg19.fasta")
    parser.add_argument('-bam', '--bam', required=True, type=str,
                        help='bam file')
    parser.add_argument('-bed', '--bed', required=True, type=str,
                        help='panel8 str bed file')
    parser.add_argument('-qc', "--qc_data", required=False, type=str,
                        help='qc data file')
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
    parser.add_argument("-specbp", action="store", dest="specbp", default=5,
                        help="left or right number of base")
    parser.add_argument("-p", "--patient", type=str,
                        help="output file prefix")
    parser.add_argument("-o", "--outdir", type=str, default=os.getcwd(),
                        help="output file dir, default is current dir")
    args = parser.parse_args()
    try:
        main(args)
        sys.exit(0)
    except Exception as e:
        traceback.print_exc()
        sys.exit(1)
