import pandas as pd
from glob import glob
import numpy as np
from scipy.stats import ranksums
from statsmodels.stats import multitest as mul
import argparse as apa


def make_bed(qc_data, positive_dir, negative_dir, bed):
    bed = pd.read_table(bed)
    target_loc = pd.read_table(qc_data)
    target_loc = target_loc.loc[target_loc["Average_Total_Reads"] >= 30]
    target_loc = target_loc.loc[target_loc["Average_Number_Peaks"] >= 1.5]
    positive_samples = glob("{0}/*_MSIscore.xls".format(positive_dir))
    negative_samples = glob("{0}/*_MSIscore.xls".format(negative_dir))

    positive_s = []
    for s in positive_samples:
        data = pd.read_table(s)
        data = data.loc[data["MSID"].isin(target_loc["MSID"].tolist())]
        data = data[["MSID", "Normalized_Number_of_Peaks"]]
        data.columns = ["MSID", s.split("/")[-1].split("_")[0] + "_positive"]
        target_loc = pd.merge(target_loc, data, on="MSID", how="inner")
        positive_s.append(s.split("/")[-1].split("_")[0] + "_positive")

    negative_s = []
    for s in negative_samples:
        data = pd.read_table(s)
        data = data.loc[data["MSID"].isin(target_loc["MSID"].tolist())]
        data = data[["MSID", "Normalized_Number_of_Peaks"]]
        data.columns = ["MSID", s.split("/")[-1].split("_")[0] + "_negative"]
        target_loc = pd.merge(target_loc, data, on="MSID", how="inner")
        negative_s.append(s.split("/")[-1].split("_")[0] + "_negative")


    target_loc["pval"] = [ranksums(i, ii).pvalue for i, ii in zip(target_loc[positive_s].as_matrix(), target_loc[negative_s].as_matrix())]
    fdr = target_loc["pval"]
    reject, pvals_corrected = mul.fdrcorrection(fdr)
    target_loc['FDR_bh'] = pvals_corrected
    target_loc = target_loc.loc[target_loc["pval"] <= 0.01]
    target_loc.to_csv("peaks.txt", sep="\t", index=False)
    bed = bed.loc[bed["MSID"].isin(target_loc["MSID"].tolist())]
    bed.to_csv("bed.txt", sep="\t", index=False)


if __name__ == "__main__":
    usage = "\n\tchoose msi loc"
    parser = apa.ArgumentParser(prog="convert")
    parser.add_argument('-qc', "--qc_data", required=False, type=str,
                        help='qc data file')
    parser.add_argument("-bed", "--bed", type=str, required=True, help="your panel bed")
    parser.add_argument('-positive_dir', '--positive_dir', required=True, type=str,
                        help='positive peaks dir')
    parser.add_argument('-negative_dir', '--negative_dir', required=True, type=str,
                        help='negative peaks dir')
    args = parser.parse_args()
    make_bed(args.qc_data, args.positive_dir, args.negative_dir, args.bed)
