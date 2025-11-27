#!/usr/bin/env python3
##########################################################################
# NSAp - Copyright (C) CEA, 2023
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
# Antoine Dufournet
##########################################################################

import argparse
import os
import numpy as np
import pandas as pd
from postgwas_tools.annot.utils import find_files

def parse_args():
    parser = argparse.ArgumentParser(description="Calculate a dice score from different locus windows (for instance in GenomicRiskLoci.txt).")
    parser.add_argument( "-p", "--paths", nargs="+", type=str, required=True,
        help="List of paths, file names, or wildcard patterns")
    parser.add_argument("--out", required=True, help="Output directory")
    parser.add_argument("--names", nargs="+", help="Names corresponding to input files")
    parser.add_argument("--position", type=int, help="Position of the names corresponding to input files (from right to left when splitting on /)")
    parser.add_argument("--verbose", action="store_true")
    return parser.parse_args()

def read_genomic_loci(genomic_loci_path):
    genomic_loci_df = pd.read_csv(genomic_loci_path, sep="\t")
    genomic_loci_df = genomic_loci_df.rename(columns=str.lower)
    genomic_loci_df = genomic_loci_df.rename(columns={
            'chr': 'CHR',
            'start': 'START',
            'end': 'END'
        })
    return genomic_loci_df


def DSC(df1, df2):
    intersection = 0
    for chromi in df1['CHR'].unique():
        tmp_df1 = df1[df1['CHR']==chromi]
        tmp_df2 = df2[df2['CHR']==chromi]
        for i in range(len(tmp_df1)):
            start = tmp_df1.iloc[i]['START']
            end = tmp_df1.iloc[i]['END']
            for j in range(len(tmp_df2)):
                if not (tmp_df2.iloc[j]['END'] < start or tmp_df2.iloc[j]['START'] > end):
                    intersection+=1
    dsc = 2*intersection/(len(df1)+len(df2))        
    return dsc

def main():
    args = parse_args()
    outdir = args.out
    files = find_files(args.paths)
    n = len(files)
    position_name = args.position

    if args.verbose:
        print(f"Found {n} files")

    if args.names:
        if len(args.names) != len(files):
            raise ValueError("Number of names must match number of files")
        names = args.names
    else:
        names = [f.split('/')[-position_name] for f in files]

    # Initialize result matrix with NaNs
    dsc_matrix = pd.DataFrame(np.nan, index=names, columns=names)
    # Outer loop: load left once, then iterate rights
    for i, left_path in enumerate(files):
        if args.verbose:
            print(f"Loading left ({i+1}/{n}): {left_path}")
        left_df = read_genomic_loci(left_path)
        for j in range(i + 1, n):
            right_path = files[j]
            if args.verbose:
                print(f"  Loading right ({j+1}/{n}): {right_path}")
            right_df = read_genomic_loci(right_path)

             # Fill symmetric positions (i,j) and (j,i)
            dsc = DSC(left_df, right_df)
            dsc_matrix.iat[i, j] = dsc
            dsc_matrix.iat[j, i] = dsc

            del right_df
        del left_df
    for k in range(n):
        dsc_matrix.iat[k, k] = 1.0

    os.makedirs(outdir, exist_ok=True)
    dsc_matrix.to_csv(f"{outdir}/dsc_matrix.csv")
    if args.verbose:
        print(f"Saved dsc matrix to {outdir}/dsc_matrix.csv")


if __name__ == "__main__":
    main()
