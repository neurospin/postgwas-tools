# coding: utf-8
##########################################################################
# NSAp - Copyright (C) CEA, 2023
# Distributed under the terms of the CeCILL-B license.
# Antoine Dufournet
##########################################################################

import pandas as pd
import matplotlib.colors as mcolors
import glob
import os

def find_files(paths):
    """
    Find all files from a list of paths or wildcard patterns.

    Parameters
    ----------
    paths : list of str
        List of file paths, directory paths, or glob patterns.

    Returns
    -------
    list of str
        Sorted list of matching file paths.
    """
    all_files = []

    for path in paths:
        # If it's a directory, take all files inside
        if os.path.isdir(path):
            matched = glob.glob(os.path.join(path, "**"), recursive=True)
            matched = [f for f in matched if os.path.isfile(f)]
        # If it contains a wildcard pattern like *.sumstats
        elif "*" in path or "?" in path or "[" in path:
            matched = glob.glob(path, recursive=True)
        # If it's a direct file path
        elif os.path.isfile(path):
            matched = [path]
        else:
            print(f"Warning: No match found for '{path}'")
            matched = []

        all_files.extend(matched)

    # Deduplicate and sort
    all_files = sorted(list(set(all_files)))
    return all_files

def adjust_color_brightness(color, factor):
    """
    Adjusts brightness of a color by blending with white (factor > 1)
    or black (factor < 1).
    """
    color = mcolors.to_rgb(color)
    return tuple(min(1, max(0, c * factor)) for c in color)


def read_sumstats(file_path, A1=None, A2=None):
    """
    Read a GWAS summary statistics file and standardize column names.

    Automatically detects separator (tab, comma, whitespace) and supports both
    plain text and gzipped (.gz) files. Maps columns:
      - 'chr', '#chr', '#chrom', 'chrom', 'chromosome' -> 'CHR'
      - 'bp', 'pos', 'posgrch37', 'position-> 'BP'
      - 'snp', 'id', 'markername', 'rs', 'rs_number', 'snpid', 'rsid' -> 'SNP'
      - 'p', 'pvalue', 'p_value', 'pval', 'p_val', 'gc_pvalue', 'p-value' -> 'P'

    Parameters
    ----------
    file_path : str
        Path to the summary statistics file (can be .txt, .csv, .tsv, or .gz).

    Returns
    -------
    pd.DataFrame
        Standardized dataframe with at least 'CHR', 'BP', 'SNP', and 'P' columns.
    """
    # --- Detect compression ---
    compression = 'gzip' if file_path.endswith('.gz') else None

    # --- Try multiple separators in order of likelihood ---
    seps = ['\t', ',', r'\s+']
    last_error = None
    df = None
    for sep in seps:
        try:
            df = pd.read_csv(file_path, sep=sep, compression=compression, engine='python')
            # Basic sanity check: at least a few columns
            if df.shape[1] >= 3:
                break
        except Exception as e:
            last_error = e
            continue

    if df is None:
        raise ValueError(f"Could not read file {file_path}. Last error: {last_error}")

    # --- Normalize column names ---
    df.columns = [col.strip().lower() for col in df.columns]

    # --- Standardize known column names ---
    col_map = {}
    for col in df.columns:
        if col in ['chr', '#chr', '#chrom', 'chrom', 'chromosome']:
            col_map[col] = 'CHR'
        elif col in ['bp', 'pos', 'posgrch37', 'position']:
            col_map[col] = 'BP'
        elif col in ['snp', 'rsid', 'markername', 'rs', 'id', 'rs_number', 'snpid']:
            col_map[col] = 'SNP'
        elif col in ['p', 'pvalue', 'p_value', 'pval', 'p_val', 'gc_pvalue', 'p-value']:
            col_map[col] = 'P'
    if A1 and A1.strip().lower() in df.columns:
        col_map[A1.strip().lower()] = 'A1'
    if A2 and A2.strip().lower() in df.columns:
        col_map[A2.strip().lower()] = 'A2'

    df = df.rename(columns=col_map)

    # --- Validation ---
    list_to_keep = ['CHR', 'BP', 'SNP', 'P']
    if A1:
        list_to_keep+=['A1']
    if A2:
        list_to_keep+=['A2']

    missing = [c for c in list_to_keep if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns {missing} in {file_path}")

    # --- Keep only the main columns ---
    keep_cols = [c for c in list_to_keep if c in df.columns]
    df = df[keep_cols]

    return df

