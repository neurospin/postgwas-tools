#!/usr/bin/env python3
##########################################################################
# NSAp - Copyright (C) CEA, 2023
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
# Antoine Dufournet
##########################################################################
"""
replication_by_rsid.py

Usage:
    python replication_by_rsid.py \
        --lead-snps path/to/leadSNPs.txt \
        --sumstats path/to/ABCD_sumstats.txt \
        --out path/to/replication_matches.txt

Description:
    For each lead SNP (from FUMA), find its corresponding line in replication
    summary stats (ABCD) by rsID.
    - If the lead SNP is missing, try each independent SNP listed in 'IndSigSNPs'.
    - Only the *first found* SNP per lead is written to the output.
    - A detailed .log file is created describing which were found,
      replaced by fallback SNPs, or missing entirely.
"""

import argparse
import gzip
import os

def open_maybe_gz(path):
    """Open text or gzipped file automatically."""
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path, 'r')


def read_lead_rsids(path):
    """Read lead SNPs and their fallback independent SNPs."""
    lead_to_all = {}
    with open_maybe_gz(path) as f:
        header = f.readline().strip().split('\t')
        # Identify columns
        try:
            rs_col = [i for i, h in enumerate(header) if h.lower() in ('rsid', 'rs_id', 'snp', 'variant', 'markername')][0]
        except IndexError:
            raise ValueError(f"Could not find rsID column in {path}. Header: {header}")

        indsig_col = None
        for i, h in enumerate(header):
            if h == 'IndSigSNPs':
                indsig_col = i
                break

        for line in f:
            if not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) <= rs_col:
                continue
            lead = parts[rs_col].strip()
            proxies = []
            if indsig_col is not None and len(parts) > indsig_col:
                proxies = [s.strip() for s in parts[indsig_col].split(';') if s.strip() and s.strip() != lead]
            lead_to_all[lead] = [lead] + proxies

    print(f"Loaded {len(lead_to_all)} lead SNPs (+ fallback independent SNPs) from {path}")
    return lead_to_all


def index_replication_file(sumstats_path):
    """Read replication file once and build a dictionary of rsID -> line."""
    print(f"Indexing replication file {sumstats_path}...")
    rep_dict = {}
    with open_maybe_gz(sumstats_path) as f:
        header = f.readline().rstrip('\n')
        sep = '\t' if '\t' in header else None
        header_parts = header.split(sep or None)
        rs_col = None
        for i, h in enumerate(header_parts):
            if h.lower() in ('rsid', 'rs_id', 'snp', 'variant', 'markername'):
                rs_col = i
                break
        if rs_col is None:
            raise ValueError(f"Could not find rsID column in replication file. Header: {header_parts}")

        for line in f:
            if not line.strip():
                continue
            parts = line.strip().split(sep)
            if len(parts) <= rs_col:
                continue
            rs = parts[rs_col].strip()
            rep_dict[rs] = line.strip()

    print(f"Indexed {len(rep_dict):,} SNPs from replication file.")
    return rep_dict, header


def match_by_lead(lead_to_all, rep_dict, header, out_path):
    """For each lead SNP, write only one best match (lead or fallback)."""
    log_path = os.path.splitext(out_path)[0] + ".log"

    n_found = 0
    found_direct = []
    replaced = {}
    not_found = []

    with open(out_path, 'w') as wf, open(log_path, 'w') as log:
        wf.write(header + '\n')

        for lead, snp_list in lead_to_all.items():
            found_snp = None
            for snp in snp_list:
                if snp in rep_dict:
                    found_snp = snp
                    break

            if found_snp:
                wf.write(rep_dict[found_snp] + '\n')
                n_found += 1
                if found_snp == lead:
                    found_direct.append(lead)
                else:
                    replaced[lead] = found_snp
            else:
                not_found.append(lead)

        # Logging summary
        total = len(lead_to_all)
        covered = len(found_direct) + len(replaced)
        log.write(f"Total lead SNPs: {total}\n")
        log.write(f"Total found (lead or fallback): {covered}\n")
        log.write(f"Directly found lead SNPs: {len(found_direct)}\n")
        log.write(f"Replaced by fallback SNPs: {len(replaced)}\n")
        log.write(f"Missing lead SNPs: {len(not_found)}\n")
        log.write(f"Coverage: {covered/total*100:.2f}%\n\n")

        log.write("=== Lead SNPs replaced by fallback SNPs ===\n")
        if replaced:
            for lead, fb in replaced.items():
                log.write(f"{lead} → {fb}\n")
        else:
            log.write("None\n")
        log.write("\n")

        log.write("=== Lead SNPs not found ===\n")
        if not_found:
            for lead in not_found:
                log.write(f"{lead}\n")
        else:
            log.write("None\n")

    print(f"Wrote {n_found} matched lines to {out_path}")
    print(f"Detailed log saved to {log_path}")
    print(f"Found {len(found_direct)} direct, {len(replaced)} replaced, {len(not_found)} missing.")
    return {
        "found_direct": found_direct,
        "replaced": replaced,
        "not_found": not_found,
        "log": log_path
    }


def main():
    parser = argparse.ArgumentParser(description="Match lead SNPs by rsID, with fallback and logging (1 match per lead).")
    parser.add_argument("--lead-snps", required=True, help="Lead SNPs file (FUMA format with rsID and IndSigSNPs columns).")
    parser.add_argument("--sumstats", required=True, help="Replication summary stats (can be hg37 or hg38). Can be .gz.")
    parser.add_argument("--out", required=True, help="Output file for matched replication SNPs.")
    args = parser.parse_args()

    lead_to_all = read_lead_rsids(args.lead_snps)
    rep_dict, header = index_replication_file(args.sumstats)
    match_by_lead(lead_to_all, rep_dict, header, args.out)


if __name__ == "__main__":
    main()
