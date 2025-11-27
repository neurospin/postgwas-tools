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
leadSNP.py

Usage example:
python3 leadSNP.py \
  --sumstats most_orig.sumstats \
  --col-a1 A1 --col-a2 A2 \
  --bfile /path/to/filt_imputed_autosomes_maf-0.01 \
  --plink "plink" \
  --out results

This script:
 - runs PLINK clumping for lead (r2=r2_lead) and independent snps (r2=r2_ind),
 - expands locus windows using all SNPs listed in SP2+lead,
 - merges loci within merge-distance (kb) on same chromosome,
 - writes GenomicRiskLoci.txt, leadSNPs.txt, IndSigSNPs.txt, params.config
"""
import argparse
import os
import subprocess
import sys
import pandas as pd
import gc
from postgwas_tools.annot.utils import read_sumstats

def parse_args():
    p = argparse.ArgumentParser(description="Run PLINK clumping and produce FUMA-like outputs with merged locus windows.")
    p.add_argument("--sumstats", required=True, help="GWAS summary stats file (will be reformatted).")
    p.add_argument("--col-a1", default="A1", help="Column name for effect allele / A1 (default: A1)")
    p.add_argument("--col-a2", default="A2", help="Column name for other allele / A2 (default: A2)")
    p.add_argument("--bfile", required=True, help="PLINK --bfile prefix (path to .bed/.bim/.fam without extension)")
    p.add_argument("--plink", default="plink", help="PLINK binary to call (default 'plink')")
    p.add_argument("--r2-lead", type=float, default=0.1, help="r2 for lead SNP clumping (default: 0.1)")
    p.add_argument("--r2-ind", type=float, default=0.6, help="r2 for independent SNPs clumping (default: 0.6)")
    p.add_argument("--clump-p1", type=float, default=5e-8, help="PLINK --clump-p1 (default 5e-8)")
    p.add_argument("--clump-p2", type=float, default=5e-6, help="PLINK --clump-p2 (default 5e-6)")
    p.add_argument("--clump-kb", type=int, default=1000, help="PLINK --clump-kb (default 1000)")
    p.add_argument("--merge-distance", type=int, default=250, help="distance in kb to merge loci (default 250)")
    p.add_argument("--out", required=True, help="Output directory")
    p.add_argument("--keep-intermediate", action="store_true", help="Keep intermediate reformatted sumstats file")
    return p.parse_args()

def write_good_format(df, outpath, p2_threshold):
    # drop NA in essential fields
    sumstats = df.dropna(subset=["SNP", "CHR", "BP", "P"])

    # Apply p2 threshold prefilter (reduce file size for PLINK)
    if p2_threshold is not None:
        sumstats = sumstats[sumstats["P"] <= float(p2_threshold)].copy()

    # Create uniqID here with alphabetically sorted alleles if A1/A2 present
    def make_uid(row):
        a1 = row.get("A1", "NA")
        a2 = row.get("A2", "NA")
        try:
            a1s, a2s = sorted([str(a1), str(a2)])
        except Exception:
            a1s, a2s = str(a1), str(a2)
        return f"{int(row['CHR'])}:{int(row['BP'])}:{a1s}:{a2s}"
    sumstats["uniqID"] = sumstats.apply(make_uid, axis=1)

    # Save and return
    sumstats.to_csv(outpath, sep="\t", index=False)
    return sumstats

def run_plink_clump(plink_bin, bfile_prefix, sumstats_file, out_prefix, p1, p2, r2, kb):
    cmd = [
        plink_bin,
        "--bfile", bfile_prefix,
        "--clump", sumstats_file,
        "--clump-p1", str(p1),
        "--clump-p2", str(p2),
        "--clump-r2", str(r2),
        "--clump-kb", str(kb),
        "--out", out_prefix
    ]
    print("Running PLINK:", " ".join(cmd))
    res = subprocess.run(" ".join(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, shell=True)
    if res.returncode != 0:
        print("PLINK failed. stdout:\n", res.stdout)
        print("PLINK stderr:\n", res.stderr, file=sys.stderr)
        raise RuntimeError("PLINK clumping failed.")
    else:
        # optionally print trimmed stdout / last lines
        if res.stdout:
            print(res.stdout.splitlines()[-5:])
        else:
            print("PLINK completed.")
    return out_prefix + ".clumped"

def read_clumped(path):
    # PLINK clumped files typically have whitespace-separated columns and a header row.
    if not os.path.exists(path):
        print(f"Warning: clumped file {path} not found; returning empty DataFrame")
        return pd.DataFrame()
    try:
        df = pd.read_csv(path, sep=r"\s+", comment="#", engine='python')
        return df
    except Exception as e:
        print(f"Issue when reading plink output {path}: {e}")
        return pd.DataFrame()
    
def sort_by_chr_pos(df):
    """
    Sort numerically by chromosome and position extracted from uniqID
    """
    if df.empty or "uniqID" not in df.columns:
        return df
    df = df.copy()
    df[['chr_sort', 'pos_sort']] = df['uniqID'].str.extract(r'^(\d+):(\d+):')
    df['chr_sort'] = pd.to_numeric(df['chr_sort'], errors='coerce')
    df['pos_sort'] = pd.to_numeric(df['pos_sort'], errors='coerce')
    df = df.sort_values(['chr_sort', 'pos_sort']).reset_index(drop=True)
    df = df.drop(columns=['chr_sort', 'pos_sort'])
    return df

def parse_SP2_field(sp2_field):
    """Parse PLINK SP2 field like 'rs1(1),rs2(1),...' -> ['rs1','rs2',...]"""
    if pd.isna(sp2_field):
        return []
    items = [x.strip() for x in str(sp2_field).split(',') if x.strip()]
    snps = []
    for it in items:
        if it.upper() == "NONE":
            continue
        if '(' in it:
            snp = it.split('(')[0].strip()
        else:
            snp = it
        if snp:
            snps.append(snp)
    return snps

def compute_locus_windows_from_clumps(lead_df, good_sumstats_df):
    """
    From PLINK clump lead (r2=r2_lead), compute locus windows (start/end) using all SNPs
    mentioned in each clump's SP2 plus lead SNP.

    Returns:
      windows_df: DataFrame with columns: CHR, LeadSNP, LeadP, Start, End, nSNPs, SNPs
      lead_to_sp2: dict mapping lead_snp -> list_of_sp2_snps (does NOT filter by independent set here)
    """
    if lead_df.empty:
        return pd.DataFrame(columns=["CHR","LeadSNP","LeadP","Start","End","nSNPs","SNPs"]), {}

    # mapping SNP -> BP and P from good sumstats
    snp_to_bp = dict(zip(good_sumstats_df["SNP"].astype(str), good_sumstats_df["BP"].astype(int)))
    snp_to_p = dict(zip(good_sumstats_df["SNP"].astype(str), good_sumstats_df["P"].astype(float)))

    records = []
    lead_to_sp2 = {}

    for _, row in lead_df.iterrows():
        lead_snp = str(row.get('SNP', row.get('snp', None)))
        if lead_snp is None or lead_snp.upper() == "NONE":
            continue
        chr_val = row.get('CHR', row.get('chr', None))
        if pd.isna(chr_val):
            continue
        try:
            chr_int = int(chr_val)
        except Exception:
            chr_int = int(str(chr_val))

        lead_p = float(row.get('P', row.get('p', float('nan'))))
        sp2_field = row.get('SP2', row.get('sp2', None))
        sp2_snps = parse_SP2_field(sp2_field)
        # make sure lead included in mapping (FUMA considers lead as potential independent)
        combined_snps = [lead_snp] + [s for s in sp2_snps if s != lead_snp]
        # remember mapping for later filtering against independent set
        lead_to_sp2[lead_snp] = combined_snps

        # Collect BP positions for available SNPs
        bp_positions = [snp_to_bp[s] for s in combined_snps if s in snp_to_bp]
        # fallback to BP from row if nothing in good_sumstats
        if not bp_positions:
            bp_from_row = row.get('BP', row.get('bp', None))
            if not pd.isna(bp_from_row):
                try:
                    bp_positions = [int(bp_from_row)]
                except Exception:
                    bp_positions = []
        if not bp_positions:
            # skip locus if no BP found
            continue
        start = min(bp_positions)
        end = max(bp_positions)
        nSNPs = len(set([s for s in combined_snps if s in snp_to_bp]))
        records.append({
            "CHR": chr_int,
            "LeadSNP": lead_snp,
            "LeadP": lead_p,
            "Start": int(start),
            "End": int(end),
            "nSNPs": nSNPs,
            "SNPs": ";".join(sorted(set([s for s in combined_snps if s in snp_to_bp])))
        })

    windows_df = pd.DataFrame(records)
    if windows_df.empty:
        return windows_df, lead_to_sp2
    windows_df = windows_df.sort_values(["CHR", "Start"]).reset_index(drop=True)
    return windows_df, lead_to_sp2

def merge_windows_by_distance(windows_df, merge_kb, lead_to_sp2):
    """
    Merge windows on the same chromosome when the distance between them is <= merge_kb*1000 bp.

    Returns merged list of dicts with keys:
      CHR, Start, End, SNPs (concatenated), LeadSNPs (list), LeadP_min, nSNPs, LeadToSP2 (mapping for the merged region)
    """
    if windows_df.empty:
        return []

    merged = []
    merge_bp = merge_kb * 1000
    windows_df = windows_df.sort_values(["CHR", "Start"]).reset_index(drop=True)

    # initialize current with first row
    first = windows_df.loc[0]
    current = {
        "CHR": int(first["CHR"]),
        "Start": int(first["Start"]),
        "End": int(first["End"]),
        "SNPs": set(str(first["SNPs"]).split(";")) if first["SNPs"] else set(),
        "LeadSNPs": [first["LeadSNP"]],
        "LeadP_min": float(first["LeadP"]),
        "nSNPs": int(first["nSNPs"]),
        "LeadToSP2": { first["LeadSNP"]: lead_to_sp2.get(first["LeadSNP"], []) }
    }

    for i in range(1, len(windows_df)):
        r = windows_df.loc[i]
        chr_i = int(r["CHR"])
        start_i = int(r["Start"])
        end_i = int(r["End"])
        snps_i = set(str(r["SNPs"]).split(";")) if r["SNPs"] else set()
        lead_i = r["LeadSNP"]
        if chr_i == current["CHR"] and start_i - current["End"] <= merge_bp:
            # merge into current
            current["End"] = max(current["End"], end_i)
            current["SNPs"].update(snps_i)
            current["LeadSNPs"].append(lead_i)
            current["LeadP_min"] = min(current["LeadP_min"], float(r["LeadP"]))
            current["nSNPs"] += int(r["nSNPs"])
            # merge lead->sp2 mapping
            current["LeadToSP2"][lead_i] = lead_to_sp2.get(lead_i, [])
        else:
            # finalize current
            current["SNPs"] = ";".join(sorted([s for s in current["SNPs"] if s and s.upper()!="NONE"]))
            merged.append(current)
            # start new current
            current = {
                "CHR": chr_i,
                "Start": start_i,
                "End": end_i,
                "SNPs": snps_i,
                "LeadSNPs": [lead_i],
                "LeadP_min": float(r["LeadP"]),
                "nSNPs": int(r["nSNPs"]),
                "LeadToSP2": { lead_i: lead_to_sp2.get(lead_i, []) }
            }
    # append last
    current["SNPs"] = ";".join(sorted([s for s in current["SNPs"] if s and s.upper()!="NONE"]))
    merged.append(current)
    return merged

# Ensure we have the sort_by_chr_pos helper (use existing one if present)
def sort_by_chr_pos_local(df):
    if df.empty or "uniqID" not in df.columns:
        return df
    d = df.copy()
    d[['chr_sort', 'pos_sort']] = d['uniqID'].str.extract(r'^(\d+):(\d+):')
    d['chr_sort'] = pd.to_numeric(d['chr_sort'], errors='coerce')
    d['pos_sort'] = pd.to_numeric(d['pos_sort'], errors='coerce')
    d = d.sort_values(['chr_sort', 'pos_sort']).reset_index(drop=True)
    d = d.drop(columns=['chr_sort', 'pos_sort'])
    return d

def build_outputs(merged_regions, good_sumstats_df, ind_df):
    """
    Build GenomicRiskLoci.txt, leadSNPs.txt, IndSigSNPs.txt

    Fixes:
     - nSNPs/nGWASSNPs for each independent SNP are computed from its own LD block
       as defined by the r2_ind PLINK clump (ind_df -> SP2).
     - Lead rows are preserved even when they have zero IndSigSNPs (show 0 and empty field).
     - Sorting is done by uniqID (chr:pos:...) and No columns are sequential after sorting.
    """
    gs = good_sumstats_df.copy()
    gs["SNP"] = gs["SNP"].astype(str)
    snp_to_bp = dict(zip(gs["SNP"], gs["BP"]))
    snp_to_p = dict(zip(gs["SNP"], gs["P"]))
    snp_to_uniq = dict(zip(gs["SNP"], gs.get("uniqID", gs["SNP"])))

    # Build independent LD map from ind_df (leader -> [lead + SP2-members])
    ind_ld_map = {}
    independent_snps = set()
    if not ind_df.empty:
        for _, row in ind_df.iterrows():
            lead = str(row["SNP"])
            sp2_members = parse_SP2_field(row.get("SP2", None))
            members = [lead] + [s for s in sp2_members if s != lead]
            ind_ld_map[lead] = members
            independent_snps.add(lead)

    # prefer global sort_by_chr_pos if exists
    try:
        sort_by_chr_pos_fn = sort_by_chr_pos  # noqa: F821
    except Exception:
        sort_by_chr_pos_fn = sort_by_chr_pos_local

    genomic_rows = []
    lead_rows = []
    ind_rows = []

    # For numbering final files (will be assigned after sorting)
    # Build data structures
    for locus_idx, reg in enumerate(merged_regions, start=1):
        chr_ = int(reg["CHR"])
        start = int(reg["Start"])
        end = int(reg["End"])
        # SNPs listed in merged region (from SP2 expansion)
        snps_in_region = [s for s in reg["SNPs"].split(";") if s]
        nSNPs_region = len(snps_in_region)

        # Determine representative lead SNP for locus (min p among GWAS SNPs in region)
        lead_snp = None
        lead_p = float("nan")
        lead_pos = None
        uniqid = None
        best_p = float("inf")
        for s in snps_in_region:
            if s in snp_to_p:
                pval = float(snp_to_p[s])
                if pval < best_p:
                    best_p = pval
                    lead_snp = s
                    lead_p = pval
                    lead_pos = snp_to_bp.get(s)
                    uniqid = snp_to_uniq.get(s)
        if lead_snp is None:
            # fallback to first lead if no p found
            if reg.get("LeadSNPs"):
                lead_snp = reg["LeadSNPs"][0]
                lead_p = float("nan")
                lead_pos = None
                uniqid = f"{chr_}:{start}:NA:NA"
            else:
                lead_snp = None
                uniqid = f"{chr_}:{start}:NA:NA"

        # Build per-lead IndSigSNPs using LeadToSP2 mapping (keeps possibility of zero hits)
        region_ind_snps = set()
        lead_to_sp2 = reg.get("LeadToSP2", {})
        per_lead_ind_map = {}
        for lead in reg.get("LeadSNPs", []):
            # get SP2 advertised for this lead (from merged mapping)
            candidates = lead_to_sp2.get(lead, [])
            if lead not in candidates:
                candidates = [lead] + candidates
            # keep only those that are independent (i.e. appear in ind_ld_map keys)
            ind_for_lead = [s for s in candidates if s in independent_snps]
            per_lead_ind_map[lead] = ind_for_lead
            region_ind_snps.update(ind_for_lead)

        nInd = len(region_ind_snps)

        # GenomicRiskLoci row
        genomic_rows.append({
            "GenomicLocus": locus_idx,
            "uniqID": uniqid,
            "rsID": lead_snp,
            "chr": chr_,
            "pos": lead_pos,
            "p": lead_p,
            "start": start,
            "end": end,
            "nSNPs": nSNPs_region,
            "nGWASSNPs": nSNPs_region,
            "nIndSigSNPs": nInd,
            "IndSigSNPs": ";".join(sorted(region_ind_snps)),
            "nLeadSNPs": len(reg.get("LeadSNPs", [])),
            "LeadSNPs": ";".join(reg.get("LeadSNPs", []))
        })

        # LeadSNP rows (preserve even when ind_for_lead is empty)
        for lead in reg.get("LeadSNPs", []):
            lead_bp = snp_to_bp.get(lead, None)
            lead_pval = snp_to_p.get(lead, float("nan"))
            lead_uniq = snp_to_uniq.get(lead, f"{chr_}:{lead_bp}:NA:NA")
            ind_for_lead = per_lead_ind_map.get(lead, [])
            # record lead row (nIndSigSNPs may be 0; IndSigSNPs empty string)
            lead_rows.append({
                "No": None,               # will set after sorting
                "GenomicLocus": locus_idx,
                "uniqID": lead_uniq,
                "rsID": lead,
                "chr": chr_,
                "pos": lead_bp,
                "p": lead_pval,
                "nIndSigSNPs": len(ind_for_lead),
                "IndSigSNPs": ";".join(ind_for_lead)
            })

        # Independent SNP rows: compute per-independent-SNP LD counts from ind_ld_map
        for s in sorted(region_ind_snps):
            ld_members = ind_ld_map.get(s, [s])  # members listed by clump_ind (lead + SP2)
            nSNPs_ld = len(ld_members)
            # count how many of those members are present in the GWAS file (good_sumstats_df)
            nGWASSNPs_ld = sum(1 for x in ld_members if x in snp_to_uniq)
            ind_rows.append({
                "No": None,   # assigned later after sorting
                "GenomicLocus": locus_idx,
                "uniqID": snp_to_uniq.get(s, f"{chr_}:{snp_to_bp.get(s)}:NA:NA"),
                "rsID": s,
                "chr": chr_,
                "pos": snp_to_bp.get(s),
                "p": snp_to_p.get(s),
                "nSNPs": nSNPs_ld,
                "nGWASSNPs": nGWASSNPs_ld
            })

    # Convert to DataFrames
    genomic_df = pd.DataFrame(genomic_rows)
    lead_df = pd.DataFrame(lead_rows)
    ind_df_out = pd.DataFrame(ind_rows)

    # Sort using uniqID-based sort and then set No sequentially
    genomic_df = sort_by_chr_pos_fn(genomic_df) if not genomic_df.empty else genomic_df
    lead_df = sort_by_chr_pos_fn(lead_df) if not lead_df.empty else lead_df
    ind_df_out = sort_by_chr_pos_fn(ind_df_out) if not ind_df_out.empty else ind_df_out

    # Assign No as 1-based row number (after sorting)
    if not lead_df.empty:
        lead_df["No"] = range(1, len(lead_df) + 1)
    if not ind_df_out.empty:
        ind_df_out["No"] =  range(1, len(ind_df_out) + 1)

    # Ensure genomic_df columns ordering matches spec
    if not genomic_df.empty:
        genomic_df = genomic_df[[
            "GenomicLocus","uniqID","rsID","chr","pos","p","start","end",
            "nSNPs","nGWASSNPs","nIndSigSNPs","IndSigSNPs","nLeadSNPs","LeadSNPs"
        ]]

    # leadSNPs must have the exact columns requested: No, GenomicLocus, uniqID, rsID, chr, pos, p, nIndSigSNPs, IndSigSNPs
    if not lead_df.empty:
        lead_df = lead_df[["No","GenomicLocus","uniqID","rsID","chr","pos","p","nIndSigSNPs","IndSigSNPs"]]

    # IndSigSNPs columns ordering
    if not ind_df_out.empty:
        ind_df_out = ind_df_out[["No","GenomicLocus","uniqID","rsID","chr","pos","p","nSNPs","nGWASSNPs"]]

    return genomic_df, lead_df, ind_df_out

def write_results(outdir, genomic_df, lead_df, ind_df_out):
    # Write files
    genomic_path = os.path.join(outdir, "GenomicRiskLoci.txt")
    lead_path = os.path.join(outdir, "leadSNPs.txt")
    ind_path = os.path.join(outdir, "IndSigSNPs.txt")

    # Always write header — even when dataframe is empty
    genomic_df.to_csv(genomic_path, sep="\t", index=False)
    lead_df.to_csv(lead_path, sep="\t", index=False)
    ind_df_out.to_csv(ind_path, sep="\t", index=False)

    print(f"GenomicRiskLoci saved to: {genomic_path}")
    print(f"Lead SNP table saved to: {lead_path}")
    print(f"Independent SNP table saved to: {ind_path}")

def remove_files(outdir, okstats_path):
        # list of intermediate files to delete if they exist
        intermediate_files = [
            okstats_path,
            os.path.join(outdir, "clump_ind.clumped"),
            os.path.join(outdir, "clump_lead.clumped"),
            os.path.join(outdir, "clump_ind.log"),
            os.path.join(outdir, "clump_lead.log"),
        ]

        for f in intermediate_files:
            try:
                if os.path.exists(f):
                    os.remove(f)
            except Exception:
                pass

def main():
    args = parse_args()
    outdir = args.out
    os.makedirs(outdir, exist_ok=True)

    # Read and standardize input sumstats (this uses your read_sumstats utils)
    print("Reading summary statistics...")
    df = read_sumstats(args.sumstats, A1=args.col_a1, A2=args.col_a2)
    okstats_path = os.path.join(outdir, "good_format_gwas.txt")

    # prefilter by clump_p2 to reduce file size/time/memory
    print("Reformatting and prefiltering summary stats...")
    good_df = write_good_format(df, okstats_path, p2_threshold=args.clump_p2)
    print(f"Reformatted (prefiltered p<={args.clump_p2}) sumstats written to {okstats_path}; {len(good_df)} rows")

    # free original large dataframe
    try:
        del df
        gc.collect()
    except Exception:
        pass

    # Run PLINK clumping for lead (r2_lead) and independents (r2_ind)
    clump_lead_prefix = os.path.join(outdir, "clump_lead")
    clump_ind_prefix = os.path.join(outdir, "clump_ind")

    clumped_lead_path = run_plink_clump(args.plink,
                                        args.bfile,
                                        okstats_path,
                                        clump_lead_prefix,
                                        args.clump_p1,
                                        args.clump_p2,
                                        args.r2_lead,
                                        args.clump_kb)
    clumped_ind_path = run_plink_clump(args.plink,
                                       args.bfile,
                                       okstats_path,
                                       clump_ind_prefix,
                                       args.clump_p1,
                                       args.clump_p2,
                                       args.r2_ind,
                                       args.clump_kb)

    # Read PLINK outputs
    lead_df = read_clumped(clumped_lead_path)
    ind_df = read_clumped(clumped_ind_path)

    # Normalize column names if needed (p -> P)
    for d in (lead_df, ind_df):
        if not d.empty and 'P' not in d.columns and 'p' in d.columns:
            d.rename(columns={'p': 'P'}, inplace=True)

    # Compute windows from r2_lead clumps using all SNPs listed (SP2)
    print("Computing locus windows from r2_lead clumps (expanding with SP2 SNPs)...")
    windows_df, lead_to_sp2 = compute_locus_windows_from_clumps(lead_df, good_df)
    if windows_df.empty:
        print("No locus windows computed (no clumps or missing data). Exiting.")
        # still write params.config
        config_path = os.path.join(outdir, "params.config")
        with open(config_path, "w") as f:
            for k, v in vars(args).items():
                f.write(f"{k}={v}\n")
        print(f"Parameters saved to: {config_path}")
        print("Remove intermediate files.")
        remove_files(outdir, okstats_path)
        genomic_columns = [
        "GenomicLocus","uniqID","rsID","chr","pos","p","start","end",
        "nSNPs","nGWASSNPs","nIndSigSNPs","IndSigSNPs","nLeadSNPs","LeadSNPs"]
        genomic_df = pd.DataFrame(columns=genomic_columns)
        lead_columns= ["No","GenomicLocus","uniqID","rsID","chr","pos","p",
                       "nIndSigSNPs","IndSigSNPs"]
        lead_df = pd.DataFrame(columns=lead_columns)
        ind_columns=["No","GenomicLocus","uniqID","rsID","chr","pos","p",
                     "nSNPs","nGWASSNPs"]
        ind_df = pd.DataFrame(columns=ind_columns)
        write_results(outdir, genomic_df, lead_df, ind_df)
        return

    # Merge windows by distance (kb)
    print(f"Merging windows closer than {args.merge_distance} kb on same chromosome...")
    merged_regions = merge_windows_by_distance(windows_df, args.merge_distance, lead_to_sp2)

    # Build and write outputs: GenomicRiskLoci.txt, leadSNPs.txt, IndSigSNPs.txt
    print("Building final tables...")
    genomic_df, lead_df, ind_df_out = build_outputs(merged_regions, good_df, ind_df)
    write_results(outdir, genomic_df, lead_df, ind_df_out)
    # Write params.config
    config_path = os.path.join(outdir, "params.config")
    with open(config_path, "w") as f:
        for k, v in vars(args).items():
            f.write(f"{k}={v}\n")
    print(f"Parameters saved to: {config_path}")

    if not args.keep_intermediate:
        remove_files(outdir, okstats_path)

if __name__ == "__main__":
    main()
