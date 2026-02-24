#!/usr/bin/env python3
"""
summarize_replication_nature.py

Aggregate region-wise replication results. 

Usage: python3 summarize_replication.py \ 
        --discovery /path/to/UKB/Champollion_V1_32 \ 
        --replication /path/to/ABCD/Champollion_V1_32 \ 
        --out /path/to/final_replication_summary.tsv

Generate publication-quality replication figure (Population: All only).
Outputs:
- region-level TSV summary
- Nature-style horizontal dot plot (PDF vector format)
"""

import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt


def summarize_region(lead_file, repl_file):
    """Return replication statistics for one region."""
    lead = pd.read_csv(lead_file, sep="\t", dtype=str)
    repl = pd.read_csv(repl_file, sep="\t", dtype=str)

    rs_col_lead = next((c for c in lead.columns if c.lower() in ("rsid","rs_id","snp","variant","markername")), None)
    rs_col_repl = next((c for c in repl.columns if c.lower() in ("rsid","rs_id","snp","variant","markername")), None)
    p_col_repl = next((c for c in repl.columns if c.lower() in ("pval","p_value","p-value","p")), None)

    if not rs_col_lead or not rs_col_repl or not p_col_repl:
        raise ValueError(f"Missing rsID or PVAL columns in {lead_file} or {repl_file}")

    lead = lead[[rs_col_lead]].drop_duplicates()
    repl = repl[[rs_col_repl, p_col_repl]].rename(columns={p_col_repl: "p_rep"})

    merged = pd.merge(lead, repl, left_on=rs_col_lead, right_on=rs_col_repl, how="inner")

    n_tested = len(merged)
    if n_tested == 0:
        return 0, 0, 0

    merged["p_rep"] = pd.to_numeric(merged["p_rep"], errors="coerce")
    threshold = 0.05 
    n_rep = (merged["p_rep"] < threshold).sum()
    perc = (n_rep / n_tested) * 100

    return n_tested, n_rep, perc


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--discovery", required=True)
    parser.add_argument("--replication", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    results = []

    for root, _, files in os.walk(args.replication):
        if "replication_matches.txt" not in files:
            continue

        # Only keep "All" population
        if "White" in root or "white" in root:
            continue

        repl_file = os.path.join(root, "replication_matches.txt")
        region_path = "/".join(os.path.relpath(repl_file, args.replication).split("/")[:2])

        lead_file = os.path.join(
            args.discovery,
            region_path,
            "32PCs/White/FUMA/leadSNPs.txt"
        )

        if not os.path.exists(lead_file):
            continue

        try:
            n_tested, n_rep, perc = summarize_region(lead_file, repl_file)

            results.append({
                "region": region_path.split("/name")[0],
                "N_tested": n_tested,
                "N_replicated": n_rep,
                "Replication_%": perc
            })

            print(f"{region_path}: {n_rep}/{n_tested} ({perc:.1f}%)")

        except Exception as e:
            print(f"{region_path} failed: {e}")

    if not results:
        print("No results found.")
        return

    df = pd.DataFrame(results)

    # Sort by replication percentage
    df = df.sort_values(by="Replication_%", ascending=True)

    # Save TSV
    tsv_out = args.out.replace(".pdf", ".tsv")
    df.to_csv(tsv_out, sep="\t", index=False)
    print(f"Saved summary: {tsv_out}")

    # -------------------------
    # Publication-quality plot
    # -------------------------

    plt.rcParams.update({
        "font.size": 8,
        "axes.linewidth": 0.8,
        "pdf.fonttype": 42,
        "ps.fonttype": 42
    })

    fig, ax = plt.subplots(figsize=(3.5, 6))  # Nature column width ≈ 3.5 inches

    ax.scatter(
        df["Replication_%"],
        df["region"],
        s=20,
        color="#0072B2"
    )

    ax.set_xlabel("Replication rate (%)")
    ax.set_ylabel("")

    ax.set_xlim(0, 100)

    # Remove top/right spines
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Subtle x grid only
    ax.grid(axis="x", linestyle="--", linewidth=0.5, alpha=0.3)

    plt.tight_layout()

    plt.savefig(args.out, format="pdf", bbox_inches="tight")
    print(f"Saved Nature-style figure: {args.out}")


if __name__ == "__main__":
    main()