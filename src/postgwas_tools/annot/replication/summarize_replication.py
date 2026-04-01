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
    n_rep_bonf = (merged["p_rep"] < threshold/n_tested).sum()
    perc = (n_rep / n_tested) * 100
    perc_bonf = (n_rep_bonf / n_tested) * 100

    return n_tested, n_rep, perc, n_rep_bonf, perc_bonf


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--discovery", required=True)
    parser.add_argument("--replication", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    results = []

    for root, _, files in os.walk(args.replication):
        if "replication_matches_8e-10.txt" not in files:
            continue

        # Only keep "All" population
        if "White" in root or "white" in root:
            continue

        repl_file = os.path.join(root, "replication_matches_8e-10.txt")
        region_path = "/".join(os.path.relpath(repl_file, args.replication).split("/")[:2])

        lead_file = os.path.join(
            args.discovery,
            region_path,
            "32PCs/White/FUMA8e-10most/leadSNPs.txt"
        )

        if not os.path.exists(lead_file):
            continue

        try:
            n_tested, n_rep, perc, n_rep_bonf, perc_bonf = summarize_region(lead_file, repl_file)

            results.append({
                "region": region_path.split("/name")[0],
                "N_tested": n_tested,
                "N_replicated": n_rep,
                "Replication_%": perc, 
                "N_replicated_bonf": n_rep_bonf,
                "Replication__bonf%": perc_bonf, 
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
    df.to_csv(args.out, sep="\t", index=False)
    print(f"Saved summary: {args.out}")

    # -------------------------
    # Publication-quality plot
    # -------------------------

    plt.rcParams.update({
        "font.size": 8,
        "axes.linewidth": 0.8,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "font.family": "sans-serif",
    })
 
    # Split hemispheres
    df_left  = df[df["region"].str.endswith("_left")].copy()
    df_right = df[df["region"].str.endswith("_right")].copy()
 
    # Strip suffix for y-axis labels
    df_left["label"]  = df_left["region"].str.replace("_left$",  "", regex=True)
    df_right["label"] = df_right["region"].str.replace("_right$", "", regex=True)
 
    # Sort each hemisphere independently by nominal replication rate
    df_left  = df_left.sort_values("Replication_%", ascending=True).reset_index(drop=True)
    df_right = df_right.sort_values("Replication_%", ascending=True).reset_index(drop=True)
 
    COLOR_NOM  = "#0072B2"   # blue  — nominal threshold (p < 0.05)
    COLOR_BONF = "#D55E00"   # orange — Bonferroni-corrected threshold (p < 0.05/n_loci)
    DOT_SIZE   = 25
    BAR_HEIGHT = 0.25        # vertical offset between the two dots per region
 
    def draw_panel(ax, df_hemi, title):
        """Draw one hemisphere panel."""
        n = len(df_hemi)
        y_positions = range(n)
 
        for i, row in df_hemi.iterrows():
            y = i
            # Horizontal line connecting the two dots
            ax.hlines(
                y + BAR_HEIGHT / 2,
                min(row["Replication_%"], row["Replication__bonf%"]),
                max(row["Replication_%"], row["Replication__bonf%"]),
                color="grey", linewidth=0.6, alpha=0.5, zorder=1
            )
            # Nominal dot (p < 0.05)
            ax.scatter(
                row["Replication_%"],
                y + BAR_HEIGHT,
                s=DOT_SIZE, color=COLOR_NOM,
                zorder=2, label="p < 0.05" if i == 0 else ""
            )
            # Bonferroni dot (p < 0.05 / n_loci)
            ax.scatter(
                row["Replication__bonf%"],
                y,
                s=DOT_SIZE, color=COLOR_BONF, marker="D",
                zorder=2, label="p < 0.05 / n loci" if i == 0 else ""
            )
 
        ax.set_yticks([i + BAR_HEIGHT / 2 for i in range(n)])
        ax.set_yticklabels(df_hemi["label"], fontsize=7)
        ax.set_xlim(0, 100)
        ax.set_ylim(-0.5, n)
        ax.set_xlabel("Replication rate (%)", fontsize=8)
        ax.set_title(title, fontsize=9, fontweight="bold", pad=6)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(axis="x", linestyle="--", linewidth=0.4, alpha=0.35)
 
    # Figure dimensions: two panels side by side
    n_left  = len(df_left)
    n_right = len(df_right)
    fig_height = max(n_left, n_right) * 0.32 + 1.2   # ~0.32 in per region + margins
 
    fig, (ax_left, ax_right) = plt.subplots(
        1, 2,
        figsize=(7.2, fig_height),   # Nature full-width ≈ 7.2 in
        sharey=False
    )
 
    draw_panel(ax_left,  df_left,  "Left hemisphere")
    draw_panel(ax_right, df_right, "Right hemisphere")
 
    # Shared legend — placed above the figure
    handles = [
        plt.scatter([], [], s=DOT_SIZE, color=COLOR_NOM,  label="p < 0.05"),
        plt.scatter([], [], s=DOT_SIZE, color=COLOR_BONF, marker="D",
                    label=r"p < 0.05 / $n_{\mathrm{loci}}$"),
    ]
    fig.legend(
        handles=handles,
        loc="upper center",
        ncol=2,
        frameon=False,
        fontsize=8,
        bbox_to_anchor=(0.5, 1.01)
    )
 
    plt.tight_layout(rect=[0, 0, 1, 0.97])
 
    png_out = args.out.replace(".tsv", "_replication.png")
    plt.savefig(png_out, format="png", dpi=300, bbox_inches="tight")
    print(f"Saved figure: {png_out}")


if __name__ == "__main__":
    main()