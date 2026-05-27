#!/usr/bin/env python3
"""
QQ plot for GWAS summary statistics.
 
Computes:
  - Observed vs expected -log10(p) QQ plot with 95% confidence band
  - Genomic inflation factor lambda GC
    on observed ~ expected chi2 scores
"""
 
import argparse
import os
 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats
from postgwas_tools.annot.utils import read_sumstats
 
# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------
 
def compute_lambda_gc(pvalues: np.ndarray) -> float:
    """
    Genomic inflation factor λ_GC.
 
    λ_GC = median(χ²_observed) / median(χ²_expected under H0)
         = median(χ²_observed) / 0.4549   [χ²(1) median]
 
    A value of 1.0 means no inflation.
    Values > 1 indicate systematic inflation (population stratification,
    cryptic relatedness, or true polygenic signal).
 
    Parameters
    ----------
    pvalues : array of p-values (already cleaned, no NaN/0/1 extremes required
              but clipping is applied internally).
 
    Returns
    -------
    lambda_gc : float
    """
    pvalues = np.clip(pvalues, 1e-300, 1.0)           # avoid -inf in chi2
    chi2_obs = stats.chi2.isf(pvalues, df=1)           # inverse survival = ppf(1-p)
    chi2_null_median = stats.chi2.ppf(0.5, df=1)       # ≈ 0.4549
    return float(np.median(chi2_obs) / chi2_null_median)

 
 
# ---------------------------------------------------------------------------
# Confidence band
# ---------------------------------------------------------------------------
 
def _qqplot_ci(n: int, alpha: float = 0.95) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Pointwise 95 % confidence band for a QQ plot of uniform order statistics,
    using the Beta distribution (exact for order statistics of U[0,1]).
 
    For rank k out of n:
        U_(k) ~ Beta(k, n-k+1)
 
    Returns
    -------
    expected   : -log10 of expected p-value midpoints
    ci_lower   : -log10 lower CI bound
    ci_upper   : -log10 upper CI bound
    """
    ranks = np.arange(1, n + 1)
    p_low  = (1 - alpha) / 2
    p_high = 1 - p_low
 
    # Beta quantiles for each order statistic
    ci_lo = stats.beta.ppf(p_low,  ranks, n - ranks + 1)
    ci_hi = stats.beta.ppf(p_high, ranks, n - ranks + 1)
    expected = ranks / (n + 1)                     # median of Beta(k, n-k+1)
 
    # Clip to avoid log10(0)
    ci_lo    = np.clip(ci_lo,    1e-300, 1)
    ci_hi    = np.clip(ci_hi,    1e-300, 1)
    expected = np.clip(expected, 1e-300, 1)
 
    return (
        -np.log10(expected[::-1]),
        -np.log10(ci_hi[::-1]),
        -np.log10(ci_lo[::-1]),
    )
 
 
# ---------------------------------------------------------------------------
# Main plotting function
# ---------------------------------------------------------------------------
 
def qqplot(
    pvalues: np.ndarray,
    title: str = "QQ plot",
    ax=None,
    point_color: str = "#2166ac",
    ci_color: str = "#b2d8f7",
    diagonal_color: str = "#d73027",
    marker: str = "o",
    markersize: float = 3.0,
    alpha_points: float = 0.7,
    lam_gc = None,
) -> plt.Axes:
    """
    Draw a GWAS QQ plot from an array of p-values.
 
    Parameters
    ----------
    pvalues        : 1-D array-like of raw p-values.
    title          : Plot title.
    ax             : Existing matplotlib Axes (created if None).
    point_color    : Colour for observed data points.
    ci_color       : Fill colour for the 95 % confidence band.
    diagonal_color : Colour for the y = x reference line.
    marker         : Marker style.
    markersize     : Marker size.
    alpha_points   : Transparency of data points.
    show_lambda    : Annotate with λ_GC.
 
    Returns
    -------
    ax : matplotlib Axes
    """
    pvalues = np.asarray(pvalues, dtype=float)
    pvalues = pvalues[np.isfinite(pvalues)]
    pvalues = np.clip(pvalues, 1e-300, 1.0)
    n = len(pvalues)
 
    if n == 0:
        raise ValueError("No valid p-values to plot.")
 
    # ---- Expected & observed -log10(p) ------------------------------------
    expected_unif = (np.arange(1, n + 1) - 0.5) / n          # uniform spacing
    expected_log  = -np.log10(np.sort(expected_unif))         # ascending x-axis
    observed_log  = -np.log10(np.sort(pvalues))               # ascending y-axis
 
    # ---- Confidence band ---------------------------------------------------
    exp_ci, ci_lower, ci_upper = _qqplot_ci(n)
 
    # ---- Statistics --------------------------------------------------------
    lam_gc = compute_lambda_gc(pvalues)
 
    # ---- Plot --------------------------------------------------------------
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 6), facecolor="w", edgecolor="k")
 
    # Confidence band (fill between)
    ax.fill_between(
        exp_ci, ci_lower, ci_upper,
        color=ci_color, alpha=0.5, label="95% CI",
        linewidth=0,
    )
 
    # Diagonal y = x
    max_val = max(expected_log.max(), observed_log.max()) * 1.05
    ax.plot([0, max_val], [0, max_val],
            color=diagonal_color, linewidth=1.2,
            linestyle="--", label="y = x", zorder=3)
 
    # Observed points
    ax.scatter(
        expected_log, observed_log,
        color=point_color, s=markersize**2,
        marker=marker, alpha=alpha_points,
        linewidths=0, zorder=4, label="Observed",
    )
 
    # ---- Annotations -------------------------------------------------------
    annotation_lines = []
    if lam_gc is not None:
        annotation_lines.append(rf"$\lambda_{{GC}}$ = {lam_gc:.4f}")
 
    if annotation_lines:
        annotation_text = "\n".join(annotation_lines)
        ax.text(
            0.05, 0.95, annotation_text,
            transform=ax.transAxes,
            fontsize=9,
            verticalalignment="top",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                      edgecolor="#cccccc", alpha=0.9),
        )
 
    # ---- Formatting --------------------------------------------------------
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.set_xlabel(r"Expected $-\log_{10}(p)$", fontsize=12)
    ax.set_ylabel(r"Observed $-\log_{10}(p)$", fontsize=12)
    ax.set_title(title, fontsize=13, fontweight="bold")
    ax.legend(loc="lower right", fontsize=9, framealpha=0.9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
 
    return ax
 
 
# ---------------------------------------------------------------------------
# CLI entry-point
# ---------------------------------------------------------------------------
 
def main():
    parser = argparse.ArgumentParser(
        description="Generate a QQ plot for GWAS summary statistics."
    )
    parser.add_argument("-p", "--path", type=str, required=True,
                        help="File path to the summary statistics file.")
    parser.add_argument("-o", "--out", type=str, default=None,
                        help="Output folder for saving the plot.")
    parser.add_argument("--title", type=str, default="QQ plot",
                        help="Plot title (default: 'QQ plot').")
    parser.add_argument("--dpi", type=int, default=150,
                        help="Resolution of saved figure (default: 150).")
    args = parser.parse_args()
 
    file_path     = args.path
    output_folder = args.out
    print(f"Working with file: {file_path}")
 
    df = read_sumstats(file_path)
    df = df.dropna(subset=["P"])
 
    pvalues = df["P"].values
 
    # ---- Compute & print stats --------------------------------------------
    lam_gc              = compute_lambda_gc(pvalues)
 
    print(f"\n{'='*45}")
    print(f"  N SNPs (after QC)   : {len(pvalues):,}")
    print(f"  λ_GC                : {lam_gc:.4f}")
    print(f"{'='*45}\n")
 
    # ---- Plot --------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(6, 6), facecolor="w", edgecolor="k")
    qqplot(pvalues=pvalues, title=args.title, ax=ax, lam_gc=lam_gc)
 
    if not output_folder:
        path_to_save = os.path.join(os.path.dirname(file_path), "QQplot.png")
    else:
        os.makedirs(output_folder, exist_ok=True)
        path_to_save = os.path.join(output_folder, "QQplot.png")
 
    fig.tight_layout()
    fig.savefig(path_to_save, dpi=args.dpi, bbox_inches="tight")
    print(f"Figure saved at: {path_to_save}")
 
 
if __name__ == "__main__":
    main()