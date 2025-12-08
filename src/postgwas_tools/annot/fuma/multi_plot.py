# coding: utf-8
##########################################################################
# NSAp - Copyright (C) CEA, 2023
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
# Antoine Dufournet
##########################################################################

from postgwas_tools.annot.utils import (
    find_files,
    adjust_color_brightness,
    read_sumstats,
)

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os


def plot_manhattan(file_paths, output_folder, y_max=None):
    # List of distinct base colors for each model
    base_colors = ['blue', 'red', 'green', 'purple', 'orange', 'brown', 'pink', 'gray', 'olive', 'cyan']

    plt.figure(figsize=(21, 10.5))

    # Initialize chrom_offsets for alignment
    chrom_offsets = None
    dic_start_end_chr = {}

    # Loop through each file and plot data
    for i, file in enumerate(file_paths):
        print(f"Working with file: {file}")
        model = file.split('/')[-2]
        
        df = read_sumstats(file)
        
        if i == 0:
            # Calculate cumulative position for x-axis alignment only once
            df.sort_values(by=['CHR', 'BP'], inplace=True)
            chromosomes = sorted(df['CHR'].unique())
            chrom_offsets = df.groupby('CHR')['BP'].max().cumsum().shift(fill_value=0)
            padding_values = [30_000_000]*11 + [3_000_000]*3 + [30_000_000]*5 + [3_000_000]*3
            padding_series = pd.Series(padding_values, index=chrom_offsets.index)
            padding_cumsum = padding_series.cumsum().shift(fill_value=0)
            # Final offsets with padding between chromosomes
            chrom_offsets += padding_cumsum

            for chrom in chromosomes:
                dic_start_end_chr[chrom] = [df[df["CHR"]==chrom]["BP"].min(),df[df["CHR"]==chrom]["BP"].max()]

        # Filter SNPs
        df = df[df["P"] < 1e-2] #1e-3
        df['neg_log10_pval'] = -np.log10(df['P'])
        df['x_val'] = df.apply(lambda row: row['BP'] + chrom_offsets.loc[row['CHR']], axis=1)

        # Loop through chromosomes to alternate colors
        chromosomes = sorted(df['CHR'].unique())
        for chrom_idx, chrom in enumerate(chromosomes):
            chrom_data = df[df['CHR'] == chrom]

            # Adjust color brightness based on chromosome index
            base_color = base_colors[i % len(base_colors)]  # Base color based on the model
            brightness_factor = 1.5 if chrom_idx % 2 == 0 else 0.5  # Stronger light or dark shades
            adjusted_color = adjust_color_brightness(base_color, brightness_factor)

            # Plot data for the chromosome
            plt.scatter(chrom_data['x_val'], chrom_data['neg_log10_pval'], 
                        color=adjusted_color, s=7,
                        label=model if chrom_idx == 0 else "")  # Add label only once per model

    # Add a horizontal significance threshold line
    plt.axhline(y=-np.log10(0.05/1000000), color='r', linestyle='--')
    plt.axhline(y=-np.log10(0.05/10000), color='g', linestyle='--')

    # Customize plot labels and title
    plt.xlabel('Chromosome', fontsize=22)
    plt.ylabel(r"$-log_{10}{(p)}$", fontsize=22)
    plt.title('Manhattan Plot', fontsize=24)

    # Add chromosome labels at their midpoints
    chromosome_ticks = [chrom_offsets[chrom] + df[df['CHR'] == chrom]['BP'].max() / 2 for chrom in sorted(df['CHR'].unique())]
    chromosome_labels = [f"{chrom}" for chrom in sorted(df['CHR'].unique())]
    plt.xticks(chromosome_ticks, chromosome_labels, fontsize=16) #rotation=45
    plt.xlim([dic_start_end_chr[1][0] + chrom_offsets[1] - 13_000_000, dic_start_end_chr[22][1] + chrom_offsets[22]+ 13_000_000])
    plt.ylim(bottom=0)
    if y_max:
        plt.ylim(top=y_max) 
    #plt.legend(loc='upper right',
    #           fontsize=16,
    #           frameon=False     )
    plt.yticks(fontsize=16)
    plt.tight_layout()
    output_path = f'manhattan_plot.png'
    plt.savefig(f"{output_folder}/{output_path}", format='png', transparent=True)
    print(f"Plot saved to: {output_folder}/{output_path}")
    #plt.show()

def plot_miami(file_paths, output_folder, y_max=None):
    base_colors = ['blue', 'red', 'green', 'purple', 'orange', 'brown', 'pink', 'gray', 'olive', 'cyan']
    chrom_offsets = None
    dic_start_end_chr = {}
    min_p_top = 1e-2
    min_p_bot = 1e-2
    if not y_max:
        y_max = 300
    
    fig, axes = plt.subplots(2, 1, figsize=(21, 10.5), sharex=True, gridspec_kw={'hspace': 0.12})
    pos = axes[0].get_position()
    axes[0].set_position([pos.x0, pos.y0 - 0.02, pos.width, pos.height])

    for i, file in enumerate(file_paths):  # Expecting exactly 2 files
        print(f"Reading {file}")
        df = read_sumstats(file)
        df["P"] = pd.to_numeric(df["P"], errors="coerce")
        df = df.dropna(subset=["P"])

        if i == 0:
            df.sort_values(by=['CHR', 'BP'], inplace=True)
            chromosomes = sorted(df['CHR'].unique())
            chrom_offsets = df.groupby('CHR')['BP'].max().cumsum().shift(fill_value=0)
            padding_values = [30_000_000]*11 + [3_000_000]*3 + [30_000_000]*5 + [3_000_000]*3
            padding_series = pd.Series(padding_values, index=chrom_offsets.index)
            padding_cumsum = padding_series.cumsum().shift(fill_value=0)
            chrom_offsets += padding_cumsum
            for chrom in chromosomes:
                dic_start_end_chr[chrom] = [
                    df[df["CHR"] == chrom]["BP"].min(),
                    df[df["CHR"] == chrom]["BP"].max()
                ]

        df = df[df["P"] < 1e-2]
        df['neg_log10_pval'] = -np.log10(df['P'])
        df['x_val'] = df.apply(lambda row: row['BP'] + chrom_offsets.loc[row['CHR']], axis=1)

        j = i % 2
        chromosomes = sorted(df['CHR'].unique())

        for chrom_idx, chrom in enumerate(chromosomes):
            chrom_data = df[df['CHR'] == chrom]
            base_color = base_colors[j % len(base_colors)]
            brightness = 1.5 if chrom_idx % 2 == 0 else 0.5
            adjusted_color = adjust_color_brightness(base_color, brightness)
            y_vals = chrom_data['neg_log10_pval'] if j == 0 else -chrom_data['neg_log10_pval']
            axes[j].scatter(chrom_data['x_val'], y_vals, color=adjusted_color, s=7)

        sig1 = -np.log10(0.05 / 1_000_000)
        sig2 = -np.log10(0.05 / 10_000)

        if j == 0:
            min_p_top = min(df["P"].min(), min_p_top)
            top_y_limit = min(y_max, -np.log10(min_p_top))
            axes[j].axhline(y=sig1, color='r', linestyle='--')
            axes[j].axhline(y=sig2, color='g', linestyle='--')
            axes[j].set_ylim(0, top_y_limit)
        else:
            min_p_bot = min(df["P"].min(), min_p_bot)
            bottom_y_limit = min(y_max, -np.log10(min_p_bot))
            bottom_y_limit = max(top_y_limit, bottom_y_limit)
            axes[j].axhline(y=-sig1, color='r', linestyle='--')
            axes[j].axhline(y=-sig2, color='g', linestyle='--')
            axes[j].set_ylim(-bottom_y_limit, 0)

        # Remove 0 label from y-axis
        yticks = axes[j].get_yticks()
        yticks = [t for t in yticks if t != 0]
        axes[j].set_yticks(yticks)
        axes[j].set_yticklabels([str(int(abs(t))) for t in yticks], fontsize=16)

    # Frame line cleanup
    axes[0].spines['bottom'].set_visible(True)
    axes[1].spines['top'].set_visible(True)

    # Chromosome labels between plots
    chromosome_ticks = [
        chrom_offsets[chrom] + df[df['CHR'] == chrom]['BP'].max() / 2
        for chrom in sorted(df['CHR'].unique())
    ]
    chromosome_labels = [str(chrom) for chrom in sorted(df['CHR'].unique())]
    axes[0].xaxis.set_ticks_position('none')
    axes[1].xaxis.set_ticks_position('top')
    axes[1].set_xticks(chromosome_ticks)
    axes[1].set_xticklabels(chromosome_labels, fontsize=15)
    axes[1].tick_params(axis='x', which='both', direction='out', pad=0, length=0)
    axes[1].set_xlabel('Chromosome', fontsize=22, labelpad=10)
    axes[1].xaxis.set_label_position('bottom')


    axes[0].set_ylabel(r"$-log_{10}(p)$", fontsize=22)
    axes[0].yaxis.set_label_coords(-0.03, 0)
    plt.suptitle("Miami Plot", fontsize=26, y=0.93)
    plt.xlim([
        dic_start_end_chr[1][0] + chrom_offsets[1] - 13_000_000,
        dic_start_end_chr[22][1] + chrom_offsets[22] + 13_000_000
    ])

    output_path = os.path.join(output_folder, "miami_plot.png")
    plt.savefig(output_path, format='png', transparent=True, bbox_inches="tight", pad_inches=0.1)
    print(f"Saved Miami subplot to {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Generate Manhattan plots for one or multiple GWAS summary statistics."
    )
    parser.add_argument(
        '-p', '--paths', nargs='+', type=str, required=True,
        help="List of paths, file names, or wildcard patterns (e.g. /data/*.sumstats)"
    )
    parser.add_argument(
        '-k', '--kind', type=str, default='manhattan',
        choices=['manhattan', 'miami'],
        help="Type of plot to generate: 'manhattan' (default) or 'miami'"
    )
    parser.add_argument(
        '-o', '--out', type=str, required=True,
        help="Output folder for saving the plot."
    )
    parser.add_argument(
        '--ymax', type=int, default=None,
        help="Value of the max -log10(p-value) to be plot"
    )
    args = parser.parse_args()

    # Find files matching any path or pattern
    file_paths = find_files(args.paths)
    plot_type = args.kind
    output_folder = args.out
    output_folder = output_folder.replace('/', '') if output_folder.endswith('/') else output_folder

    if not file_paths:
        print("No files found. Please check your paths or patterns.")
        return

    print(f"Found {len(file_paths)} file(s):")
    for f in file_paths:
        print(f"   - {f}")

    if plot_type == 'manhattan':
        plot_manhattan(file_paths, output_folder, y_max=args.ymax)
    elif plot_type == 'miami':
        plot_miami(file_paths, output_folder, y_max=args.ymax)

if __name__ == "__main__":
    main()
