#!/usr/bin/env python3

# Replacement for circtools_reconstruct_coverage_graph.R
# Loads a circRNA coverage profile and plots a smoothed coverage graph per transcript.

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def smoothing(x):
    """
    Smooth a coverage array by taking local means in a window of size ~1% of
    the total length, centred on each position (mirrors the R implementation).
    """
    x = np.asarray(x, dtype=float)
    n = len(x)
    w = max(1, round(n * 0.01))
    smoothed = np.empty(n)
    for i in range(n):
        lo = max(0, i - w)
        hi = min(n, i + w + 1)
        smoothed[i] = np.nanmean(x[lo:hi])
    return smoothed


def plot_coverage(coverage_file, output_folder):
    # Derive identifiers from the filename
    basename = os.path.basename(coverage_file)           # e.g. chr1:100|200.NM_001.txt
    parts = basename.split('.')
    circle_id = parts[0]
    transcript_name = parts[1] if len(parts) > 2 else ''

    df = pd.read_csv(coverage_file, sep='\t')

    smoothed = smoothing(df['coverage'].values)

    # Colour each bar by exon number (exon column + 1 maps to matplotlib colour cycle)
    exon_vals = df['exon'].values
    unique_exons = np.unique(exon_vals)
    cmap = plt.get_cmap('tab10')
    colours = [cmap((e) % 10) for e in exon_vals]

    title_str = f"{circle_id}\n{transcript_name}"
    xlabel_str = f"Exon: {exon_vals.min()} - Exon: {exon_vals.max()}"
    out_pdf = os.path.join(output_folder, f"{circle_id}_{transcript_name}.pdf")

    with PdfPages(out_pdf) as pdf:
        fig, ax = plt.subplots(figsize=(10, 5))
        x = np.arange(len(smoothed))
        ax.bar(x, smoothed, color=colours, width=1.0, linewidth=0)
        ax.set_title(title_str)
        ax.set_xlabel(xlabel_str)
        ax.set_ylabel('number of reads')
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

    print(f"Written: {out_pdf}")


def main():
    parser = argparse.ArgumentParser(
        description='Plot smoothed coverage graph for a circRNA coverage profile.'
    )
    parser.add_argument('coverage_file', help='Coverage profile .txt file (tab-separated with header)')
    parser.add_argument('output_folder', help='Folder where the PDF will be written')
    args = parser.parse_args()

    os.makedirs(args.output_folder, exist_ok=True)
    plot_coverage(args.coverage_file, args.output_folder)


if __name__ == '__main__':
    main()
