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
from scipy.ndimage import gaussian_filter1d



def smoothing(x):
    x = np.asarray(x, dtype=float)
    n = len(x)
    sigma = max(1, round(n * 0.15))
    return gaussian_filter1d(x, sigma=sigma)

def plot_coverage(coverage_file, output_folder):
    basename = os.path.basename(coverage_file)          
    parts = basename.split('.')
    circle_id = parts[0]
    transcript_name = parts[1] if len(parts) > 2 else ''

    df = pd.read_csv(coverage_file, sep='\t')

    smoothed = smoothing(df['coverage'].values)

   
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
