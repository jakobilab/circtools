#!/usr/bin/env python3

# Replacement for circtools_reconstruct_summary_graphs.R
#
# Reads mate_status.genes.txt files for four brain-region / age combinations
# and produces boxplots of single / double / rolling-circle fractions and
# circle lengths — matching the layout of the original R script.
#
# Usage (hardcoded filenames mirror the original R script):
#   python circtools_reconstruct_summary_graphs.py [output.pdf]
#
# Optionally supply an output PDF path; default is "summary_graphs.pdf".

import argparse
import sys
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


SAMPLE_FILES = {
    'Old Cerebellum':   'old_cerebellum.mate_status.genes.txt',
    'Young Cerebellum': 'young_cerebellum.mate_status.genes.txt',
    'Old Hippocampus':  'old_hippocampus.mate_status.genes.txt',
    'Young Hippocampus':'young_hippocampus.mate_status.genes.txt',
}

SHORT_NAMES = ['oCB', 'yCB', 'oHC', 'yHC']   # same order as SAMPLE_FILES


def load(path):
    """Load a mate_status.genes.txt, return DataFrame or None."""
    if not os.path.isfile(path):
        print(f'WARNING: {path} not found — skipping')
        return None
    return pd.read_csv(path, sep='\t')


def fraction_boxplot(ax, datasets, labels, col, total_col, title, hlines=(0.25, 0.5, 0.75)):
    """Boxplot of col/total_col ratios, one box per dataset."""
    data = []
    for df in datasets:
        if df is not None and col in df.columns and total_col in df.columns:
            data.append((df[col] / df[total_col]).dropna().values)
        else:
            data.append(np.array([]))

    ax.boxplot(data, labels=labels, patch_artist=True,
               boxprops=dict(facecolor='lightblue'))
    for h in hlines:
        ax.axhline(h, linestyle='--', color='grey', lw=0.8)
    ax.set_title(title)
    ax.set_ylabel('fraction')


def length_boxplot(ax, datasets, labels, col, title, log_y=True):
    """Boxplot of a length column across datasets."""
    data = []
    for df in datasets:
        if df is not None and col in df.columns:
            data.append(df[col].dropna().values)
        else:
            data.append(np.array([]))

    ax.boxplot(data, labels=labels, patch_artist=True,
               boxprops=dict(facecolor='lightblue'))
    if log_y:
        ax.set_yscale('log')
        ax.axhline(500, linestyle='--', color='grey', lw=0.8)
    ax.set_title(title)
    ax.set_ylabel(col)


def main():
    parser = argparse.ArgumentParser(
        description='Summary boxplots for FUCHS mate-status files (four brain samples).'
    )
    parser.add_argument('output_pdf', nargs='?', default='summary_graphs.pdf',
                        help='Output PDF path (default: summary_graphs.pdf)')
    args = parser.parse_args()

    frames = {name: load(path) for name, path in SAMPLE_FILES.items()}
    dfs   = list(frames.values())                 # same order as SAMPLE_FILES
    names = list(frames.keys())
    short = SHORT_NAMES

    with PdfPages(args.output_pdf) as pdf:

        # ---- Page 1: single / double / rolling fractions (2×2) -------------
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        axes = axes.flatten()
        for ax, (name, df) in zip(axes, frames.items()):
            dlist = [df] if df is not None else [None]
            # single panel with three boxes: single, double, undefined
            cols  = ['single', 'double', 'undefined']
            bdata = []
            for c in cols:
                if df is not None and c in df.columns and 'num_reads' in df.columns:
                    bdata.append((df[c] / df['num_reads']).dropna().values)
                else:
                    bdata.append(np.array([]))
            bp = ax.boxplot(bdata, labels=['single', 'double', 'rolling'], patch_artist=True,
                            boxprops=dict(facecolor='lightblue'))
            for h in (0.25, 0.5, 0.75):
                ax.axhline(h, linestyle='--', color='grey', lw=0.8)
            ax.set_title(name)
            ax.set_ylabel('fraction of reads')

        fig.suptitle('Read-fraction by breakpoint class', fontsize=14)
        fig.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

        # ---- Page 2: min_length boxplot -------------------------------------
        fig, ax = plt.subplots(figsize=(8, 5))
        length_boxplot(ax, dfs, short, 'min_length', 'Circle length (min)', log_y=True)
        fig.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

        # ---- Page 3: max_length boxplot -------------------------------------
        fig, ax = plt.subplots(figsize=(8, 5))
        length_boxplot(ax, dfs, short, 'max_length', 'Circle length (max)', log_y=True)
        fig.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

    print(f'Written: {args.output_pdf}')


if __name__ == '__main__':
    main()
