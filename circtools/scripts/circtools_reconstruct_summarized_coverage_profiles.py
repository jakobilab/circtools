#!/usr/bin/env python3

# Replacement for circtools_reconstruct_summarized_coverage_profiles.R
#
# For every coverage-profile .txt file in a folder:
#   - normalises positions to a [0, 1] grid (101 bins)
#   - builds a summary table
#   - plots average coverage for all / short (<500 bp) / medium (500-999 bp) / long (>=1000 bp) circles
#   - k-means clusters each size class and writes cluster assignments + means
#
# Usage:
#   python circtools_reconstruct_summarized_coverage_profiles.py <folder>

import argparse
import glob
import os
import sys
import warnings

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.cluster.vq import kmeans, vq, whiten


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def choose_centers(n):
    """Mirror the R choose_centers() logic."""
    centers = 1
    if n > 2:   centers = 2
    if n >= 10: centers = 3
    if n >= 20: centers = 4
    if n >= 100:
        centers = round(n / 20)
    if centers >= 10:
        centers = 10
    return centers


def load_profile(filepath):
    """
    Read one coverage-profile file and return (circle_id, length, coverage_101)
    where coverage_101 is a numpy array of 101 values on a normalised grid,
    or None if the file cannot be reduced to exactly 101 bins.
    """
    try:
        df = pd.read_csv(filepath, sep='\t')
    except Exception:
        return None

    if 'relative_pos_in_circle' not in df.columns or 'coverage' not in df.columns:
        return None

    pos = df['relative_pos_in_circle'].values.astype(float)
    if pos[-1] == 0:
        return None
    pos = np.round(pos / pos[-1], 2)

    # tapply equivalent: mean coverage per rounded position
    df_tmp = pd.DataFrame({'pos': pos, 'coverage': df['coverage'].values})
    agg = df_tmp.groupby('pos')['coverage'].mean()

    if len(agg) != 101:
        return None

    basename = os.path.basename(filepath)
    circle_id = basename.split('.')[0]
    length = len(df)

    return circle_id, length, agg.values


def build_summary_table(folder):
    files = sorted(glob.glob(os.path.join(folder, '*.txt')))
    rows = []
    for f in files:
        result = load_profile(f)
        if result is None:
            continue
        circle_id, length, cov101 = result
        rows.append({'circle_id': circle_id, 'length': length, **{f'b{i}': v for i, v in enumerate(cov101)}})

    if not rows:
        return pd.DataFrame()

    df = pd.DataFrame(rows).drop_duplicates()
    df = df.dropna(subset=[c for c in df.columns if c.startswith('b')])
    return df


def cov_cols(df):
    return [c for c in df.columns if c.startswith('b')]


# ---------------------------------------------------------------------------
# plotting helpers
# ---------------------------------------------------------------------------

def plot_mean_profile(ax, df, title):
    cols = cov_cols(df)
    means = df[cols].mean()
    ax.plot(means.values, lw=2, color='dodgerblue')
    ax.set_title(title, fontsize=14)
    ax.set_xlabel('relative position', fontsize=12)
    ax.set_ylabel('avg. coverage', fontsize=12)


def plot_cluster(ax, center, size):
    ax.plot(center, lw=5, color='dodgerblue')
    ax.set_title(f'Cluster Size: {size}', fontsize=14)
    ax.set_xlabel('relative position', fontsize=12)
    ax.set_ylabel('avg. coverage', fontsize=12)


def run_and_save_clustering(df, label, folder):
    """
    Cluster df (if enough rows), save PDFs + TSVs.
    Returns (cluster_df, centers_df) or (None, None).
    """
    n = len(df)
    centers = choose_centers(n)

    if centers <= 1:
        print(f'ERROR not enough {label} circles in table to perform a clustering ({n} circles)')
        return None, None

    cols = cov_cols(df)
    X = df[cols].values.astype(float)
    X_w = whiten(X)

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        cluster_centers, _ = kmeans(X_w, centers, iter=50)
        labels, _ = vq(X_w, cluster_centers)

    pdf_path = os.path.join(folder, f'coverage.clusters.{label}.pdf')
    with PdfPages(pdf_path) as pdf:
        for i in range(centers):
            size = int((labels == i).sum())
            if size > 0:
                fig, ax = plt.subplots(figsize=(8, 4))
                plot_cluster(ax, cluster_centers[i], size)
                pdf.savefig(fig, bbox_inches='tight')
                plt.close(fig)

    assoc = pd.DataFrame({
        'circle_id': df['circle_id'].values,
        'length':    df['length'].values,
        'cluster_id': labels + 1   # 1-based to match R output
    })
    assoc.to_csv(os.path.join(folder, f'cluster_association.{label}.tsv'), sep='\t', index=False)

    centers_df = pd.DataFrame(cluster_centers, columns=cols)
    centers_df.index = centers_df.index + 1
    centers_df.to_csv(os.path.join(folder, f'cluster_means.{label}.tsv'), sep='\t')

    return assoc, centers_df



def main():
    parser = argparse.ArgumentParser(
        description='Summarise and cluster circRNA coverage profiles from a folder of .txt files.'
    )
    parser.add_argument('folder', help='Folder containing coverage-profile .txt files')
    args = parser.parse_args()

    folder = args.folder
    if not os.path.isdir(folder):
        sys.exit(f'ERROR: {folder} is not a directory')

    print('Loading coverage profiles …')
    summary = build_summary_table(folder)

    if summary.empty:
        sys.exit('ERROR: no valid coverage profiles found in the folder')

    short  = summary[summary['length'] <  500]
    medium = summary[(summary['length'] >= 500) & (summary['length'] < 1000)]
    long_  = summary[summary['length'] >= 1000]

    # ---- average-profile PDF ------------------------------------------------
    avg_pdf_path = os.path.join(folder, 'coverage_profiles.all_circles.pdf')
    with PdfPages(avg_pdf_path) as pdf:
        for subset, title in [
            (summary, f'All circles ({len(summary)})'),
            (short,   f'Short circles ({len(short)})'),
            (medium,  f'Medium circles ({len(medium)})'),
            (long_,   f'Long circles ({len(long_)})'),
        ]:
            if len(subset) == 0:
                print(f'ERROR no circles in subset: {title}')
                continue
            fig, ax = plt.subplots(figsize=(8, 4))
            plot_mean_profile(ax, subset, title)
            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)

    print(f'Written: {avg_pdf_path}')

    # ---- clustering ---------------------------------------------------------
    for subset, label in [
        (summary, 'all_circles'),
        (short,   'short_circles'),
        (medium,  'medium_circles'),
        (long_,   'long_circles'),
    ]:
        run_and_save_clustering(subset, label, folder)
        if len(subset) > 1:
            print(f'Clustering done for {label}')

    print('Done.')


def run(folder):
    import sys
    _saved = sys.argv
    sys.argv = ['circtools_reconstruct_summarized_coverage_profiles', folder]
    main()
    sys.argv = _saved


if __name__ == '__main__':
    main()
