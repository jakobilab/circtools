#!/usr/bin/env python3

# Replacement for circtools_reconstruct_visualization.R
#
# Reads FUCHS / circtools reconstruct output files from the current working
# directory (*.exon_counts.bed, *.mate_status.txt, *alternative_splicing.txt,
# *.coverage_profiles/cluster_association.all_circles.tsv) and produces a
# multi-page results.pdf with the same plots as the original R script.
#
# Usage:
#   python circtools_reconstruct_visualization.py \
#       <reconstruct_data_path> \
#       <group_indices>          \   # comma-separated 1-based ints, e.g. "1,1,2,2"
#       <condition_list>         \   # comma-separated labels,      e.g. "Control,Treated"
#       <colour_mode>                # "colour" or "bw"
#
# Significance brackets are drawn with statannotations (if installed) or
# omitted gracefully when the package is absent.

import argparse
import glob
import os
import subprocess
import sys
import warnings
from itertools import combinations

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages
from scipy import stats

try:
    from statannotations.Annotator import Annotator
    HAS_STATANNOTATIONS = True
except ImportError:
    HAS_STATANNOTATIONS = False


# ---------------------------------------------------------------------------
# data loading (mirrors the awk / sed / grep shell calls in the R script)
# ---------------------------------------------------------------------------

def load_exon_counts():
    """
    awk '{split($11,array,","); sum=0; for(i=1;i<=length(array);i++){sum+=array[i]};
          print FILENAME"\t"$3-$2"\t"$10"\t"sum}' *.exon_counts.bed |
    sed 's/\.exon_counts\.bed//g' | egrep -v '\W0\W'
    Columns: Library, Length_total, Exons, Length_exons
    """
    rows = []
    for path in glob.glob('*.exon_counts.bed'):
        lib = path.replace('.exon_counts.bed', '')
        with open(path) as fh:
            for line in fh:
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 11:
                    continue
                try:
                    length_total = int(parts[2]) - int(parts[1])
                    exons = parts[9]
                    length_exons = sum(int(x) for x in parts[10].split(',') if x)
                except (ValueError, IndexError):
                    continue
                # egrep -v '\W0\W' — skip lines whose non-word-bounded token is "0"
                joined = f'\t{length_total}\t{exons}\t{length_exons}'
                if '\t0\t' in joined or joined.endswith('\t0'):
                    continue
                rows.append({'Library': lib, 'Length_total': length_total,
                             'Exons': exons, 'Length_exons': length_exons})
    return pd.DataFrame(rows)


def load_mate_status():
    """
    awk '{print FILENAME"\t"$6"\t"$7"\t"$8}' *.mate_status.txt |
    sed 's/\.mate_status\.txt//g' | grep -v single |
    awk '{print $0"\t"$3/($2+$3+$4)}'
    Columns: Library, Single, Double, Unknown, Ratio
    """
    rows = []
    for path in glob.glob('*.mate_status.txt'):
        lib = path.replace('.mate_status.txt', '')
        with open(path) as fh:
            for line in fh:
                if 'single' in line.lower():
                    continue
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 8:
                    continue
                try:
                    s, d, u = int(parts[5]), int(parts[6]), int(parts[7])
                    total = s + d + u
                    ratio = d / total if total else 0.0
                except (ValueError, IndexError):
                    continue
                rows.append({'Library': lib, 'Single': s, 'Double': d,
                             'Unknown': u, 'Ratio': ratio})
    return pd.DataFrame(rows)


def load_isoforms():
    """
    grep -H -o -n ',' *alternative_splicing.txt | cut -d: -f1,2 | uniq -c |
    sed 's/\.alternative.*//g' | awk '{print $2"\t"$1}'
    Columns: Library, num
    """
    rows = []
    for path in glob.glob('*alternative_splicing.txt'):
        lib = path.split('.alternative')[0]
        counts = {}
        with open(path) as fh:
            for i, line in enumerate(fh, 1):
                n = line.count(',')
                if n:
                    key = (lib, i)
                    counts[key] = n
        for (l, _), n in counts.items():
            rows.append({'Library': l, 'num': n})
    return pd.DataFrame(rows)


def load_circ_length():
    """
    awk '{if($2>0){print FILENAME"\t"$2}}' *.coverage_profiles/cluster_association.all_circles.tsv |
    sed 's/\.coverage_profiles\/cluster_association\.all_circles\.tsv//g' | grep -v length
    Columns: Library, length
    """
    rows = []
    for path in glob.glob('*.coverage_profiles/cluster_association.all_circles.tsv'):
        lib = path.replace('.coverage_profiles/cluster_association.all_circles.tsv', '')
        try:
            df = pd.read_csv(path, sep='\t')
        except Exception:
            continue
        for _, row in df.iterrows():
            try:
                if int(row['length']) > 0:
                    rows.append({'Library': lib, 'length': int(row['length'])})
            except (ValueError, KeyError):
                pass
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# plotting helpers
# ---------------------------------------------------------------------------

PALETTE_COLOUR = plt.rcParams['axes.prop_cycle'].by_key()['color']
PALETTE_BW     = [str(v) for v in np.linspace(0.1, 0.9, 10)]


def get_palette(mode, n):
    base = PALETTE_BW if mode == 'bw' else PALETTE_COLOUR
    return [base[i % len(base)] for i in range(n)]


def assign_groups(df, lib_col, lib_names, condition_names):
    """Add a 'group' column by mapping library names to condition labels."""
    mapping = {lib: condition_names[i] for i, lib in enumerate(lib_names)}
    df = df.copy()
    df['group'] = df[lib_col].map(mapping)
    return df


def boxplot_page(pdf, df, y_col, group_col, title, subtitle, ylabel,
                 colour_mode, comparisons, log_y=False):
    """Draw one boxplot page and save to the PdfPages object."""
    groups = df[group_col].dropna().unique()
    palette = get_palette(colour_mode, len(groups))
    colour_map = {g: palette[i] for i, g in enumerate(sorted(groups))}

    fig, ax = plt.subplots(figsize=(8, 5.5))
    group_data = [df.loc[df[group_col] == g, y_col].dropna().values
                  for g in sorted(groups)]
    bp = ax.boxplot(group_data, labels=sorted(groups), patch_artist=True)
    for patch, grp in zip(bp['boxes'], sorted(groups)):
        patch.set_facecolor(colour_map[grp])

    if log_y:
        ax.set_yscale('log')

    # significance annotations (statannotations)
    if HAS_STATANNOTATIONS and len(comparisons) > 0:
        try:
            annotator = Annotator(ax, comparisons, data=df, x=group_col, y=y_col)
            annotator.configure(test='Mann-Whitney', text_format='star', loc='inside')
            annotator.apply_and_annotate()
        except Exception:
            pass  # graceful degradation

    ax.set_title(f'{title}\n{subtitle}', fontsize=11)
    ax.set_ylabel(ylabel)
    ax.set_xlabel('')
    ax.tick_params(axis='x', bottom=False, labelbottom=False)

    # legend patches
    patches = [mpatches.Patch(color=colour_map[g], label=g) for g in sorted(groups)]
    ax.legend(handles=patches, title='Sample', loc='upper right')

    fig.tight_layout()
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)


def quantile_page(pdf, circ_length, group_col, colour_mode, max_val):
    """One density + quantile fill plot per group."""
    groups = sorted(circ_length[group_col].dropna().unique())
    palette = get_palette(colour_mode, 4)   # 4 quantile bands
    quantile_probs = [0.90, 0.95, 0.99]

    for grp in groups:
        sub = circ_length.loc[circ_length[group_col] == grp, 'length'].dropna().values
        if len(sub) < 2:
            continue
        from scipy.stats import gaussian_kde
        kde = gaussian_kde(sub)
        xs = np.linspace(0, max_val, 500)
        ys = kde(xs)
        quantiles = np.quantile(sub, quantile_probs)

        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(xs, ys, color='black', lw=0.5)

        band_labels = ['< 90%', '90–95%', '95–99%', '> 99%']
        edges = [0] + list(quantiles) + [max_val + 1]
        for band_i in range(4):
            mask = (xs >= edges[band_i]) & (xs < edges[band_i + 1])
            ax.fill_between(xs, 0, ys, where=mask,
                            color=palette[band_i], label=band_labels[band_i], alpha=0.8)

        ax.set_title(f'Circular RNA reconstruction results\n'
                     f'CircRNA length – quantile plot for sample "{grp}"')
        ax.set_xlabel('CircRNA length')
        ax.set_ylabel('Density')
        ax.set_xlim(0, max_val)
        ax.legend(title='Quantiles', loc='upper right')
        fig.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Produce reconstruction summary PDF (replacement for circtools_reconstruct_visualization.R).'
    )
    parser.add_argument('reconstruct_data', help='Path to reconstruct data directory (sets cwd)')
    parser.add_argument('group_indices',
                        help='Comma-separated 1-based group indices, e.g. "1,1,2,2"')
    parser.add_argument('condition_list',
                        help='Comma-separated condition names, e.g. "Control,Treated"')
    parser.add_argument('colour_mode', choices=['colour', 'bw'],
                        help='Colour scheme: "colour" or "bw"')
    args = parser.parse_args()

    os.chdir(args.reconstruct_data)

    group_indices   = [int(x) for x in args.group_indices.split(',')]
    condition_list  = args.condition_list.split(',')
    colour_mode     = args.colour_mode

    # Build library → condition name mapping
    # The R script uses: names[group_indices[x]] for each library in order
    # We infer library order from isoforms (as the R script does)
    exon_counts = load_exon_counts()
    mate_status  = load_mate_status()
    isoforms     = load_isoforms()
    circ_length  = load_circ_length()

    for df, name in [(exon_counts, 'exon_counts'), (mate_status, 'mate_status'),
                     (isoforms, 'isoforms'), (circ_length, 'circ_length')]:
        if df.empty:
            print(f'WARNING: no data loaded for {name}')

    # Sequencing library names (alphabetical to match R levels())
    all_libs = sorted(isoforms['Library'].unique()) if not isoforms.empty else []
    if len(all_libs) != len(group_indices):
        print(f'WARNING: {len(all_libs)} libraries found but {len(group_indices)} '
              f'group indices provided. Mapping may be incorrect.')

    detailed_names = [condition_list[group_indices[i] - 1]
                      for i in range(len(all_libs))]

    def add_group(df, lib_col='Library'):
        return assign_groups(df, lib_col, all_libs, detailed_names)

    exon_counts = add_group(exon_counts) if not exon_counts.empty else exon_counts
    mate_status  = add_group(mate_status) if not mate_status.empty else mate_status
    isoforms     = add_group(isoforms)   if not isoforms.empty else isoforms
    circ_length  = add_group(circ_length) if not circ_length.empty else circ_length

    # circRNA count per library
    if not circ_length.empty:
        count_df = circ_length.groupby('Library').size().reset_index(name='count')
        count_df = add_group(count_df)
    else:
        count_df = pd.DataFrame(columns=['Library', 'count', 'group'])

    comparisons = list(combinations(sorted(condition_list), 2))

    max_len = int(circ_length['length'].max()) if not circ_length.empty else 10000

    # ---- write PDF ----------------------------------------------------------
    out_pdf = 'results.pdf'
    with PdfPages(out_pdf) as pdf:

        # 1. Absolute circRNA count (normal scale)
        if not count_df.empty:
            boxplot_page(pdf, count_df, 'count', 'group',
                         'Circular RNA reconstruction results',
                         'Absolute count – normal scale',
                         'Total number of circular RNAs',
                         colour_mode, comparisons, log_y=False)

        # 2. Absolute circRNA count (log10 scale)
        if not count_df.empty:
            boxplot_page(pdf, count_df, 'count', 'group',
                         'Circular RNA reconstruction results',
                         'Absolute count – log10 scale',
                         'Total number of circular RNAs (log10)',
                         colour_mode, comparisons, log_y=True)

        # 3. Quantile length plots (one page per group)
        if not circ_length.empty:
            quantile_page(pdf, circ_length, 'group', colour_mode, max_len)

        # 4. Isoforms per host gene
        if not isoforms.empty:
            boxplot_page(pdf, isoforms, 'num', 'group',
                         'Circular RNA reconstruction results',
                         '# circRNA isoforms per host gene',
                         'Number of isoforms',
                         colour_mode, comparisons, log_y=False)

        # 5. Total length (exons + introns), log scale
        if not exon_counts.empty:
            boxplot_page(pdf, exon_counts, 'Length_total', 'group',
                         'Circular RNA reconstruction results',
                         'Total length of circRNAs',
                         'Total length (exons + introns)',
                         colour_mode, comparisons, log_y=True)

        # 6. Exon-based length, log scale
        if not exon_counts.empty:
            boxplot_page(pdf, exon_counts, 'Length_exons', 'group',
                         'Circular RNA reconstruction results',
                         'Exon-based length of circRNAs',
                         'Length (exons only)',
                         colour_mode, comparisons, log_y=True)

        # 7. Number of exons per circRNA
        if not exon_counts.empty:
            boxplot_page(pdf, exon_counts, 'Exons', 'group',
                         'Circular RNA reconstruction results',
                         'Number of exons per circRNA',
                         '# exons per circRNA',
                         colour_mode, comparisons, log_y=False)

        # 8. Single breakpoints (absolute), log scale
        if not mate_status.empty:
            boxplot_page(pdf, mate_status, 'Single', 'group',
                         'Circular RNA reconstruction results',
                         'Single breakpoint circRNAs (absolute)',
                         '# single breakpoints',
                         colour_mode, comparisons, log_y=True)

        # 9. Double breakpoints (absolute), log scale
        if not mate_status.empty:
            boxplot_page(pdf, mate_status, 'Double', 'group',
                         'Circular RNA reconstruction results',
                         'Double breakpoint circRNAs (absolute)',
                         '# double breakpoints',
                         colour_mode, comparisons, log_y=True)

        # 10. Ratio of double breakpoints
        if not mate_status.empty:
            boxplot_page(pdf, mate_status, 'Ratio', 'group',
                         'Circular RNA reconstruction results',
                         'Ratio of double breakpoints',
                         'Ratio double breakpoints',
                         colour_mode, comparisons, log_y=False)

    print(f'Written: {out_pdf}')


if __name__ == '__main__':
    main()
