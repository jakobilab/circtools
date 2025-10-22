#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Python equivalent of circtools_quickcheck_wrapper.R
Recreates the three-page QC report using matplotlib + seaborn.
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from datetime import datetime

sns.set(style="whitegrid", context="talk")


def parse_args():
    parser = argparse.ArgumentParser(description="circtools quickcheck (Python equivalent)")
    parser.add_argument("dcc_dir", help="Path to DCC output directory")
    parser.add_argument("star_dir", help="Path to STAR mapping directory")
    parser.add_argument("output_prefix", help="Output PDF prefix")
    parser.add_argument("conditions", help="Comma-separated list of condition names")
    parser.add_argument("grouping", help="Comma-separated list of numeric group indices")
    parser.add_argument("colour_mode", choices=["colour", "bw"], help="Colour or black/white mode")
    parser.add_argument("cleanup_string", help="Regex cleanup pattern for sample names")
    parser.add_argument("starfolder_suffix", help="STAR folder suffix or 0")
    parser.add_argument("remove_suffix_chars", type=int)
    parser.add_argument("remove_prefix_chars", type=int)
    parser.add_argument("remove_columns", help="Comma-separated list of columns to remove (or 0)")
    return parser.parse_args()


def read_unique_mappings(star_folder):
    """Read STAR Log.final.out and extract uniquely mapped reads."""
    log_file = os.path.join(star_folder, "Log.final.out")
    if not os.path.exists(log_file):
        return np.nan

    with open(log_file, "r") as f:
        lines = f.readlines()

    # Extract line 7 (R script assumes row 7, column 2)
    try:
        unique_line = lines[6]
        val = unique_line.split("|")[1].strip().replace("\t", "")
        return float(val)
    except Exception:
        return np.nan


def main():
    args = parse_args()
    cond_list = args.conditions.split(",")
    grouping = [int(x) for x in args.grouping.split(",")]

    # Read DCC data
    circ = pd.read_csv(os.path.join(args.dcc_dir, "CircRNACount"), sep="\t", header=0)
    linear = pd.read_csv(os.path.join(args.dcc_dir, "LinearCount"), sep="\t", header=0)

    # Drop columns if specified
    if args.remove_columns != "0":
        remove_cols = [int(x) for x in args.remove_columns.split(",")]
        circ.drop(circ.columns[remove_cols], axis=1, inplace=True)
        linear.drop(linear.columns[remove_cols], axis=1, inplace=True)

    # Clean sample names
    circ.columns = [c.replace(args.cleanup_string, "") for c in circ.columns]
    linear.columns = [c.replace(args.cleanup_string, "") for c in linear.columns]

    # Trim/remove characters from sample names
    def clean_names(cols):
        new_cols = []
        for c in cols:
            new_c = c[args.remove_prefix_chars:]
            if args.remove_suffix_chars > 0:
                new_c = new_c[: -args.remove_suffix_chars]
            new_cols.append(new_c)
        return new_cols

    circ.columns = list(circ.columns[:3]) + clean_names(list(circ.columns[3:]))
    linear.columns = list(linear.columns[:3]) + clean_names(list(linear.columns[3:]))

    samples = circ.columns[3:]
    num_samples = len(samples)
    print(f"Found {num_samples} data columns")

    # Sum counts per library
    circ_sum = circ.iloc[:, 3:].sum()
    lin_sum = linear.iloc[:, 3:].sum()

    # STAR runs
    star_dirs = [os.path.join(args.star_dir, d) for d in os.listdir(args.star_dir)
                 if os.path.isdir(os.path.join(args.star_dir, d)) and "mate" not in d]
    if args.starfolder_suffix != "0":
        star_dirs = [d for d in star_dirs if d.endswith(args.starfolder_suffix)]

    # Read unique mappings
    unique_reads = [read_unique_mappings(d) for d in star_dirs[:num_samples]]
    if len(unique_reads) < num_samples:
        unique_reads += [np.nan] * (num_samples - len(unique_reads))

    # Compute number of circles (x > 1)
    num_circles = (circ.iloc[:, 3:] > 1).sum()

    # Grouping
    if len(grouping) < num_samples:
        grouping = (grouping * (num_samples // len(grouping)))[:num_samples]
    groups = [cond_list[g - 1] for g in grouping]

    # Build data frames for plotting
    df_raw = pd.DataFrame({
        "circ": circ_sum.values,
        "linear": lin_sum.values,
        "group": groups
    }, index=samples)

    df_circles = pd.DataFrame({
        "unique": unique_reads,
        "circles": num_circles.values,
        "group": groups
    }, index=samples)

    df_ratio = pd.DataFrame({
        "name": samples,
        "num": num_circles.values / (np.array(unique_reads) / 1e6),
        "group": groups
    }).sort_values("num")

    # Plot setup
    palette = "colorblind" if args.colour_mode == "colour" else "gray"
    pdf_path = f"{args.output_prefix}.pdf"
    date_str = datetime.now().strftime("%Y-%m-%d %H:%M")

    with PdfPages(pdf_path) as pdf:

        # Page 1 — Circular vs Linear
        plt.figure(figsize=(11.69, 8.27))
        sns.scatterplot(data=df_raw, x="circ", y="linear", hue="group", style="group",
                        palette=palette, s=80, edgecolor="black")
        plt.xscale("log"); plt.yscale("log")
        plt.xlabel("Circular RNA read count (log scale)")
        plt.ylabel("Linear RNA read count (log scale)")
        plt.title("Circular vs. linear read counts throughout selected libraries")
        plt.suptitle(f"based on data from {len(star_dirs)} mapped libraries\ncreated: {date_str}")
        plt.tight_layout()
        pdf.savefig(); plt.close()

        # Page 2 — Unique mapped reads vs circles
        plt.figure(figsize=(11.69, 8.27))
        sns.scatterplot(data=df_circles, x="unique", y="circles", hue="group", style="group",
                        palette=palette, s=80, edgecolor="black")
        plt.xscale("log"); plt.yscale("log")
        plt.xlabel("Number of mapped reads (log scale)")
        plt.ylabel("Number of detected circles (log scale)")
        plt.title("Number of mapped reads vs. number of detected circles per library")
        plt.suptitle(f"based on data from {len(star_dirs)} mapped libraries\ncreated: {date_str}")
        plt.tight_layout()
        pdf.savefig(); plt.close()

        # Page 3 — Circles per million reads
        plt.figure(figsize=(11.69, 8.27))
        sns.barplot(data=df_ratio, x="name", y="num", hue="group", palette=palette, edgecolor="black")
        plt.xticks(rotation=45, ha="right")
        plt.xlabel("Sequencing library")
        plt.ylabel("Number of detected circular RNAs")
        plt.title("Detected circular RNAs per million unique mapped reads")
        plt.suptitle(f"based on data from {len(star_dirs)} mapped libraries\ncreated: {date_str}")
        plt.tight_layout()
        pdf.savefig(); plt.close()

    print(f"✅ QuickCheck report created: {pdf_path}")


if __name__ == "__main__":
    main()
