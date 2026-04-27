#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
circtools_quickcheck_wrapper.py
Python equivalent of the R-based circtools_quickcheck_wrapper
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

sns.set(style="whitegrid")

def parse_args():
    parser = argparse.ArgumentParser(description="circtools quickcheck (Python equivalent)")

    parser.add_argument("dcc_dir", help="Path to circtools detect output directory")
    parser.add_argument("star_dir", help="Path to STAR mapping directory")
    parser.add_argument("output_prefix", help="Output file prefix (PDF will be created)")
    parser.add_argument("conditions", help="Comma-separated list of condition names")
    parser.add_argument("grouping", help="Comma-separated list of numeric group indices")
    parser.add_argument("colour_mode", choices=["colour", "bw"], help="Plot colour mode")
    parser.add_argument("cleanup_string", help="Regex cleanup pattern for sample names")
    parser.add_argument("starfolder_suffix", help="STAR folder suffix to use or 0")
    parser.add_argument("remove_suffix_chars", type=int)
    parser.add_argument("remove_prefix_chars", type=int)
    parser.add_argument("remove_columns", help="Comma-separated column indices to remove (or 0)")

    return parser.parse_args()


def read_star_unique_mappings(star_folder):
    log_file = os.path.join(star_folder, "Log.final.out")
    if not os.path.exists(log_file):
        return np.nan
    with open(log_file) as f:
        lines = f.readlines()
        for line in lines:
            if "Uniquely mapped reads %" in line:
                # extract percentage
                return float(line.split()[-1].replace("%", ""))
    return np.nan


def main():
    args = parse_args()

    # prepare arguments
    cond_list = args.conditions.split(",")
    grouping = list(map(int, args.grouping.split(",")))

    # parse remove_columns safely — ignore '0'
    remove_cols = []
    if args.remove_columns and args.remove_columns.strip() not in ["0", ""]:
        try:
            remove_cols = [int(x) for x in args.remove_columns.split(",") if x.strip().isdigit()]
        except ValueError:
            remove_cols = []


    # read DCC data
    circ = pd.read_csv(os.path.join(args.dcc_dir, "CircRNACount"), sep="\t")
    linear = pd.read_csv(os.path.join(args.dcc_dir, "LinearCount"), sep="\t")
    
    # Determine where numeric data starts in each file
    circ_base_cols = circ.columns[:5].tolist()
    linear_base_cols = linear.columns[:5].tolist()

    circ_data_start = 4 if any(col.lower() == "strand" for col in circ_base_cols) else 3
    linear_data_start = 4 if any(col.lower() == "strand" for col in linear_base_cols) else 3

    print(f"DEBUG circ_data_start = {circ_data_start}, linear_data_start = {linear_data_start}")

        
    

    if remove_cols:
        circ.drop(circ.columns[remove_cols], axis=1, inplace=True)
        linear.drop(linear.columns[remove_cols], axis=1, inplace=True)

    # clean column names
    circ.columns = [c.replace(args.cleanup_string, "") for c in circ.columns]
    linear.columns = [c.replace(args.cleanup_string, "") for c in linear.columns]

    circ.columns = list(circ.columns[:circ_data_start]) + [
        c[args.remove_prefix_chars:len(c) - args.remove_suffix_chars]
        for c in circ.columns[circ_data_start:]
    ]
    linear.columns = list(linear.columns[:linear_data_start]) + [
        c[args.remove_prefix_chars:len(c) - args.remove_suffix_chars]
        for c in linear.columns[linear_data_start:]
    ]

    circ.iloc[:, circ_data_start:] = circ.iloc[:, circ_data_start:].apply(pd.to_numeric, errors="coerce").fillna(0)
    linear.iloc[:, linear_data_start:] = linear.iloc[:, linear_data_start:].apply(pd.to_numeric, errors="coerce").fillna(0)


    # sum counts per library
    circ_sums = circ.iloc[:, circ_data_start:].sum()
    lin_sums = linear.iloc[:, linear_data_start:].sum()
    samples = circ.columns[circ_data_start:]
    num_samples = len(samples)

        
    # get STAR mappings
    star_subdirs = [os.path.join(args.star_dir, d) for d in os.listdir(args.star_dir)
                    if os.path.isdir(os.path.join(args.star_dir, d))]
    if args.starfolder_suffix != "0":
        star_subdirs = [d for d in star_subdirs if d.endswith(args.starfolder_suffix)]

    unique_mappings = []
    for subdir in star_subdirs[:num_samples]:
        unique_mappings.append(read_star_unique_mappings(subdir))
    if len(unique_mappings) < num_samples:
        unique_mappings.extend([np.nan] * (num_samples - len(unique_mappings)))

    # derive groups and colours
    # derive groups and colours
    if len(grouping) < num_samples:
        repeats = (num_samples // len(grouping)) + 1
        grouping = (grouping * repeats)[:num_samples]

    # safety: clamp out-of-range group indices
    groups = []
    for i, g in enumerate(grouping[:num_samples]):
        if isinstance(g, (int, float)) and 1 <= g <= len(cond_list):
            groups.append(cond_list[int(g) - 1])
        else:
            print(f"⚠️  Invalid group index {g} for sample {samples[i] if i < len(samples) else 'unknown'}, "
                  f"falling back to last condition '{cond_list[-1]}'")
            groups.append(cond_list[-1])

    # If fewer groups than samples, repeat until matching
    if len(groups) < num_samples:
        repeats = (num_samples // len(groups)) + 1
        groups = (groups * repeats)[:num_samples]

    print(f"\nDEBUG condition list: {cond_list} (len={len(cond_list)})")
    print(f"DEBUG grouping indices: {grouping} (len={len(grouping)})")
    print(f"DEBUG resolved groups: {groups} (len={len(groups)})")
    print(f"DEBUG sample names: {list(samples)} (len={len(samples)})\n")

    # --- compute detected circles ---
    num_circles = (circ.iloc[:, circ_data_start:] > 0).sum()
    print(f"DEBUG num_circles length = {len(num_circles)}, head:\n{num_circles.head()}\n")

    # --- prepare dataframes for plotting ---
    df_raw = pd.DataFrame({
        "sample": samples,
        "circ": pd.to_numeric(circ_sums.values, errors="coerce"),
        "linear": pd.to_numeric(lin_sums.values, errors="coerce"),
        "group": [str(g) for g in groups]
    })

    df_circles = pd.DataFrame({
        "sample": samples,
        "unique": pd.to_numeric(unique_mappings, errors="coerce"),
        "circles": pd.to_numeric(num_circles.values, errors="coerce"),
        "group": [str(g) for g in groups]
    })

    # avoid division by zero or NaN in ratio
    unique_safe = np.array(unique_mappings, dtype=float)
    unique_safe[unique_safe == 0] = np.nan
    ratio_values = pd.to_numeric(num_circles.values, errors="coerce") / (unique_safe / 1e6)

    df_ratio = pd.DataFrame({
        "sample": samples,
        "num": ratio_values,
        "group": [str(g) for g in groups]
    }).sort_values("num")
    
        # --- PCA on circRNA counts ---
    circ_matrix = circ.iloc[:, circ_data_start:].T  # samples x circRNAs

    # safety: numeric + fill
    circ_matrix = circ_matrix.apply(pd.to_numeric, errors="coerce").fillna(0)

    do_pca = circ_matrix.shape[0] >= 2 and circ_matrix.sum(axis=1).gt(0).all()

    if do_pca:
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(circ_matrix)

        pca = PCA(n_components=2)
        pcs = pca.fit_transform(X_scaled)

        pca_df = pd.DataFrame(
            pcs,
            columns=["PC1", "PC2"],
            index=samples
        )
        pca_df["group"] = groups

        explained = pca.explained_variance_ratio_
    else:
        print("⚠️  Skipping PCA — insufficient or zero-count samples")


    # --- debug prints for verification ---
    print("DEBUG df_raw info:")
    print(df_raw.info())
    print(df_raw.head(), "\n")

    print("DEBUG df_circles info:")
    print(df_circles.info())
    print(df_circles.head(), "\n")

    print("DEBUG df_ratio info:")
    print(df_ratio.info())
    print(df_ratio.head(), "\n")


    # --- plotting section ---
    pdf_name = args.output_prefix + ".pdf"
    with PdfPages(pdf_name) as pdf:
        def clean_label(name: str) -> str:
            """Remove STAR/DCC suffixes for prettier sample labels."""
            return name.replace(".Chimeric.out.junction", "").replace("_Chimeric.out.junction", "")

        # --- Plot 1: Circular vs Linear Counts ---
        plt.figure(figsize=(11.69, 8.27))
        sns.scatterplot(data=df_raw, x="circ", y="linear", hue="group", style="group", s=80, edgecolor="black")
        plt.xscale("log")
        plt.yscale("log")
        plt.xlabel("Circular RNA read count (log scale)")
        plt.ylabel("Linear RNA read count (log scale)")
        plt.title("Circular vs. Linear Read Counts per Library", fontweight="bold")
        plt.grid(True, linestyle="--", linewidth=0.5, alpha=0.6)
        plt.legend(title="Group", loc="upper right", frameon=True)
        for i, row in df_raw.iterrows():
            plt.text(row["circ"], row["linear"], clean_label(row["sample"]),
                    fontsize=9, ha="left", va="bottom")
        plt.tight_layout()
        pdf.savefig()
        plt.close()

        # --- Plot 2: Unique Reads vs. Detected circRNAs (>1 count) ---
        print("DEBUG generating Plot 2 (Unique Reads vs Detected circRNAs)")

        df_circles_clean = pd.DataFrame({
            "sample": [clean_label(s) for s in samples],
            "unique": pd.to_numeric(unique_mappings, errors="coerce"),
            "circles": pd.to_numeric(num_circles.values, errors="coerce"),
            "group": [str(g) for g in groups]
        }).dropna(subset=["unique", "circles"])

        if df_circles_clean.empty:
            print("⚠️  Skipping Plot 2 — no valid 'unique' or 'circles' values\n")
        else:
            plt.figure(figsize=(11.69, 8.27))
            sns.scatterplot(data=df_circles_clean, x="unique", y="circles", hue="group", style="group", s=80, edgecolor="black")
            for i, row in df_circles_clean.iterrows():
                plt.text(row["unique"], row["circles"], row["sample"], fontsize=9, ha="left", va="bottom")
            plt.xlabel("Unique Mapped Reads")
            plt.ylabel("Detected circRNAs (>1 count)")
            plt.title("Unique Reads vs. Detected circRNAs", fontweight="bold")
            plt.grid(True, linestyle="--", linewidth=0.5, alpha=0.6)
            plt.legend(title="Group", loc="upper right", frameon=True)
            plt.tight_layout()
            pdf.savefig()
            plt.close()

        # --- Plot 3: Detected circRNAs per Million Unique Reads ---
        print("DEBUG generating Plot 3 (Detected circRNAs per Million Unique Reads)")

        unique_array = np.array(unique_mappings, dtype=float)
        unique_array[unique_array == 0] = np.nan
        norm_counts = num_circles.values / (unique_array / 1e6)

        df_ratio_clean = pd.DataFrame({
            "sample": [clean_label(s) for s in samples],
            "norm": norm_counts,
            "group": [str(g) for g in groups]
        }).dropna(subset=["norm"])

        if df_ratio_clean.empty:
            print("⚠️  Skipping Plot 3 — no valid normalized counts\n")
        else:
            plt.figure(figsize=(11.69, 8.27))
            order = df_ratio_clean.sort_values("norm")
            bars = plt.bar(order["sample"], order["norm"], color="#3874c8", edgecolor="black")
            plt.xticks(rotation=45, ha="right")
            plt.xlabel("Sample")
            plt.ylabel("Detected circRNAs per Million Unique Reads")
            plt.title("Detected circRNAs per Million Unique Reads", fontweight="bold")
            plt.grid(axis="y", linestyle="--", linewidth=0.5, alpha=0.6)
            plt.legend(["Group"], loc="upper right", frameon=True)
            for idx, val in enumerate(order["norm"]):
                plt.text(idx, val, f"{val:.0f}", ha="center", va="bottom", fontsize=8)
            plt.tight_layout()
            pdf.savefig()
            plt.close()
            
                # --- Plot 4: PCA of circRNA counts ---
        if do_pca:
            print("DEBUG generating Plot 4 (PCA of circRNA counts)")

            plt.figure(figsize=(11.69, 8.27))
            sns.scatterplot(
                data=pca_df,
                x="PC1",
                y="PC2",
                hue="group",
                style="group",
                s=100,
                edgecolor="black"
            )

            for sample, row in pca_df.iterrows():
                plt.text(
                    row["PC1"],
                    row["PC2"],
                    clean_label(sample),
                    fontsize=9,
                    ha="left",
                    va="bottom"
                )

            plt.xlabel(f"PC1 ({explained[0]*100:.1f}%)")
            plt.ylabel(f"PC2 ({explained[1]*100:.1f}%)")
            plt.title("PCA of circRNA Read Counts", fontweight="bold")
            plt.grid(True, linestyle="--", linewidth=0.5, alpha=0.6)
            plt.legend(title="Group", loc="best", frameon=True)
            plt.tight_layout()
            pdf.savefig()
            plt.close()


    print(f"✅ QuickCheck report saved: {pdf_name}\n")





if __name__ == "__main__":
    main()
