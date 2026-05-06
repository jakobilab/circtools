#!/usr/bin/env python3

import pandas as pd
import os
from scripts.circtest_functions import circ_filter, circ_test
from openpyxl import Workbook
from openpyxl.utils.dataframe import dataframe_to_rows
import logging
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf as pdf_backend
import numpy as np
import warnings
import matplotlib.patches as mpatches


logger = logging.getLogger(__name__)
MAX_LINES = 50000




COLOUR_PALETTE = ["#F4A582", "#2DBDB6"]   
BW_PALETTE     = ["#AAAAAA", "#333333"]  


def circ_ratio_plot(circ_row, linear_row, coord_row, group_indicators,
                    lab_legend, y_axis_range, colour_mode, ax):

    circ_vals   = np.array(circ_row.values,   dtype=float)
    linear_vals = np.array(linear_row.values, dtype=float)
    total       = circ_vals + linear_vals
    ratio       = np.where(total > 0, circ_vals / total, np.nan)

    groups      = np.array(group_indicators)
    uniq_grp    = list(dict.fromkeys(groups))   # preserve order
    palette     = COLOUR_PALETTE if colour_mode == "colour" else BW_PALETTE
    bar_width   = 0.45
    x_positions = np.arange(len(uniq_grp))

    for g_idx, grp in enumerate(uniq_grp):
        mask      = groups == grp
        grp_ratio = ratio[mask]
        mean_val  = np.nanmean(grp_ratio)
        n_valid   = np.sum(~np.isnan(grp_ratio))
        sem_val   = np.nanstd(grp_ratio, ddof=1) / np.sqrt(n_valid) if n_valid > 1 else 0.0
        color     = palette[g_idx % len(palette)]

        ax.bar(
            x_positions[g_idx], mean_val,
            width=bar_width,
            color=color,
            edgecolor="black",
            linewidth=0.8,
            zorder=2
        )

        if sem_val > 0:
            ax.errorbar(
                x_positions[g_idx], mean_val,
                yerr=sem_val,
                fmt="none",
                ecolor="black",
                elinewidth=1.2,
                capsize=4,
                capthick=1.2,
                zorder=3
            )

    gene  = coord_row.get("Gene", coord_row.get("gene", "unknown"))
    chrom = coord_row.iloc[0]
    start = coord_row.iloc[1]
    end   = coord_row.iloc[2]
    ax.set_title(
        f"Annotation: {gene}\nChr {chrom}, {start}, {end}",
        fontsize=11, loc="left", pad=8
    )

    y_max = y_axis_range if y_axis_range > 0 else 1.0
    ax.set_ylim(0, y_max)
    ax.set_xlim(-0.5, len(uniq_grp) - 0.5)
    ax.set_xticks([])
    ax.set_ylabel("circRNA/(circRNA + Linear RNA)", fontsize=10)
    ax.set_xlabel("")

    y_ticks = np.arange(0, y_max + 0.001, 0.25)
    ax.set_yticks(y_ticks)
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:.2f}"))

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    legend_patches = [
        mpatches.Patch(
            facecolor=palette[i % len(palette)],
            edgecolor="black",
            linewidth=0.6,
            label=grp
        )
        for i, grp in enumerate(uniq_grp)
    ]
    ax.legend(
        handles=legend_patches,
        title=lab_legend if lab_legend else None,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.08),
        ncol=len(uniq_grp),
        frameon=False,
        fontsize=10,
        handlelength=1.2,
        handleheight=1.2
    )


def run_circ_test_wrapper(args):

    print("Loading CircRNACount")
    CircRNACount = pd.read_csv(
        os.path.join(args.dcc_data, "CircRNACount"), sep="\t", header=0
    )

    print("Loading LinearCount")
    LinearCount = pd.read_csv(
        os.path.join(args.dcc_data, "LinearCount"), sep="\t", header=0
    )

    print("Loading CircCoordinates")
    CircCoordinates = pd.read_csv(
        os.path.join(args.dcc_data, "CircCoordinates"), sep="\t", header=0
    )

    cond_cols_1based = [int(x) for x in args.condition_columns.split(",")]
    cond_cols = [x - 1 for x in cond_cols_1based]        

    circle_description = list(range(3))                   

    cond_names = args.condition_list.split(",")            
    group_indices = [int(x) for x in args.groups.split(",")] 


    final_grouping = [cond_names[g - 1] for g in group_indices]

    print(f"[INFO] Condition columns (0-based): {cond_cols}")
    print(f"[INFO] Group labels (per sample):   {final_grouping}")


    all_cols = circle_description + cond_cols

    circ_slice  = CircRNACount.iloc[:, all_cols]
    linear_slice = LinearCount.iloc[:, all_cols]


    print("Filtering circRNA counts")
    CircRNACount_filtered = circ_filter(
        circ=circ_slice,
        linear=linear_slice,
        Nreplicates=args.replicates,
        filter_sample=args.filter_sample,
        filter_count=args.filter_count,
        percentage=args.percent_filter,
        circle_description=circle_description
    )


    print("Filtering circRNA coordinates")
    CircCoordinates_filtered = CircCoordinates.loc[CircRNACount_filtered.index]

    print("Filtering linear RNA counts")
    LinearCount_filtered = linear_slice.loc[CircRNACount_filtered.index]

    print("Running circTest")
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning, 
                                message="invalid value encountered in sqrt")
        results = circ_test(
            Circ=CircRNACount_filtered,
            Linear=LinearCount_filtered,
            CircCoordinates=CircCoordinates_filtered,
            group=final_grouping,
            alpha=args.max_fdr,
            circle_description=circle_description
        )

    summary = results["summary_table"]

    if len(summary) == 0:
        print("No candidates to plot, exiting.")
        return


    print("Generating plots")
    max_n = min(args.max_plots, len(summary))

    if args.only_negative:
        plot_rows = [
            i for i in summary.index
            if "direction" in summary.columns and summary.loc[i, "direction"] < 0
        ]
    else:
        plot_rows = list(summary.index[:max_n])

    if plot_rows:
        pdf_path = args.output_name + ".pdf"
        with pdf_backend.PdfPages(pdf_path) as pdf:
            plots_per_page = 2
            for page_start in range(0, len(plot_rows), plots_per_page):
                batch = plot_rows[page_start : page_start + plots_per_page]
                fig, axes = plt.subplots(plots_per_page, 1, figsize=(8.2, 11.69))
                if plots_per_page == 1:
                    axes = [axes]

                for ax, circ_id in zip(axes, batch):
                    n_desc = len(circle_description)
                    circ_row   = CircRNACount_filtered.loc[circ_id].iloc[n_desc:]
                    linear_row = LinearCount_filtered.loc[circ_id].iloc[n_desc:]
                    coord_row  = CircCoordinates_filtered.loc[circ_id]

                    circ_ratio_plot(
                        circ_row, linear_row, coord_row,
                        group_indicators=final_grouping,
                        lab_legend=args.output_label,
                        y_axis_range=args.y_range,
                        colour_mode=args.colour_mode,
                        ax=ax
                    )

                for ax in axes[len(batch):]:
                    ax.set_visible(False)

                pdf.savefig(fig, bbox_inches="tight")
                plt.close(fig)

        print(f"Plots saved to {pdf_path}")

    print("Saving CSV")
    csv_data = summary.head(MAX_LINES).dropna()     
    csv_data.to_csv(
        args.output_name + ".csv",
        sep="\t",
        index=False,
        header=args.add_header
    )

    print("Saving Excel")
    wb = Workbook()

    ws1 = wb.active
    ws1.title = "Significant circles"
    sig_df = summary.head(MAX_LINES)
    ws1.append(list(sig_df.columns))
    for row in sig_df.itertuples(index=False):
        ws1.append(list(row))

    sig_idx = sig_df.index
    ws2 = wb.create_sheet("Circle Counts")
    for row in dataframe_to_rows(
        circ_slice.loc[sig_idx], index=False, header=True
    ):
        ws2.append(row)

    ws3 = wb.create_sheet("Linear Counts")
    for row in dataframe_to_rows(
        linear_slice.loc[sig_idx], index=False, header=True
    ):
        ws3.append(row)

    wb.save(args.output_name + ".xlsx")
    print("Done.")
    
    
def main():
    import argparse
    parser = argparse.ArgumentParser(description="CircTest wrapper")
    parser.add_argument("-d", "--dcc-data", required=True, dest="dcc_data")
    parser.add_argument("-c", "--condition-list", required=True, dest="condition_list")
    parser.add_argument("-l", "--condition-columns", required=True, dest="condition_columns")
    parser.add_argument("-g", "--groups", required=True, dest="groups")
    parser.add_argument("-o", "--output-name", required=True, dest="output_name")
    parser.add_argument("-r", "--replicates", type=int, default=2, dest="replicates")
    parser.add_argument("-f", "--max-fdr", type=float, default=0.05, dest="max_fdr")
    parser.add_argument("-m", "--max-plots", type=int, default=10, dest="max_plots")
    parser.add_argument("--filter-sample", type=int, default=1, dest="filter_sample")
    parser.add_argument("--filter-count", type=int, default=5, dest="filter_count")
    parser.add_argument("--output-label", default="", dest="output_label")
    parser.add_argument("--percent-filter", type=float, default=0.1, dest="percent_filter")
    parser.add_argument("--only-negative", action="store_true", dest="only_negative")
    parser.add_argument("--add-header", action="store_true", dest="add_header")
    parser.add_argument("--y-range", type=float, default=1.0, dest="y_range")
    parser.add_argument("--colour-mode", default="colour", dest="colour_mode")

    args = parser.parse_args()
    run_circ_test_wrapper(args)