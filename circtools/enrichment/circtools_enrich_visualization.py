#!/usr/bin/env python3

import sys
import warnings
import textwrap

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.cm import get_cmap

warnings.filterwarnings("ignore")

COLNAMES = [
    "RBP", "Annotation", "chr", "start", "stop", "strand",
    "p_val_circular", "raw_count_circ_rna", "observed_input_peaks_circ_rna",
    "length_circ_rna", "length_normalized_count_circ_rna",
    "number_of_features_intersecting_circ", "circ_rna_confidence_interval_0.05",
    "p_val_linear", "raw_count_host_gene", "observed_input_peaks_host_gene",
    "length_host_gene_without_circ_rna", "length_normalized_count_host_gene",
    "number_of_features_intersecting_linear", "host_gene_confidence_interval_0.05",
    "distance_normalized_counts",
]



def commapos(x, pos=None):
    """Formatter: absolute value with comma thousands separator."""
    return f"{abs(int(x)):,}"


def read_data(path, rbp_name=None):
    """Read a tab-separated enrichment file and assign column names.

    Handles two formats:
    - enrichment_check.py output: has a header row, no leading RBP column.
      Pass rbp_name= to inject a synthetic RBP column (e.g. output_filename).
    - Original multi-RBP TSV (standalone use): no header, first column is RBP.
      Leave rbp_name=None and the first column is used as-is.
    """
    if rbp_name is not None:
        # enrichment_check format: has a header, no RBP column
        df = pd.read_csv(path, sep="\t", header=0, comment="#")
        df.columns = COLNAMES[1 : df.shape[1] + 1]
        df.insert(0, "RBP", rbp_name)
    else:
        # standalone multi-RBP format: no header, RBP is the first column
        df = pd.read_csv(path, sep="\t", header=None, comment="#")
        df.columns = COLNAMES[: df.shape[1]]
    return df


def filter_data(df, pval, use_only_circ):
    """Apply p-value filter matching the R logic."""
    if use_only_circ:
        return df[df["p_val_circular"] < pval].copy()
    else:
        return df[(df["p_val_circular"] < pval) & (df["p_val_linear"] > pval)].copy()


def get_colormap(n, bw=False):
    """Return a list of n colours."""
    if bw:
        grays = np.linspace(0.0, 0.85, n)
        return [str(g) for g in grays]
    cmap = get_cmap("tab20" if n <= 20 else "hsv")
    return [cmap(i / max(n - 1, 1)) for i in range(n)]


def wrap_label(label, width=18):
    """Wrap a long string for axis tick labels."""
    return "\n".join(textwrap.wrap(str(label), width))


def plot_circrnas_per_rbp(pdf, df1, df2, label1, label2, max_rbps, pval, bw):
    """Bar chart: # unique circRNAs per RBP (isoforms collapsed)."""

    def count_rbps(df):
        tmp = df[["RBP", "Annotation"]].drop_duplicates()
        counts = tmp.groupby("RBP").size().sort_values(ascending=False)
        return counts.head(max_rbps)

    rbps1 = count_rbps(df1)
    all_rbps = list(rbps1.index)

    fig, ax = plt.subplots(figsize=(11.69, 8.2))

    x = np.arange(len(all_rbps))
    width = 0.4 if df2 is not None else 0.6
    colours = get_colormap(len(all_rbps), bw)

    if df2 is not None:
        rbps2 = count_rbps(df2).reindex(all_rbps, fill_value=0)
        vals1 = rbps1.values
        vals2 = rbps2.values
        ax.bar(x - width / 2, vals1, width=width, color=colours, edgecolor="black", linewidth=0.4, label=label1)
        ax.bar(x + width / 2, -vals2, width=width, color=colours, edgecolor="black", linewidth=0.4, label=label2)
        # sample labels
        ax.text(len(all_rbps) - 1.5, vals1.max() * 0.9, label1,
                ha="right", va="top", bbox=dict(fc="white", ec="gray", boxstyle="round,pad=0.3"))
        ax.text(len(all_rbps) - 1.5, -vals2.max() * 0.9, label2,
                ha="right", va="bottom", bbox=dict(fc="white", ec="gray", boxstyle="round,pad=0.3"))
    else:
        vals1 = rbps1.values
        ax.bar(x, vals1, width=width, color=colours, edgecolor="black", linewidth=0.4)
        ax.text(len(all_rbps) - 1.5, vals1.max() * 0.9, label1,
                ha="right", va="top", bbox=dict(fc="white", ec="gray", boxstyle="round,pad=0.3"))

    ax.set_xticks(x)
    ax.set_xticklabels([wrap_label(r) for r in all_rbps], rotation=90, ha="center", fontsize=10)
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(commapos))
    ax.set_xlabel("RNA binding protein", fontsize=12)
    ax.set_ylabel("Number of circular RNAs", fontsize=12)
    ax.set_title(f"{label1}: Number of circular RNAs per RBP\n"
                 f"Counting circRNAs (including different isoforms) with significantly enriched RBPs (p < {pval})",
                 fontsize=11, fontweight="bold")
    ax.axhline(0, color="black", linewidth=0.8)
    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


def plot_top_circs_by_rbp_hits(pdf, df1, df2, label1, label2, max_circs, max_rbps, pval, bw):
    """Bar chart: # distinct RBPs per circRNA."""

    def count_circs(df):
        counts = df.groupby("Annotation")["RBP"].nunique().sort_values(ascending=False)
        return counts

    circs1 = count_circs(df1)

    if df2 is not None:
        circs2 = count_circs(df2)
        combined = pd.DataFrame({"A": circs1, "B": circs2}).fillna(0)
        combined["Total"] = combined["A"] + combined["B"]
        combined = combined.sort_values("A", ascending=False)
    else:
        combined = pd.DataFrame({"A": circs1}).fillna(0)
        combined["Total"] = combined["A"]
        combined = combined.sort_values("A", ascending=False)

    combined = combined[combined["A"] > 0].head(max_circs)
    circs_list = list(combined.index)
    colours = get_colormap(len(circs_list), bw)

    fig, ax = plt.subplots(figsize=(11.69, 8.2))
    x = np.arange(len(circs_list))
    width = 0.6

    ax.bar(x, combined["A"].values, width=width, color=colours, edgecolor="black", linewidth=0.3)
    label_pos = combined["A"].max()
    ax.text(len(circs_list) - 4.5, label_pos * 0.95, label1,
            ha="right", va="top", bbox=dict(fc="white", ec="gray", boxstyle="round,pad=0.3"))

    if df2 is not None:
        ax.bar(x, -combined["B"].values, width=width, color=colours, edgecolor="black", linewidth=0.3)
        ax.text(len(circs_list) - 4.5, -label_pos * 0.95, label2,
                ha="right", va="bottom", bbox=dict(fc="white", ec="gray", boxstyle="round,pad=0.3"))

    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels([wrap_label(c, 14) for c in circs_list], rotation=45, ha="right", fontsize=8)
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(commapos))
    ax.set_xlabel("CircRNA", fontsize=12)
    ax.set_ylabel("Number of different RBPs found for circRNA", fontsize=12)
    ax.set_title(f"{label1}: Top CircRNAs (by RBP hits)\n"
                 f"plotting colour-coded RBPs per circRNA, ordered by number of distinct RBP hits (p < {pval})",
                 fontsize=11, fontweight="bold")
    num_circs_total = len(combined)
    num_rbps_total = df1["RBP"].nunique()
    ax.set_xlabel(
        f"CircRNA\n(based on data from {num_circs_total} circRNAs and {num_rbps_total} RBPs, "
        f"showing top {max_circs} circRNAs)",
        fontsize=10,
    )
    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)



def plot_polar_rbp_per_isoform(pdf, df1, df2, label1, label2, max_circs, max_rbps, pval, bw):
    """Pie charts showing RBP composition per circRNA isoform, matching R coord_polar style."""

    def top_circrnas(df, n):
        counts = df.groupby("Annotation").size().sort_values(ascending=False)
        return list(counts.head(n).index)

    def pie_chart(ax, mini_df, title, subtitle, pval, bw):
        """Draw a pie chart matching R geom_bar + coord_polar style."""
        mini_df = mini_df.sort_values("observed_input_peaks_circ_rna", ascending=False).head(max_rbps)
        if mini_df.empty:
            ax.set_visible(False)
            return

        values = mini_df["observed_input_peaks_circ_rna"].values.astype(float)
        labels = [f"{r}: {int(p)}" for r, p in
                  zip(mini_df["RBP"], mini_df["observed_input_peaks_circ_rna"])]
        colours = get_colormap(len(values), bw)

        # Draw pie — matches R coord_polar(geom_bar, width=1)
        wedges, texts = ax.pie(
            values,
            labels=labels,
            colors=colours,
            startangle=90,
            counterclock=False,
            wedgeprops=dict(edgecolor="white", linewidth=0.5),
            textprops=dict(fontsize=6),
        )

        ax.set_title(f"{title}\n{subtitle}", fontsize=7, fontweight="bold", pad=10)
        ax.set_xlabel(
            f"Top circRNAs enriched for RBP peaks compared to their host gene ( p < {pval} )",
            fontsize=6
        )

    # -------------------------------------------------------------------------
    sample_list = [(df1, label1)]
    if df2 is not None:
        sample_list.append((df2, label2))

    for df, label in sample_list:
        top_circs = top_circrnas(df, max_circs)

        for circ in top_circs:
            isoform_df = df[df["Annotation"] == circ][["start", "stop"]].drop_duplicates()
            n_iso = len(isoform_df)
            if n_iso == 0:
                continue

            ncols = min(n_iso, 2)
            nrows = int(np.ceil(n_iso / ncols))
            fig, axes = plt.subplots(nrows, ncols,
                                     figsize=(11.69, max(4, nrows * 4.5)))
            axes = np.array(axes).flatten()

            for idx, (_, row) in enumerate(isoform_df.iterrows()):
                iso_data = df[
                    (df["Annotation"] == circ) &
                    (df["start"] == row["start"]) &
                    (df["stop"] == row["stop"]) &
                    (df["observed_input_peaks_circ_rna"] > 0)
                ]
                chr_val = iso_data["chr"].iloc[0] if not iso_data.empty else "?"
                title = f"{label}:\nComposition of RBP landscape for circRNA {circ}"
                subtitle = (f"Isoform {idx+1}: Chromsome {chr_val}, "
                            f"{int(row['start']):,}->{int(row['stop']):,}")

                ax_pos = axes[idx].get_position()
                axes[idx].remove()
                ax_pie = fig.add_axes(ax_pos, aspect="equal")
                pie_chart(ax_pie, iso_data.copy(), title, subtitle, pval, bw)

            for j in range(n_iso, len(axes)):
                axes[j].set_visible(False)

            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)



def plot_accumulated_peaks_per_circ(pdf, df1, df2, label1, label2, max_circs, max_rbps, pval, bw):
    """Stacked bar chart: accumulated eCLIP peaks per circRNA, split by RBP."""

    def max_peaks_per_rbp_circ(df):
        """For each (RBP, Annotation) pair, keep the row with the max observed peaks."""
        idx = (
            df.groupby(["RBP", "Annotation"])["observed_input_peaks_circ_rna"]
            .idxmax()
        )
        return df.loc[idx, ["RBP", "Annotation", "observed_input_peaks_circ_rna"]].copy()

    data1 = max_peaks_per_rbp_circ(df1)

    if df2 is not None:
        data2 = max_peaks_per_rbp_circ(df2)
        merged = pd.merge(data1, data2, on=["Annotation", "RBP"], how="outer",
                          suffixes=("_A", "_B")).fillna(0)
        merged.columns = ["Annotation", "RBP", "A", "B"]
        merged["rbp_sum"] = merged["A"] + merged["B"]
    else:
        merged = data1.rename(columns={"observed_input_peaks_circ_rna": "A"})
        merged["B"] = 0
        merged["rbp_sum"] = merged["A"]

    # compute total per circRNA for sorting
    circ_totals = merged.groupby("Annotation")["A"].sum().sort_values(ascending=False)
    top_circs = list(circ_totals.head(max_circs).index)
    merged = merged[merged["Annotation"].isin(top_circs)]

    all_rbps = sorted(merged["RBP"].unique())
    colours_map = dict(zip(all_rbps, get_colormap(len(all_rbps), bw)))

    fig, ax = plt.subplots(figsize=(11.69, 8.2))
    circ_order = [c for c in top_circs if c in merged["Annotation"].values]

    bottom_pos = np.zeros(len(circ_order))
    bottom_neg = np.zeros(len(circ_order))

    x = np.arange(len(circ_order))
    circ_idx = {c: i for i, c in enumerate(circ_order)}

    for rbp in all_rbps:
        sub = merged[merged["RBP"] == rbp].set_index("Annotation")
        vals_A = np.array([sub.loc[c, "A"] if c in sub.index else 0 for c in circ_order])
        ax.bar(x, vals_A, bottom=bottom_pos, color=colours_map[rbp],
               edgecolor="black", linewidth=0.1, label=rbp)
        bottom_pos += vals_A

        if df2 is not None:
            vals_B = np.array([sub.loc[c, "B"] if c in sub.index else 0 for c in circ_order])
            ax.bar(x, -vals_B, bottom=bottom_neg, color=colours_map[rbp],
                   edgecolor="black", linewidth=0.1)
            bottom_neg -= vals_B

    label_pos = bottom_pos.max() / 2 if bottom_pos.max() > 0 else 1
    ax.text(len(circ_order) - 4.5, label_pos * 0.9, label1,
            ha="right", va="top", bbox=dict(fc="white", ec="gray", boxstyle="round,pad=0.3"))
    if df2 is not None:
        ax.text(len(circ_order) - 4.5, -label_pos * 0.9, label2,
                ha="right", va="bottom", bbox=dict(fc="white", ec="gray", boxstyle="round,pad=0.3"))

    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels([wrap_label(c, 14) for c in circ_order], rotation=45, ha="right", fontsize=8)
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(commapos))
    ax.set_ylabel("Number of enriched RBP binding sites detected in the CLIP data set", fontsize=10)

    num_circs = len(top_circs)
    num_rbps = len(all_rbps)
    ax.set_xlabel(
        f"CircRNA\n(based on data from {num_circs} circRNAs and {num_rbps} RBPs, "
        f"showing top {max_circs} circRNAs)",
        fontsize=10,
    )
    ax.set_title(
        f"{label1}: Assignment of RBP CLIP peaks to circRNAs\n"
        f"plotting colour-coded RBPs per circRNA, ordered by accumulated RBP eCLIP peaks (p < {pval})",
        fontsize=11, fontweight="bold",
    )

    # legend – two columns, small font
    handles, labels = ax.get_legend_handles_labels()
    seen = {}
    unique_handles, unique_labels = [], []
    for h, l in zip(handles, labels):
        if l not in seen:
            seen[l] = True
            unique_handles.append(h)
            unique_labels.append(l)
    ax.legend(unique_handles, unique_labels, ncol=2, fontsize=7,
              loc="upper right", framealpha=0.8)

    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)



def main():
    if len(sys.argv) < 9:
        print(__doc__)
        sys.exit(1)

    data_file_1        = sys.argv[1]
    pval               = float(sys.argv[2])
    max_circrnas       = int(sys.argv[3])
    max_rbps           = int(sys.argv[4])
    output_file        = sys.argv[5]
    label_sample_1     = sys.argv[6]
    colour_mode        = sys.argv[7]   # "colour" or "bw"
    use_only_circ_data = sys.argv[8]   # "True" or "False"
    label_sample_2     = sys.argv[9]  if len(sys.argv) > 9  else "Sample 2"
    data_file_2        = sys.argv[10] if len(sys.argv) > 10 else None

    if use_only_circ_data not in ("True", "False"):
        print("Please specify the data mode as 'True' or 'False'")
        sys.exit(1)
    if colour_mode not in ("colour", "bw"):
        print("Please specify the colour mode as 'colour' or 'bw'")
        sys.exit(1)

    use_only_circ = use_only_circ_data == "True"
    bw = colour_mode == "bw"

    print(f"Reading input file 1: {data_file_1}")
    df1 = filter_data(read_data(data_file_1), pval, use_only_circ)

    df2 = None
    if data_file_2:
        print(f"Reading input file 2: {data_file_2}")
        df2 = filter_data(read_data(data_file_2), pval, use_only_circ)

    print("Plotting data…")
    with PdfPages(output_file) as pdf:

        # page metadata
        from datetime import datetime
        d = pdf.infodict()
        d["Title"] = f"circtools RBP enrichment analysis - {label_sample_1}"
        d["Author"] = "circtools enrich_visualization.py"
        d["CreationDate"] = datetime.today()

        # --- Plot 1: circRNAs per RBP ---
        plot_circrnas_per_rbp(
            pdf, df1, df2, label_sample_1, label_sample_2,
            max_rbps, pval, bw,
        )

        # --- Plot 2: top circRNAs by # distinct RBPs ---
        plot_top_circs_by_rbp_hits(
            pdf, df1, df2, label_sample_1, label_sample_2,
            max_circrnas, max_rbps, pval, bw,
        )

        # --- Plot 3: pie chart RBP composition per isoform ---
        plot_polar_rbp_per_isoform(
                pdf, df1, df2, label_sample_1, label_sample_2,
                max_circrnas, max_rbps, pval, bw,
            )

        # --- Plot 4: accumulated eCLIP peaks per circRNA ---
        plot_accumulated_peaks_per_circ(
            pdf, df1, df2, label_sample_1, label_sample_2,
            max_circrnas, max_rbps, pval, bw,
        )

    print(f"Printing to {output_file}")
    print("Done")


if __name__ == "__main__":
    main()