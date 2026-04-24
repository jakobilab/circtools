#!/usr/bin/env python3

import sys
import os
import math
import json
import html as html_module

import pandas as pd


# ---------------------------------------------------------------------------
# Colour helpers  (mirrors R's RColorBrewer RdBu n=5, reversed)
# ---------------------------------------------------------------------------
_PALETTE = ["#2166AC", "#92C5DE", "#F7F7F7", "#F4A582", "#D6604D"]


def _symmetric_range(column: pd.Series, default: float):
    numeric = pd.to_numeric(column, errors="coerce").dropna()
    if numeric.empty:
        return default - 1, default + 1
    top    = numeric.max() - default
    bottom = default - numeric.min()
    span   = max(top, bottom, 1)
    return default - span, default + span


def _color_bar_html(value, vmin: float, vmax: float) -> str:
    try:
        fval = float(value)
    except (TypeError, ValueError):
        return html_module.escape(str(value))
    span = vmax - vmin if vmax != vmin else 1
    frac = max(0.0, min(1.0, (fval - vmin) / span))
    idx  = min(int(frac * len(_PALETTE)), len(_PALETTE) - 1)
    color = _PALETTE[idx]
    pct  = round(frac * 100)
    return (
        f'<span class="cbar" style="'
        f'background:linear-gradient(90deg,{color} {pct}%,transparent {pct}%)">'
        f'{fval:.2f}</span>'
    )


# ---------------------------------------------------------------------------
# HTML page template
# ---------------------------------------------------------------------------
_HTML_HEAD = """\
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>circtools &mdash; {experiment_name}</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&family=JetBrains+Mono:wght@400;500&display=swap" rel="stylesheet">
<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css">
<script src="https://ajax.googleapis.com/ajax/libs/jquery/3.3.1/jquery.min.js"></script>
<script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/js/bootstrap.min.js"></script>

<script>
/* SVG_DATA is injected below by the formatter */
</script>

<script>
$(document).ready(function() {{

    $(".has-svg").each(function() {{
        var key = $(this).data("svg-key").toString();
        var svg = (typeof SVG_DATA !== "undefined" && SVG_DATA[key])
                  ? SVG_DATA[key] : "<em>No diagram</em>";
        $(this).popover({{
            html: true, trigger: "hover", placement: "auto",
            container: "body",
            title: "circRNA structure &amp; primer positions",
            content: svg
        }});
    }});

    $('[data-toggle="popover"]').popover({{
        html: true, trigger: "hover", placement: "auto", container: "body"
    }});
}});
</script>

<style>
  /* ── Design tokens ──────────────────────────────────────────── */
  :root {{
    --bg:          #f0f4f8;
    --surface:     #ffffff;
    --border:      #e2e8f0;
    --border-dark: #cbd5e1;
    --text:        #1e293b;
    --text-muted:  #64748b;
    --accent:      #6366f1;
    --accent-dark: #4f46e5;
    --ok:          #059669;
    --ok-soft:     #d1fae5;
    --hit:         #dc2626;
    --hit-soft:    #fee2e2;
    --header-top:  #1e293b;
    --header-sub:  #334155;
    --row-alt:     #f8fafc;
    --row-hover:   #eef2ff;
    --radius:      10px;
    --shadow:      0 4px 24px rgba(0,0,0,.08);
    --font:        "Inter", system-ui, sans-serif;
    --mono:        "JetBrains Mono", "Fira Code", monospace;
  }}

  /* ── Reset / base ───────────────────────────────────────────── */
  *, *::before, *::after {{ box-sizing: border-box; margin: 0; padding: 0; }}

  body {{
    background: var(--bg);
    color: var(--text);
    font-family: var(--font);
    font-size: 13px;
    line-height: 1.5;
    padding: 32px 24px 64px;
  }}

  /* ── Page header ────────────────────────────────────────────── */
  .page-header {{
    margin-bottom: 28px;
  }}
  .page-header .badge-label {{
    display: inline-block;
    background: var(--accent);
    color: #fff;
    font-size: 10px;
    font-weight: 600;
    letter-spacing: .08em;
    text-transform: uppercase;
    padding: 3px 10px;
    border-radius: 20px;
    margin-bottom: 10px;
  }}
  .page-header h1 {{
    font-size: 26px;
    font-weight: 700;
    color: var(--text);
    letter-spacing: -.3px;
  }}
  .page-header h1 span {{
    color: var(--accent);
  }}
  .page-header p {{
    color: var(--text-muted);
    margin-top: 4px;
    font-size: 13px;
  }}

  /* ── Table card ─────────────────────────────────────────────── */
  .table-card {{
    background: var(--surface);
    border-radius: var(--radius);
    box-shadow: var(--shadow);
    overflow: hidden;
    border: 1px solid var(--border);
  }}

  /* ── Table ──────────────────────────────────────────────────── */
  table {{
    width: 100%;
    border-collapse: collapse;
    font-size: 12.5px;
  }}
  th, td {{
    padding: 9px 14px;
    border-bottom: 1px solid var(--border);
    vertical-align: middle;
    white-space: nowrap;
  }}
  td {{ color: var(--text); }}

  /* group header row */
  thead tr.group-header th {{
    background: var(--header-top);
    color: #fff;
    text-align: center;
    font-size: 11px;
    font-weight: 600;
    letter-spacing: .06em;
    text-transform: uppercase;
    padding: 10px 14px;
    border-bottom: none;
  }}
  thead tr.group-header th:first-child {{
    border-right: 1px solid #334155;
  }}

  /* column header row */
  thead tr.col-header th {{
    background: var(--header-sub);
    color: #94a3b8;
    font-size: 11px;
    font-weight: 500;
    letter-spacing: .04em;
    text-transform: uppercase;
    border-bottom: 2px solid var(--border-dark);
  }}

  /* body rows */
  tbody tr {{ transition: background .12s; }}
  tbody tr:nth-child(even) {{ background: var(--row-alt); }}
  tbody tr:hover {{ background: var(--row-hover); }}
  tbody tr:last-child td {{ border-bottom: none; }}

  /* ── Colour-bar spans ───────────────────────────────────────── */
  .cbar {{
    display: inline-block;
    width: 100%;
    padding: 3px 7px;
    border-radius: 4px;
    font-variant-numeric: tabular-nums;
    font-size: 12px;
  }}

  /* ── Primer sequence badge ──────────────────────────────────── */
  .primer-btn {{
    display: inline-block;
    padding: 4px 10px;
    border-radius: 6px;
    font-family: var(--mono);
    font-size: 11.5px;
    font-weight: 500;
    letter-spacing: .01em;
    cursor: pointer;
    transition: filter .15s, transform .1s;
    border: none;
    outline: none;
  }}
  .primer-btn:hover {{ filter: brightness(1.12); transform: translateY(-1px); }}
  .primer-btn.ok  {{ background: var(--ok-soft);  color: var(--ok);  border: 1px solid #a7f3d0; }}
  .primer-btn.hit {{ background: var(--hit-soft); color: var(--hit); border: 1px solid #fca5a5; }}

  /* ── BLAST count badge ──────────────────────────────────────── */
  .blast-badge {{
    display: inline-flex;
    align-items: center;
    justify-content: center;
    min-width: 26px;
    height: 22px;
    padding: 0 8px;
    border-radius: 20px;
    font-size: 11px;
    font-weight: 600;
    cursor: pointer;
    transition: filter .15s;
  }}
  .blast-badge:hover {{ filter: brightness(1.1); }}
  .blast-badge.ok  {{ background: var(--ok-soft);  color: var(--ok); }}
  .blast-badge.hit {{ background: var(--hit-soft); color: var(--hit); }}

  /* ── Strand badge ───────────────────────────────────────────── */
  .strand {{
    display: inline-block;
    padding: 2px 8px;
    border-radius: 4px;
    font-size: 11px;
    font-weight: 600;
    background: #e0e7ff;
    color: var(--accent-dark);
  }}

  /* ── Bootstrap popover overrides ───────────────────────────── */
  .popover {{
    max-width: 560px;
    border-radius: 12px;
    border: 1px solid var(--border);
    box-shadow: 0 8px 32px rgba(0,0,0,.14);
    font-family: var(--font);
  }}
  .popover-title {{
    background: var(--header-top);
    color: #fff;
    font-size: 12px;
    font-weight: 600;
    letter-spacing: .04em;
    text-transform: uppercase;
    border-radius: 11px 11px 0 0;
    padding: 10px 16px;
    border-bottom: none;
  }}
  .popover-content {{
    padding: 12px;
    background: var(--surface);
    border-radius: 0 0 11px 11px;
  }}
  .popover-content svg {{
    width: 500px !important;
    height: 500px !important;
    display: block;
    border-radius: 6px;
  }}

  /* ── Scrollable table wrapper ───────────────────────────────── */
  .table-scroll {{ overflow-x: auto; }}
</style>
</head>
<body>

<div class="page-header">
  <div class="badge-label">circtools</div>
  <h1>Primer Design &mdash; <span>{experiment_name}</span></h1>
  <p>Hover over a primer sequence to view the circRNA diagram &bull; Hover over a BLAST badge to see off-target hits</p>
</div>

<div class="table-card">
<div class="table-scroll">
"""

_HTML_FOOT = """\
</div>
</div>

</body>
</html>"""


# ---------------------------------------------------------------------------
# Cell renderers
# ---------------------------------------------------------------------------

def _primer_cell(sequence: str, blast_count: int, svg_key: str, high: int = 0) -> str:
    cls     = "hit" if blast_count > high else "ok"
    seq_esc = html_module.escape(str(sequence))
    if svg_key:
        return (
            f'<span class="primer-btn {cls} has-svg" '
            f'data-svg-key="{html_module.escape(svg_key)}">'
            f'{seq_esc}</span>'
        )
    return f'<span class="primer-btn {cls}">{seq_esc}</span>'


def _blast_badge(count: int, blast_fmt: str, high: int = 0) -> str:
    cls          = "hit" if count > high else "ok"
    safe_content = blast_fmt.replace('"', "&quot;")
    return (
        f'<span class="blast-badge {cls}" '
        f'data-toggle="popover" '
        f'data-title="BLAST hits" '
        f'data-content="{safe_content}">'
        f'{count}</span>'
    )




def main():
    if len(sys.argv) < 3:
        sys.stderr.write(
            "Usage: circtools_primex_formatter.py "
            "<data_file> <experiment_name> [output_dir]\n"
        )
        sys.exit(1)

    data_file = sys.argv[1]
    exp_name  = sys.argv[2]

    if len(sys.argv) >= 4 and os.path.isdir(sys.argv[3]):
        output_dir = os.path.normpath(sys.argv[3])
    else:
        output_dir = os.path.dirname(os.path.abspath(data_file))

    os.makedirs(output_dir, exist_ok=True)

    csv_out  = os.path.join(output_dir, f"{exp_name}_primers.csv")
    xlsx_out = os.path.join(output_dir, f"{exp_name}_primers.xlsx")

    df = pd.read_csv(data_file, header=None, sep="\t",
                     keep_default_na=False, dtype=str)

    # Guarantee 18 columns
    while df.shape[1] < 18:
        df[df.shape[1]] = ""

    df.columns = list(range(df.shape[1]))

    # Restore whitespace escaped for TSV transport
    df[17] = df[17].str.replace("&#10;", "\n", regex=False).str.replace("&#9;", "\t", regex=False)

    df = df.rename(columns={
        0:  "Annotation", 1: "Chr",      2: "Start",    3: "Stop",
        4:  "Strand",     5: "idx",
        6:  "Left_seq",   7: "Right_seq",
        8:  "Left_pos",   9: "Right_pos",
        10: "TM_left",   11: "TM_right", 12: "GC_left", 13: "GC_right",
        14: "Product_size",
        15: "BLAST_left", 16: "BLAST_right",
        17: "SVG",
    })

    # Numeric conversions
    for col in ("TM_left", "TM_right", "GC_left", "GC_right", "Product_size"):
        df[col] = pd.to_numeric(df[col], errors="coerce")

    # ------------------------------------------------------------------
    # Clean location columns (lone "0" means FASTA-mode, no genomic coords)
    # ------------------------------------------------------------------
    df["Chr"]    = df["Chr"].str.replace(r"^\s*0\s*$", "",  regex=True)
    df["Start"]  = df["Start"].str.replace(r"^\s*0\s*$", "", regex=True)
    df["Stop"]   = df.apply(lambda r: "" if r["Start"] == "" else r["Stop"], axis=1)
    df["Strand"] = df["Strand"].str.replace(r"^\s*0\s*$", "", regex=True)


    HIGH = 0

    def _count(s):
        s = str(s).strip()
        if s in ("", "NA", "No hits", "Not blasted, no primer pair found"):
            return 0
        return s.count(";") + 1

    df["blast_left_n"]    = df["BLAST_left"].apply(_count)
    df["blast_right_n"]   = df["BLAST_right"].apply(_count)
    df["BLAST_left_fmt"]  = df["BLAST_left"].str.replace(";", "<br><br>", regex=False)
    df["BLAST_right_fmt"] = df["BLAST_right"].str.replace(";", "<br><br>", regex=False)

    svg_lookup = {}
    df["_svg_key"] = ""

    for i, row in df.iterrows():
        svg = str(row["SVG"]).strip()
        if svg:
            key = str(i)
            svg_lookup[key] = svg
            df.at[i, "_svg_key"] = key

    # ------------------------------------------------------------------
    # Export CSV / Excel (without the SVG column)
    # ------------------------------------------------------------------
    export_cols = [
        "Annotation", "Chr", "Start", "Stop", "Strand",
        "Left_seq", "Right_seq",
        "TM_left", "TM_right", "GC_left", "GC_right", "Product_size",
        "BLAST_left", "BLAST_right",
    ]
    df[export_cols].to_csv(csv_out, index=False)
    df[export_cols].to_excel(xlsx_out, sheet_name="PrimerResults", index=False)
    sys.stderr.write(f"Wrote: {csv_out}\n      {xlsx_out}\n")

    tm_min, tm_max = _symmetric_range(df["TM_left"],      60.0)
    gc_min, gc_max = _symmetric_range(df["GC_left"],      50.0)
    pr_min, pr_max = _symmetric_range(df["Product_size"], 100.0)

    thead = (
        '<thead>'
        '<tr class="group-header">'
          '<th colspan="5">Input circRNAs</th>'
          '<th colspan="9">Designed Primers</th>'
        '</tr>'
        '<tr class="col-header">'
          '<th>Annotation</th><th>Chr</th><th>Start</th><th>Stop</th><th>Strand</th>'
          '<th>TM&nbsp;fwd</th><th>TM&nbsp;rev</th>'
          '<th>GC%&nbsp;fwd</th><th>GC%&nbsp;rev</th>'
          '<th>Product&nbsp;size</th>'
          '<th>Forward primer</th><th>BLAST</th>'
          '<th>Reverse primer</th><th>BLAST</th>'
        '</tr>'
        '</thead>'
    )

    rows = []
    for _, row in df.iterrows():

        def safe(v):
            try:
                return "" if (v is None or (isinstance(v, float) and math.isnan(v))) else v
            except Exception:
                return v

        svg_key    = str(row["_svg_key"])
        left_n     = int(row["blast_left_n"])
        right_n    = int(row["blast_right_n"])
        blast_lfmt = str(row["BLAST_left_fmt"])
        blast_rfmt = str(row["BLAST_right_fmt"])

        strand_val = str(safe(row["Strand"]))
        strand_cell = f'<span class="strand">{html_module.escape(strand_val)}</span>' if strand_val else ""

        cells = [
            html_module.escape(str(safe(row["Annotation"]))),
            html_module.escape(str(safe(row["Chr"]))),
            html_module.escape(str(safe(row["Start"]))),
            html_module.escape(str(safe(row["Stop"]))),
            strand_cell,

            _color_bar_html(row["TM_left"],     tm_min, tm_max),
            _color_bar_html(row["TM_right"],    tm_min, tm_max),
            _color_bar_html(row["GC_left"],     gc_min, gc_max),
            _color_bar_html(row["GC_right"],    gc_min, gc_max),
            _color_bar_html(row["Product_size"], pr_min, pr_max),

            _primer_cell(row["Left_seq"],  left_n,  svg_key, HIGH),
            _blast_badge(left_n,  blast_lfmt, HIGH),
            _primer_cell(row["Right_seq"], right_n, svg_key, HIGH),
            _blast_badge(right_n, blast_rfmt, HIGH),
        ]

        rows.append(
            "<tr>" + "".join(f"<td>{c}</td>" for c in cells) + "</tr>"
        )

    table = (
        '<table class="table table-bordered table-hover">'
        + thead
        + "<tbody>\n" + "\n".join(rows) + "\n</tbody>"
        + "</table>"
    )

    svg_js = f"var SVG_DATA = {json.dumps(svg_lookup)};"

    html_out = _HTML_HEAD.format(
        experiment_name=html_module.escape(exp_name)
    ).replace(
        "/* SVG_DATA is injected below by the formatter */",
        svg_js
    )

    print(html_out + table + _HTML_FOOT)


if __name__ == "__main__":
    main()