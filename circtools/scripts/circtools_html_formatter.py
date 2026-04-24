#!/usr/bin/env python3


import sys
import os
import math
import json
import html as html_module

import pandas as pd


# ── Colour palette (RdBu 5-colour, reversed so blue=low) ─────────────────
_PALETTE_DIVERG = ["#2166AC", "#92C5DE", "#F7F7F7", "#F4A582", "#D6604D"]

# ── Threshold colours for siRNA (above/below threshold) ──────────────────
_COL_GOOD = "#d1fae5"
_COL_BAD  = "#fee2e2"
_TXT_GOOD = "#456d62"
_TXT_BAD  = "#991b1b"

# ── Padlock ligation-junction colours ────────────────────────────────────
_COL_PREFERRED = "#bbf7d0"
_COL_NEUTRAL   = "#fef9c3"
_TXT_PREFERRED = "#166534"
_TXT_NEUTRAL   = "#713f12"


# ─────────────────────────────────────────────────────────────────────────
# Shared CSS / JS head (same modern design as primex formatter)
# ─────────────────────────────────────────────────────────────────────────
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
/* SVG_DATA injected by formatter */
{svg_data_js}
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
            title: "{popover_title}",
            content: svg
        }});
    }});

    $('[data-toggle="popover"]').popover({{
        html: true, trigger: "hover", placement: "auto", container: "body"
    }});

    // siRNA SVG modal
    $(document).on("click", ".svg-modal-btn", function() {{
        var key = $(this).data("svg-key").toString();
        var svg = (typeof SVG_DATA !== "undefined" && SVG_DATA[key])
                  ? SVG_DATA[key] : "<em>No diagram available</em>";
        var title = $(this).data("svg-title") || "circRNA diagram";
        $("#svgModalLabel").text(title);
        $("#svgModalBody").html(svg);
        $("#svgModal").modal("show");
    }});
}});
</script>

<style>
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
    --warn:        #d97706;
    --warn-soft:   #fef3c7;
    --header-top:  #1e293b;
    --header-sub:  #334155;
    --row-alt:     #f8fafc;
    --row-hover:   #eef2ff;
    --radius:      10px;
    --shadow:      0 4px 24px rgba(0,0,0,.08);
    --font:        "Inter", system-ui, sans-serif;
    --mono:        "JetBrains Mono", "Fira Code", monospace;
  }}

  *, *::before, *::after {{ box-sizing: border-box; margin: 0; padding: 0; }}

  body {{
    background: var(--bg);
    color: var(--text);
    font-family: var(--font);
    font-size: 16px;
    line-height: 1.7;
    padding: 48px 40px 80px;
  }}

  .page-header {{ margin-bottom: 40px; }}
  .page-header .badge-label {{
    display: inline-block;
    background: var(--accent); color: #fff;
    font-size: 13px; font-weight: 600;
    letter-spacing: .07em; text-transform: uppercase;
    padding: 5px 14px; border-radius: 20px; margin-bottom: 14px;
  }}
  .page-header h1 {{
    font-size: 34px; font-weight: 700;
    color: var(--text); letter-spacing: -.5px; line-height: 1.2;
  }}
  .page-header h1 span {{ color: var(--accent); }}
  .page-header p {{ color: var(--text-muted); margin-top: 8px; font-size: 16px; line-height: 1.6; }}

  .table-card {{
    background: var(--surface);
    border-radius: var(--radius);
    box-shadow: var(--shadow);
    overflow: hidden;
    border: 1px solid var(--border);
  }}
  .table-scroll {{ overflow-x: auto; }}

  table {{ width: 100%; border-collapse: collapse; font-size: 15px; }}
  th, td {{
    padding: 14px 20px;
    border-bottom: 1px solid var(--border);
    vertical-align: middle;
    white-space: nowrap;
  }}
  td {{ color: var(--text); font-size: 15px; }}

  thead tr.group-header th {{
    background: var(--header-top); color: #fff;
    text-align: center; font-size: 14px; font-weight: 600;
    letter-spacing: .07em; text-transform: uppercase;
    padding: 16px 20px; border-bottom: none;
  }}
  thead tr.group-header th:first-child {{ border-right: 1px solid #334155; }}
  thead tr.col-header th {{
    background: var(--header-sub); color: #94a3b8;
    font-size: 13px; font-weight: 600;
    letter-spacing: .05em; text-transform: uppercase;
    border-bottom: 2px solid var(--border-dark);
    padding: 13px 20px;
  }}

  tbody tr {{ transition: background .12s; }}
  tbody tr:nth-child(even) {{ background: var(--row-alt); }}
  tbody tr:hover {{ background: var(--row-hover); }}
  tbody tr:last-child td {{ border-bottom: none; }}

  /* colour bar */
  .cbar {{
    display: inline-block; width: 100%;
    padding: 6px 12px; border-radius: 6px;
    font-variant-numeric: tabular-nums; font-size: 15px;
    font-weight: 500;
  }}

  /* sequence badges */
  .seq-badge {{
    display: inline-block;
    padding: 6px 14px; border-radius: 8px;
    font-family: var(--mono); font-size: 14px; font-weight: 500;
    cursor: pointer; letter-spacing: .02em;
    transition: filter .15s, transform .1s;
  }}
  .seq-badge:hover {{ filter: brightness(1.1); transform: translateY(-1px); }}
  .seq-badge.ok  {{ background: var(--ok-soft);  color: var(--ok);  border: 1px solid #a7f3d0; }}
  .seq-badge.hit {{ background: var(--hit-soft); color: var(--hit); border: 1px solid #fca5a5; }}

  /* count badge */
  .count-badge {{
    display: inline-flex; align-items: center; justify-content: center;
    min-width: 36px; height: 30px; padding: 0 12px;
    border-radius: 20px; font-size: 14px; font-weight: 700; cursor: pointer;
    transition: filter .15s;
  }}
  .count-badge:hover {{ filter: brightness(1.1); }}
  .count-badge.ok  {{ background: var(--ok-soft);  color: var(--ok); }}
  .count-badge.hit {{ background: var(--hit-soft); color: var(--hit); }}

  /* threshold badges (siRNA scores) */
  .thresh-badge {{
    display: inline-block; padding: 6px 14px; border-radius: 8px;
    font-size: 15px; font-weight: 600;
    font-variant-numeric: tabular-nums;
  }}
  .thresh-badge.good {{ background: var(--ok-soft);   color: var(--ok); }}
  .thresh-badge.bad  {{ background: var(--hit-soft);  color: var(--hit); }}
  .thresh-badge.warn {{ background: var(--warn-soft); color: var(--warn); }}

  /* ligation junction badge */
  .junc-badge {{
    display: inline-block; padding: 6px 14px; border-radius: 8px;
    font-size: 14px; font-weight: 600;
  }}
  .junc-badge.preferred {{ background: {col_preferred}; color: {txt_preferred}; }}
  .junc-badge.neutral   {{ background: {col_neutral};   color: {txt_neutral}; }}

  /* strand badge */
  .strand {{
    display: inline-block; padding: 4px 12px; border-radius: 6px;
    font-size: 14px; font-weight: 600;
    background: #e0e7ff; color: var(--accent-dark);
    letter-spacing: .03em;
  }}

  /* popovers */
  .popover {{
    max-width: 600px; border-radius: 14px;
    border: 1px solid var(--border);
    box-shadow: 0 8px 32px rgba(0,0,0,.14);
    font-family: var(--font);
  }}
  .popover-title {{
    background: var(--header-top); color: #fff;
    font-size: 14px; font-weight: 600;
    letter-spacing: .05em; text-transform: uppercase;
    border-radius: 13px 13px 0 0;
    padding: 14px 20px; border-bottom: none;
  }}
  .popover-content {{
    padding: 16px; background: var(--surface);
    border-radius: 0 0 13px 13px;
    font-size: 14px; line-height: 1.75;
  }}
  .popover-content svg {{
    width: 500px !important; height: 500px !important;
    display: block; border-radius: 6px;
  }}

  /* magnify button (siRNA) */
  .svg-modal-btn {{
    display: inline-flex; align-items: center; justify-content: center;
    width: 28px; height: 28px; margin-left: 8px;
    border-radius: 50%; border: 1px solid var(--border-dark);
    background: var(--surface); color: var(--text-muted);
    font-size: 15px; cursor: pointer; vertical-align: middle;
    transition: background .15s, color .15s, transform .1s;
    flex-shrink: 0;
  }}
  .svg-modal-btn:hover {{
    background: var(--accent); color: #fff;
    border-color: var(--accent); transform: scale(1.15);
  }}

  /* SVG modal */
  #svgModal .modal-dialog {{ width: 720px; max-width: 95vw; }}
  #svgModal .modal-header {{
    background: var(--header-top); color: #fff;
    border-radius: 6px 6px 0 0; padding: 18px 24px;
  }}
  #svgModal .modal-title {{
    font-family: var(--font); font-size: 15px;
    font-weight: 600; letter-spacing: .05em; text-transform: uppercase;
  }}
  #svgModal .modal-header .close {{
    color: #fff; opacity: .7; font-size: 24px; margin-top: -2px;
  }}
  #svgModal .modal-header .close:hover {{ opacity: 1; }}
  #svgModal .modal-body {{
    padding: 24px; display: flex;
    align-items: center; justify-content: center;
    background: var(--bg);
  }}
  #svgModal .modal-body svg {{
    width: 640px !important; height: 640px !important;
    display: block; border-radius: 10px;
    background: var(--surface);
    box-shadow: 0 2px 12px rgba(0,0,0,.08);
  }}
  #svgModal .modal-footer {{
    border-top: 1px solid var(--border); padding: 14px 24px;
  }}
</style>
</head>
<body>
<div class="page-header">
  <div class="badge-label">circtools</div>
  <h1>{page_title_prefix} &mdash; <span>{experiment_name}</span></h1>
  <p>{page_subtitle}</p>
</div>
<div class="table-card">
<div class="table-scroll">
"""

_HTML_FOOT = """\
</div></div>

<!-- SVG diagram modal (used by siRNA view buttons) -->
<div class="modal fade" id="svgModal" tabindex="-1" role="dialog" aria-labelledby="svgModalLabel">
  <div class="modal-dialog" role="document">
    <div class="modal-content">
      <div class="modal-header">
        <button type="button" class="close" data-dismiss="modal" aria-label="Close">
          <span aria-hidden="true">&times;</span>
        </button>
        <h4 class="modal-title" id="svgModalLabel">circRNA diagram</h4>
      </div>
      <div class="modal-body" id="svgModalBody"></div>
      <div class="modal-footer">
      </div>
    </div>
  </div>
</div>

</body></html>"""


def _symmetric_range(column: pd.Series, default: float):
    numeric = pd.to_numeric(column, errors="coerce").dropna()
    if numeric.empty:
        return default - 1, default + 1
    span = max(numeric.max() - default, default - numeric.min(), 1)
    return default - span, default + span


def _color_bar(value, vmin: float, vmax: float) -> str:
    try:
        fval = float(value)
    except (TypeError, ValueError):
        return html_module.escape(str(value))
    span = vmax - vmin if vmax != vmin else 1
    frac = max(0.0, min(1.0, (fval - vmin) / span))
    idx  = min(int(frac * len(_PALETTE_DIVERG)), len(_PALETTE_DIVERG) - 1)
    color = _PALETTE_DIVERG[idx]
    pct   = round(frac * 100)
    return (
        f'<span class="cbar" style="'
        f'background:linear-gradient(90deg,{color} {pct}%,transparent {pct}%)">'
        f'{fval:.2f}</span>'
    )


def _seq_badge(sequence: str, blast_count: int, svg_key: str = "", high: int = 0) -> str:
    cls     = "hit" if blast_count > high else "ok"
    seq_esc = html_module.escape(str(sequence))
    if svg_key:
        return (
            f'<span class="seq-badge {cls} has-svg" '
            f'data-svg-key="{html_module.escape(svg_key)}">'
            f'{seq_esc}</span>'
        )
    return f'<span class="seq-badge {cls}">{seq_esc}</span>'


def _count_badge(count: int, blast_fmt: str, high: int = 0) -> str:
    cls  = "hit" if count > high else "ok"
    safe = blast_fmt.replace('"', "&quot;")
    return (
        f'<span class="count-badge {cls}" '
        f'data-toggle="popover" data-title="BLAST hits" '
        f'data-content="{safe}">{count}</span>'
    )


def _thresh_badge(value, threshold: float, higher_is_bad: bool = True) -> str:
    try:
        fval = float(value)
    except (TypeError, ValueError):
        return html_module.escape(str(value))
    bad = (fval > threshold) if higher_is_bad else (fval < threshold)
    cls = "bad" if bad else "good"
    return f'<span class="thresh-badge {cls}">{fval:.2f}</span>'


def _junc_badge(value: str) -> str:
    cls = "preferred" if str(value).lower() == "preferred" else "neutral"
    return f'<span class="junc-badge {cls}">{html_module.escape(str(value))}</span>'


def _safe(v):
    try:
        return "" if (v is None or (isinstance(v, float) and math.isnan(v))) else v
    except Exception:
        return v


def _strand_cell(v) -> str:
    s = str(_safe(v))
    return f'<span class="strand">{html_module.escape(s)}</span>' if s else ""


def _count_blast(s: str) -> int:
    s = str(s).strip()
    if s in ("", "NA", "No hits", "Not blasted, no primer pair found", "N/A"):
        return 0
    return s.count(";") + 1


def _clean_location(df: pd.DataFrame) -> pd.DataFrame:
    for col in ("Chr", "Start", "Strand"):
        df[col] = df[col].astype(str).str.replace(r"^\s*0\s*$", "", regex=True)
    df["Stop"] = df.apply(lambda r: "" if r["Start"] == "" else r["Stop"], axis=1)
    return df


def _build_svg_lookup(df: pd.DataFrame, svg_col: str = "SVG") -> tuple[dict, pd.Series]:
    lookup: dict[str, str] = {}
    keys = []
    for i, row in df.iterrows():
        svg = str(row.get(svg_col, "")).strip()
        if svg:
            key = str(i)
            lookup[key] = svg
            keys.append(key)
        else:
            keys.append("")
    return lookup, pd.Series(keys, index=df.index)


def _thead(*groups) -> str:

    group_pairs, col_names = groups
    group_cells = "".join(f'<th colspan="{n}">{g}</th>' for g, n in group_pairs)
    col_cells   = "".join(f"<th>{c}</th>" for c in col_names)
    return (
        "<thead>"
        f'<tr class="group-header">{group_cells}</tr>'
        f'<tr class="col-header">{col_cells}</tr>'
        "</thead>"
    )


def _table(thead_html: str, rows: list[str]) -> str:
    return (
        '<table class="table table-bordered table-hover">'
        + thead_html
        + "<tbody>\n" + "\n".join(rows) + "\n</tbody>"
        + "</table>"
    )


def _render_page(mode: str, exp_name: str, table_html: str,
                 svg_lookup: dict) -> str:
    cfg = {
        "primex":  ("Primer Design",       "circRNA structure &amp; primer positions",
                    "Hover over a primer to view the circRNA diagram · Hover a BLAST badge to see off-target hits"),
        "padlock": ("Padlock Probe Design", "circRNA structure &amp; probe positions",
                    "Hover over a probe sequence to view the circRNA diagram · Hover a BLAST badge to see off-target hits"),
        "sirna":   ("siRNA Design",         "circRNA structure &amp; siRNA position",
                    "Click magnifying glass next to siRNA sequence to view its position on the circRNA"),
    }
    prefix, popover_title, subtitle = cfg.get(mode, ("Results", "Diagram", ""))
    svg_js = f"var SVG_DATA = {json.dumps(svg_lookup)};"
    head = _HTML_HEAD.format(
        experiment_name  = html_module.escape(exp_name),
        svg_data_js      = svg_js,
        popover_title    = popover_title,
        page_title_prefix= prefix,
        page_subtitle    = subtitle,
        col_preferred    = _COL_PREFERRED,
        txt_preferred    = _TXT_PREFERRED,
        col_neutral      = _COL_NEUTRAL,
        txt_neutral      = _TXT_NEUTRAL,
    )
    return head + table_html + _HTML_FOOT

#  formatter class


class CirctoolsHTMLFormatter:


    def __init__(self, mode: str):
        if mode not in ("primex", "padlock", "sirna"):
            raise ValueError(f"Unknown mode '{mode}'. Use primex, padlock or sirna.")
        self.mode = mode



    def format_file(self, data_file: str, exp_name: str,
                    output_dir: str = "", extra_flag: str = "") -> str:

        if self.mode == "primex":
            return self._format_primex(data_file, exp_name, output_dir)
        elif self.mode == "padlock":
            return self._format_padlock(data_file, exp_name, output_dir, extra_flag)
        elif self.mode == "sirna":
            blast_was_run = extra_flag.strip().lower() in ("true", "1", "yes")
            return self._format_sirna(data_file, exp_name, output_dir, blast_was_run)


    def _format_primex(self, data_file: str, exp_name: str,
                       output_dir: str) -> str:
        df = pd.read_csv(data_file, header=None, sep="\t",
                         keep_default_na=False, dtype=str)
        while df.shape[1] < 18:
            df[df.shape[1]] = ""
        df.columns = list(range(df.shape[1]))
        df[17] = df[17].str.replace("&#10;", "\n", regex=False)\
                        .str.replace("&#9;",  "\t", regex=False)
        df = df.rename(columns={
            0:"Annotation",1:"Chr",2:"Start",3:"Stop",4:"Strand",5:"idx",
            6:"Left_seq",7:"Right_seq",8:"Left_pos",9:"Right_pos",
            10:"TM_left",11:"TM_right",12:"GC_left",13:"GC_right",
            14:"Product_size",15:"BLAST_left",16:"BLAST_right",17:"SVG",
        })
        for c in ("TM_left","TM_right","GC_left","GC_right","Product_size"):
            df[c] = pd.to_numeric(df[c], errors="coerce")
        df = _clean_location(df)

        df["blast_left_n"]    = df["BLAST_left"].apply(_count_blast)
        df["blast_right_n"]   = df["BLAST_right"].apply(_count_blast)
        df["BLAST_left_fmt"]  = df["BLAST_left"].str.replace(";","<br><br>",regex=False)
        df["BLAST_right_fmt"] = df["BLAST_right"].str.replace(";","<br><br>",regex=False)

        svg_lookup, svg_keys = _build_svg_lookup(df)
        df["_svg_key"] = svg_keys

        self._export(df,
                     ["Annotation","Chr","Start","Stop","Strand",
                      "Left_seq","Right_seq","TM_left","TM_right",
                      "GC_left","GC_right","Product_size",
                      "BLAST_left","BLAST_right"],
                     exp_name, output_dir)

        tm_min,tm_max = _symmetric_range(df["TM_left"],      60.0)
        gc_min,gc_max = _symmetric_range(df["GC_left"],      50.0)
        pr_min,pr_max = _symmetric_range(df["Product_size"], 100.0)

        thead_html = _thead(
            [("Input circRNAs", 5), ("Designed Primers", 9)],
            ["Annotation","Chr","Start","Stop","Strand",
             "TM&nbsp;fwd","TM&nbsp;rev","GC%&nbsp;fwd","GC%&nbsp;rev",
             "Product&nbsp;size","Forward primer","BLAST","Reverse primer","BLAST"],
        )

        rows = []
        for _, row in df.iterrows():
            svg_key = str(row["_svg_key"])
            cells = [
                html_module.escape(str(_safe(row["Annotation"]))),
                html_module.escape(str(_safe(row["Chr"]))),
                html_module.escape(str(_safe(row["Start"]))),
                html_module.escape(str(_safe(row["Stop"]))),
                _strand_cell(row["Strand"]),
                _color_bar(row["TM_left"],     tm_min, tm_max),
                _color_bar(row["TM_right"],    tm_min, tm_max),
                _color_bar(row["GC_left"],     gc_min, gc_max),
                _color_bar(row["GC_right"],    gc_min, gc_max),
                _color_bar(row["Product_size"],pr_min, pr_max),
                _seq_badge(row["Left_seq"],  int(row["blast_left_n"]),  svg_key),
                _count_badge(int(row["blast_left_n"]),  str(row["BLAST_left_fmt"])),
                _seq_badge(row["Right_seq"], int(row["blast_right_n"]), svg_key),
                _count_badge(int(row["blast_right_n"]), str(row["BLAST_right_fmt"])),
            ]
            rows.append("<tr>" + "".join(f"<td>{c}</td>" for c in cells) + "</tr>")

        return _render_page("primex", exp_name,
                            _table(thead_html, rows), svg_lookup)


    def _format_padlock(self, data_file: str, exp_name: str,
                        output_dir: str, no_svg_flag: str) -> str:
        df = pd.read_csv(data_file, header=None, sep="\t",
                         keep_default_na=False, dtype=str)
        while df.shape[1] < 16:
            df[df.shape[1]] = ""
        df.columns = list(range(df.shape[1]))

        # Restore whitespace escaped for TSV transport
        df[15] = df[15].str.replace("&#10;", "\n", regex=False)\
                        .str.replace("&#9;",  "\t", regex=False)

        df = df.rename(columns={
            0:"Annotation",1:"Chr",2:"Start",3:"Stop",4:"Strand",
            5:"RBD5",6:"RBD3",
            7:"TM_left",8:"TM_right",9:"TM_Full",
            10:"GC_left",11:"GC_right",
            12:"Ligation_Junction",
            13:"BLAST_left",14:"BLAST_right",
            15:"SVG",
        })
        for c in ("TM_left","TM_right","TM_Full","GC_left","GC_right"):
            df[c] = pd.to_numeric(df[c], errors="coerce")
        df = _clean_location(df)

        df["blast_left_n"]    = df["BLAST_left"].apply(_count_blast)
        df["blast_right_n"]   = df["BLAST_right"].apply(_count_blast)
        df["BLAST_left_fmt"]  = df["BLAST_left"].str.replace(";","<br><br>",regex=False)
        df["BLAST_right_fmt"] = df["BLAST_right"].str.replace(";","<br><br>",regex=False)

        svg_lookup, svg_keys = _build_svg_lookup(df)
        df["_svg_key"] = svg_keys

        self._export(df,
                     ["Annotation","Chr","Start","Stop","Strand",
                      "RBD5","RBD3","TM_left","TM_right","TM_Full",
                      "GC_left","GC_right","Ligation_Junction",
                      "BLAST_left","BLAST_right"],
                     exp_name, output_dir)

        # colour ranges fixed per R script
        tm_rbd_min, tm_rbd_max = 50, 70
        tm_full_min, tm_full_max = 68, 82
        gc_min, gc_max = 35, 60

        thead_html = _thead(
            [("Input circRNAs", 5), ("Designed Probes", 10)],
            ["Annotation","Chr","Start","Stop","Strand",
             "TM&nbsp;RBD5","TM&nbsp;RBD3","TM&nbsp;Full",
             "GC%&nbsp;RBD5","GC%&nbsp;RBD3","Ligation&nbsp;Junction",
             "RBD5","BLAST","RBD3","BLAST"],
        )

        rows = []
        for _, row in df.iterrows():
            svg_key = str(row["_svg_key"])
            cells = [
                html_module.escape(str(_safe(row["Annotation"]))),
                html_module.escape(str(_safe(row["Chr"]))),
                html_module.escape(str(_safe(row["Start"]))),
                html_module.escape(str(_safe(row["Stop"]))),
                _strand_cell(row["Strand"]),
                _color_bar(row["TM_left"],  tm_rbd_min, tm_rbd_max),
                _color_bar(row["TM_right"], tm_rbd_min, tm_rbd_max),
                _color_bar(row["TM_Full"],  tm_full_min, tm_full_max),
                _color_bar(row["GC_left"],  gc_min, gc_max),
                _color_bar(row["GC_right"], gc_min, gc_max),
                _junc_badge(row["Ligation_Junction"]),
                _seq_badge(row["RBD5"], int(row["blast_left_n"]),  svg_key),
                _count_badge(int(row["blast_left_n"]),  str(row["BLAST_left_fmt"])),
                _seq_badge(row["RBD3"], int(row["blast_right_n"]), svg_key),
                _count_badge(int(row["blast_right_n"]), str(row["BLAST_right_fmt"])),
            ]
            rows.append("<tr>" + "".join(f"<td>{c}</td>" for c in cells) + "</tr>")

        return _render_page("padlock", exp_name,
                            _table(thead_html, rows), svg_lookup)


    def _format_sirna(self, data_file: str, exp_name: str,
                      output_dir: str, blast_was_run: bool) -> str:
        df = pd.read_csv(data_file, skiprows=1, header=None,
                         keep_default_na=False, dtype=str)

        while df.shape[1] < 14:
            df[df.shape[1]] = ""
        df.columns = list(range(df.shape[1]))

        df[13] = df[13].str.replace("&#10;", "\n", regex=False)\
                        .str.replace("&#9;",  "\t", regex=False)

        df = df.rename(columns={
            1:"Annotation",2:"Chr",3:"Start",4:"Stop",5:"Strand",
            6:"siRNA",7:"newsiRNA",
            8:"Silencing_Score",9:"Rule",10:"Blast",
            11:"Seed_Duplex_Stability",12:"Thermodynamic_Stability",
            13:"SVG",
        })
        for c in ("Silencing_Score","Blast","Seed_Duplex_Stability","Thermodynamic_Stability"):
            df[c] = pd.to_numeric(df[c], errors="coerce")
        df = _clean_location(df)

        svg_lookup, svg_keys = _build_svg_lookup(df)
        df["_svg_key"] = svg_keys

        self._export(df,
                     ["Annotation","Chr","Start","Stop","Strand",
                      "siRNA","Silencing_Score","Rule","Blast",
                      "Seed_Duplex_Stability","Thermodynamic_Stability"],
                     exp_name, output_dir)

        blast_min, blast_max = _symmetric_range(df["Blast"], 0.0)

        thead_html = _thead(
            [("Input circRNAs", 5), ("siRNA Candidates", 7)],
            ["Annotation","Chr","Start","Stop","Strand",
             "siRNA (Guide / Passenger)","Silencing&nbsp;Score","Rule",
             "BLAST","Seed-Duplex&nbsp;Stability","Thermodynamic&nbsp;Stability"],
        )

        rows = []
        for _, row in df.iterrows():
            svg_key = str(row["_svg_key"])
            try:
                blast_n = int(float(row["Blast"])) if not math.isnan(float(row["Blast"])) else 0
            except Exception:
                blast_n = 0
            cls = "hit" if blast_n > 0 else "ok"

            sirna_seq = html_module.escape(str(_safe(row["newsiRNA"])))
            if svg_key:
                magnify_btn = (
                    f'<span class="svg-modal-btn" '
                    f'data-svg-key="{html_module.escape(svg_key)}" '
                    f'data-svg-title="circRNA siRNA position" '
                    f'title="View diagram">&#128269;</span>'
                )
                sirna_cell = (
                    f'<span style="display:inline-flex;align-items:center;gap:4px">'
                    f'<span class="seq-badge {cls}" style="white-space:normal;max-width:340px">'
                    f'{sirna_seq}</span>{magnify_btn}</span>'
                )
            else:
                sirna_cell = (
                    f'<span class="seq-badge {cls}" '
                    f'style="white-space:normal;max-width:340px">'
                    f'{sirna_seq}</span>'
                )

            blast_cell = (
                _color_bar(row["Blast"], blast_min, blast_max)
                if blast_was_run
                else '<span style="color:var(--text-muted)">N/A</span>'
            )

            cells = [
                html_module.escape(str(_safe(row["Annotation"]))),
                html_module.escape(str(_safe(row["Chr"]))),
                html_module.escape(str(_safe(row["Start"]))),
                html_module.escape(str(_safe(row["Stop"]))),
                _strand_cell(row["Strand"]),
                sirna_cell,
                _thresh_badge(row["Silencing_Score"], 50.0, higher_is_bad=False),
                html_module.escape(str(_safe(row["Rule"]))),
                blast_cell,
                _thresh_badge(row["Seed_Duplex_Stability"],  21.5, higher_is_bad=True),
                _thresh_badge(row["Thermodynamic_Stability"], 0.0, higher_is_bad=True),
            ]
            rows.append("<tr>" + "".join(f"<td>{c}</td>" for c in cells) + "</tr>")

        return _render_page("sirna", exp_name,
                            _table(thead_html, rows), svg_lookup)


    def _export(self, df: pd.DataFrame, cols: list[str],
                exp_name: str, output_dir: str):
        if not output_dir:
            return
        os.makedirs(output_dir, exist_ok=True)
        safe = exp_name.replace(" ", "_")
        df[cols].to_csv(  os.path.join(output_dir, f"{safe}_results.csv"),  index=False)
        df[cols].to_excel(os.path.join(output_dir, f"{safe}_results.xlsx"),
                          sheet_name="Results", index=False)
        sys.stderr.write(f"Wrote CSV/XLSX to {output_dir}\n")




def main():
    if len(sys.argv) < 4:
        sys.stderr.write(
            "Usage: circtools_html_formatter.py <mode> <data_file> "
            "<experiment_name> [output_dir] [extra_flag]\n"
        )
        sys.exit(1)

    mode       = sys.argv[1]
    data_file  = sys.argv[2]
    exp_name   = sys.argv[3]
    output_dir = sys.argv[4] if len(sys.argv) > 4 else ""
    extra_flag = sys.argv[5] if len(sys.argv) > 5 else ""

    html = CirctoolsHTMLFormatter(mode).format_file(
        data_file, exp_name, output_dir, extra_flag
    )
    print(html)


if __name__ == "__main__":
    main()