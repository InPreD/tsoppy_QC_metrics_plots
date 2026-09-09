#!/usr/bin/env python3
"""
Plot RNA sample QC metrics from a TSO500 LocalApp metrics table.

Replicates the behaviour of plot_sample_metrics_RNA.R using plotnine.

Usage
-----
    python3 01_plot_sample_metrics_RNA.py \\
        --input_tsv <path> --sample_id <id> --output_pdf <path> \\
        --inpred_nomenclature [True|False]

Arguments
---------
    --input_tsv            Path to the input TSV metrics table.
    --sample_id            Sample ID to extract from the input table (also used as plot title).
    --output_pdf           Path for the output PDF figure.
    --inpred_nomenclature  "True" | "False" — decode InPreD sample-type codes (optional, default: False).

Note: Output PNG is automatically saved alongside the PDF with .png extension.
"""

import argparse
import sys

import numpy as np
import pandas as pd
from plotnine import (
    ggplot,
    aes,
    geom_rect,
    geom_text,
    geom_vline,
    scale_fill_identity,
    scale_color_identity,
    scale_size_identity,
    scale_x_continuous,
    scale_y_continuous,
    labs,
    coord_cartesian,
    theme,
    theme_void,
    element_text,
)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

VERSION = "0.1.0"

COLOR_FAIL  = "#E69F00"
COLOR_PASS  = "#009E73"
COLOR_RED   = "red"
COLOR_BLACK = "black"
COLOR_WHITE = "white"

# InPreD nomenclature: the 13th character (0-based index 12) of a 19-char ID
# encodes the sample type.
INPRED_SAMPLE_TYPES: dict[str, str] = {
    "M": "(metastasis)",
    "T": "(tumor)",
    "N": "(normal)",
    "C": "(cell-line)",
    "L": "(liquid sample)",
    "P": "(primary tumor, naive)",
    "p": "(primary tumor, post-treatment)",
    "R": "(regional metastasis, naive)",
    "r": "(regional metastasis, post-treatment)",
    "D": "(distal metastasis, naive)",
    "d": "(distal metastasis, post-treatment)",
}

# ---------------------------------------------------------------------------
# QC metrics
# ---------------------------------------------------------------------------
# 8 variables: index 0 = boolean (pipeline completed);
# indices 1-7 = continuous RNA metrics mapping one-to-one to PLOT_ROWS[1..7].
#
# RNA_TOTAL_ON_TARGET_READS (index 2) and RNA_TOTAL_PF_READS (index 7) are
# stored as raw counts in the TSV; compute_plot_obs_values divides them by
# 1 000 000 and rounds to 1 decimal, so the thresholds below (9.0 and 30.0)
# are expressed in millions.

# Float metrics that should display with 1 decimal place; all others as integers
FLOAT_METRICS = {"RNA_PCT_CHIMERIC_READS"}

QC_METRICS: list[dict] = [
    # variable                                  lowest    threshold   highest
    {"variable": "COMPLETED_ALL_STEPS",              "lowest": np.nan, "threshold": np.nan, "highest": np.nan},
    {"variable": "RNA_MEDIAN_CV_GENE_500X",          "lowest": 0.0,    "threshold": 93.0,   "highest": 130.0},
    {"variable": "RNA_TOTAL_ON_TARGET_READS",        "lowest": 0.0,    "threshold": 9.0,    "highest": 30.0},
    {"variable": "RNA_MEDIAN_INSERT_SIZE",           "lowest": 0.0,    "threshold": 80.0,   "highest": 200.0},
    {"variable": "RNA_PCT_CHIMERIC_READS",           "lowest": 0.0,    "threshold": 0.0,    "highest": 0.3},
    {"variable": "RNA_PCT_ON_TARGET_READS",          "lowest": 0.0,    "threshold": 0.0,    "highest": 100.0},
    {"variable": "RNA_SCALED_MEDIAN_GENE_COVERAGE",  "lowest": 0.0,    "threshold": 0.0,    "highest": 7000.0},
    {"variable": "RNA_TOTAL_PF_READS",               "lowest": 0.0,    "threshold": 0.0,    "highest": 30.0},
]

# ---------------------------------------------------------------------------
# Plot row definitions (8 rows: 1 boolean + 7 continuous)
# ---------------------------------------------------------------------------
# qc_type values:
#   "TrueFalse" – boolean YES/NO box (green = pass, amber = fail)
#   "FailPass"  – continuous bar; left half = fail (amber), right = pass (green)
#   "PassFail"  – continuous bar; left half = pass (green), right = fail (amber)
#   "NoThold"   – continuous bar drawn entirely in white (no reference threshold)

PLOT_ROWS: list[dict] = [
    {
        "label":           "PIPELINE COMPLETED",
        "threshold_label": "",
        "qc_type":         "TrueFalse",
    },
    {
        "label":           "CV FOR GENES AT MEDIAN COVERAGE >500",
        "threshold_label": "<= 93   [fusions, splice variants]",
        "qc_type":         "PassFail",    # lower CV is better
    },
    {
        "label":           "TOTAL ON-TARGET READS (MILLIONS)",
        "threshold_label": ">= 9   [fusions, splice variants]",
        "qc_type":         "FailPass",
    },
    {
        "label":           "MEDIAN INSERT SIZE",
        "threshold_label": ">= 80   [fusions, splice variants]",
        "qc_type":         "FailPass",
    },
    {
        "label":           "CHIMERIC READS (%)",
        "threshold_label": "No reference value",
        "qc_type":         "NoThold",
    },
    {
        "label":           "ON-TARGET READS (%)",
        "threshold_label": "No reference value",
        "qc_type":         "NoThold",
    },
    {
        "label":           "SCALED MEDIAN GENE COVERAGE",
        "threshold_label": "No reference value",
        "qc_type":         "NoThold",
    },
    {
        "label":           "TOTAL READS PASSING FILTERS (MILLIONS)",
        "threshold_label": "No reference value",
        "qc_type":         "NoThold",
    },
]

N_ROWS = len(PLOT_ROWS)

# ---------------------------------------------------------------------------
# Coordinate layout constants
# ---------------------------------------------------------------------------

# -- metric bar area (x in [0, 2]) --
X_BAR_LEFT  = 0.0
X_BAR_RIGHT = 2.0   # also the normalisation target

# -- boolean (YES/NO) boxes --
X_BOOL_LEFT   = 0.8
X_BOOL_RIGHT  = 1.2
X_BOOL_CENTER = 1.0

# -- separator and conclusion boxes --
X_DIVIDER   = 2.5
X_BOX_LEFT  = 2.8
X_BOX_RIGHT = 3.2

# -- annotation text positions --
X_MIN_LABEL = -0.1    # left of bar: shows the min-range value
X_MAX_LABEL  = 2.2    # right of bar: shows the max-range value
X_ROW_LABEL  = -0.35  # right-aligned row labels (left margin)

# -- overall plot x limits (left sized to centre content, right has 0.3 unit margin) --
X_PLOT_LEFT  = -2.5
X_PLOT_RIGHT =  3.5

# -- y-offset within each row slot (row height = 1 unit; bar height = 0.5) --
BAR_HEIGHT    = 0.5
LABEL_OFFSET  = 0.30   # y above row_bottom for the main (bold) label
THOLD_OFFSET  = 0.0    # y above row_bottom for the threshold description
CENTER_OFFSET = 0.25   # y above row_bottom for centred text inside bars
OBS_OFFSET    = 0.75   # kept for reference; obs labels now use CENTER_OFFSET

# -- observation marker --
MARKER_HALF = 0.01    # half-width of the thin vertical marker rect
MARKER_TOP  = BAR_HEIGHT + 0.1   # extends slightly above the bar top
MARKER_BOT  = -0.1               # extends slightly below the bar bottom

# -- alternating row background colours (zebra-stripe) --
ROW_BG_EVEN = "#F2F2F2"
ROW_BG_ODD  = "#FFFFFF"


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot RNA QC metrics for a single TSO500 sample.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--input_tsv",
        required=True,
        help="Path to the input TSV metrics table.",
    )
    parser.add_argument(
        "--sample_id",
        required=True,
        help="Sample ID to extract from the input table (also used in plot title).",
    )
    parser.add_argument(
        "--output_pdf",
        required=True,
        help="Path for the output PDF figure; PNG is automatically saved alongside with .png extension.",
    )
    parser.add_argument(
        "--inpred_nomenclature",
        choices=["True", "False"],
        default="False",
        help='Decode InPreD sample-type codes from sample_id (optional, default: False).',
    )
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Sample type helper
# ---------------------------------------------------------------------------

def get_sample_type(output_id: str, use_inpred: bool) -> str:
    """Return a sample-type descriptor from InPreD nomenclature when requested."""
    if not use_inpred:
        return "Sample type unspecified"
    if len(output_id) == 19:
        code = output_id[12]
        return INPRED_SAMPLE_TYPES.get(code, "Sample type unknown")
    return "Sample type unspecified"


# ---------------------------------------------------------------------------
# Data loading and value extraction
# ---------------------------------------------------------------------------

def _to_float(val) -> float:
    """
    Convert a raw table value to float.
    Handles 'TRUE'/'FALSE' strings, 'NA' strings, and numeric strings.
    """
    if val is None:
        return np.nan
    if isinstance(val, float) and np.isnan(val):
        return np.nan
    if isinstance(val, str):
        v = val.strip().upper()
        if v in ("NA", "NAN", "N/A", ""):
            return np.nan
        if v == "TRUE":
            return 1.0
        if v == "FALSE":
            return 0.0
        try:
            return float(v)
        except ValueError:
            return np.nan
    try:
        return float(val)
    except (TypeError, ValueError):
        return np.nan


def load_and_extract(tsv_path: str, sample_id: str) -> dict[str, float]:
    """
    Load the metrics TSV and return a dict of {variable: float_value}
    for each QC variable, filtered to *sample_id*.
    """
    df = pd.read_csv(tsv_path, sep="\t", dtype=str)

    mask = df["SAMPLE_ID"] == sample_id
    if not mask.any():
        sys.exit(
            f"ERROR: sample ID '{sample_id}' was not found in '{tsv_path}'.\n"
            f"Available IDs: {df['SAMPLE_ID'].tolist()}"
        )

    row = df.loc[mask].iloc[0]
    return {m["variable"]: _to_float(row.get(m["variable"])) for m in QC_METRICS}


def compute_plot_obs_values(raw: dict[str, float]) -> list[float]:
    """
    Derive the 8 observed values used by the plot.

    Index 0   : COMPLETED_ALL_STEPS  (1.0 = TRUE, 0.0 = FALSE)
    Indices 1-7: Seven continuous RNA metrics (QC_METRICS indices 1-7).
                 Read counts at indices 2 and 7 are divided by 1 000 000
                 and rounded to 1 decimal place before plotting.
    """
    obs: list[float] = [np.nan] * N_ROWS

    obs[0] = raw["COMPLETED_ALL_STEPS"]

    # Continuous metrics map one-to-one: QC_METRICS[j] → PLOT_ROWS[j]
    for j, metric in enumerate(QC_METRICS[1:], start=1):
        obs[j] = raw[metric["variable"]]

    # Scale raw read counts to millions (mirrors R: / 1000000, round(1))
    if not np.isnan(obs[2]):
        obs[2] = round(obs[2] / 1_000_000, 1)
    if not np.isnan(obs[7]):
        obs[7] = round(obs[7] / 1_000_000, 1)

    return obs


# ---------------------------------------------------------------------------
# Low-level coordinate helpers
# ---------------------------------------------------------------------------

def _row_bottom(row_index: int) -> float:
    """Y coordinate of the bottom of a row slot (rows are numbered top-down)."""
    return float(N_ROWS - 1 - row_index)


def _scale_to_bar(value: float, highest: float) -> float:
    """
    Linearly scale *value* from [0, highest] to [0, X_BAR_RIGHT].
    Values above *highest* are clamped to X_BAR_RIGHT (= 2.0).
    """
    if np.isnan(highest) or highest == 0.0:
        return 0.0
    return min(value / highest * X_BAR_RIGHT, X_BAR_RIGHT)


def _format_number(value: float, variable: str = None) -> str:
    """
    Format a float for axis annotation:
    - For FLOAT_METRICS: round to 1 decimal place
    - For others: round to integer (no decimals)
    """
    if variable and variable in FLOAT_METRICS:
        # Float metrics: always show 1 decimal
        rounded = round(value, 1)
        return str(rounded)
    else:
        # All others: show as integer
        return str(int(round(value)))


def _is_value_in_fail_region(obs: float, metric: dict, qc_type: str) -> bool:
    """
    Check if an observed value falls in the fail (orange) region.
    - FailPass: fail if value < threshold
    - PassFail: fail if value > threshold
    """
    if np.isnan(obs):
        return False
    threshold = metric["threshold"]
    if np.isnan(threshold):
        return False
    if qc_type == "FailPass":
        return obs < threshold
    elif qc_type == "PassFail":
        return obs > threshold
    return False


# ---------------------------------------------------------------------------
# Geometry builders
# ---------------------------------------------------------------------------

def build_rect_data(obs_values: list[float]) -> dict[str, pd.DataFrame]:
    """
    Build five DataFrames of rectangle geometry keyed by rendering category:

    "row_bgs" – alternating gray/white full-width row backgrounds
    "boxes"   – empty conclusion boxes (white fill, dark border)
    "bools"   – boolean YES/NO coloured boxes (dark border)
    "bands"   – coloured (or white for NoThold) metric bands (dark border)
    "markers" – thin vertical observation marker (no border)

    All DataFrames share the columns: xmin, xmax, ymin, ymax, fill.
    """
    row_bgs: list[dict] = []
    boxes:   list[dict] = []
    bools:   list[dict] = []
    bands:   list[dict] = []
    markers: list[dict] = []

    # Alternating row backgrounds centred on the bar with equal padding.
    _half_pad = (1.0 - BAR_HEIGHT) / 2
    for i in range(N_ROWS):
        yb = _row_bottom(i)
        row_bgs.append({
            "xmin": X_BAR_LEFT - 0.05,   "xmax": X_BOX_RIGHT + 0.05,
            "ymin": yb - _half_pad,       "ymax": yb + BAR_HEIGHT + _half_pad,
            "fill": ROW_BG_EVEN if i % 2 == 0 else ROW_BG_ODD,
        })

    for i, row in enumerate(PLOT_ROWS):
        yb      = _row_bottom(i)
        yt      = yb + BAR_HEIGHT
        qc_type = row["qc_type"]
        obs     = obs_values[i]

        # Empty conclusion box (always drawn)
        boxes.append({
            "xmin": X_BOX_LEFT, "xmax": X_BOX_RIGHT,
            "ymin": yb,         "ymax": yt,
            "fill": COLOR_WHITE,
        })

        if qc_type == "TrueFalse":
            # Pipeline-completed row: 1 = pass (green), 0 = fail (amber)
            box_fill = COLOR_PASS if obs == 1.0 else COLOR_FAIL
            bools.append({
                "xmin": X_BOOL_LEFT,  "xmax": X_BOOL_RIGHT,
                "ymin": yb,           "ymax": yt,
                "fill": box_fill,
            })

        else:
            # Continuous metric: QC_METRICS[i] maps directly to PLOT_ROWS[i]
            metric  = QC_METRICS[i]
            highest = metric["highest"]
            thr     = metric["threshold"]

            if np.isnan(highest) or highest == 0.0:
                continue

            thr_scaled = _scale_to_bar(thr, highest)
            obs_x      = _scale_to_bar(obs if not np.isnan(obs) else 0.0, highest)

            if qc_type == "FailPass":
                color_left, color_right = COLOR_FAIL, COLOR_PASS
            elif qc_type == "PassFail":
                color_left, color_right = COLOR_PASS, COLOR_FAIL
            else:  # "NoThold": no reference threshold; entire bar is white
                color_left, color_right = COLOR_WHITE, COLOR_WHITE

            # Left region [0, threshold] — omitted when threshold is at the origin
            if thr_scaled > X_BAR_LEFT:
                bands.append({
                    "xmin": X_BAR_LEFT, "xmax": thr_scaled,
                    "ymin": yb,         "ymax": yt,
                    "fill": color_left,
                })
            # Right region [threshold, 2]
            bands.append({
                "xmin": thr_scaled,  "xmax": X_BAR_RIGHT,
                "ymin": yb,          "ymax": yt,
                "fill": color_right,
            })

            # Observation marker: lighter red for visibility against bands and threshold lines.
            # Becomes bright red only when value exceeds the plotted range.
            marker_fill = (
                COLOR_RED if (not np.isnan(obs) and obs > highest) else "#E64C4C"
            )
            markers.append({
                "xmin": obs_x - MARKER_HALF, "xmax": obs_x + MARKER_HALF,
                "ymin": yb + MARKER_BOT,     "ymax": yt + MARKER_TOP - BAR_HEIGHT,
                "fill": marker_fill,
            })

    return {
        "row_bgs": pd.DataFrame(row_bgs),
        "boxes":   pd.DataFrame(boxes),
        "bools":   pd.DataFrame(bools),
        "bands":   pd.DataFrame(bands),
        "markers": pd.DataFrame(markers),
    }


def build_text_data(obs_values: list[float]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Build three DataFrames of text annotations:

    bold_df    – main row labels (bold, right-aligned)
    thold_df   – threshold description labels (regular, right-aligned, smaller)
    center_df  – all other text: YES/NO, min/max range, observed values, and red X for fails
    """
    bold_records:   list[dict] = []
    thold_records:  list[dict] = []
    center_records: list[dict] = []

    for i, row in enumerate(PLOT_ROWS):
        yb  = _row_bottom(i)
        obs = obs_values[i]

        # Main row label (bold, right-aligned)
        bold_records.append({
            "x": X_ROW_LABEL,
            "y": yb + LABEL_OFFSET,
            "label": row["label"],
        })

        # Threshold description (muted, right-aligned, smaller)
        if row["threshold_label"]:
            thold_records.append({
                "x": X_ROW_LABEL,
                "y": yb + THOLD_OFFSET,
                "label": row["threshold_label"],
            })

        qc_type = row["qc_type"]

        if qc_type == "TrueFalse":
            center_records.append({
                "x": X_BOOL_CENTER,
                "y": yb + CENTER_OFFSET,
                "label": "YES" if obs == 1.0 else "NO",
                "color": "#333333",
                "size": 7,
            })

        else:
            # Continuous metric: QC_METRICS[i] maps directly to PLOT_ROWS[i]
            metric  = QC_METRICS[i]
            highest = metric["highest"]
            lowest  = metric["lowest"]

            if np.isnan(highest) or highest == 0.0:
                continue

            obs_x = _scale_to_bar(obs if not np.isnan(obs) else 0.0, highest)

            # Min label (left of bar)
            if not np.isnan(lowest):
                center_records.append({
                    "x": X_MIN_LABEL,
                    "y": yb + CENTER_OFFSET,
                    "label": _format_number(lowest, metric["variable"]),
                    "color": "#333333",
                    "size": 7,
                })

            # Max label (right of bar)
            center_records.append({
                "x": X_MAX_LABEL,
                "y": yb + CENTER_OFFSET,
                "label": _format_number(highest, metric["variable"]),
                "color": "#333333",
                "size": 7,
            })

            # Observed value (centred vertically inside the bar, at its x position)
            if not np.isnan(obs):
                center_records.append({
                    "x": obs_x,
                    "y": yb + CENTER_OFFSET,
                    "label": _format_number(obs, metric["variable"]),
                    "color": "#333333",
                    "size": 7,
                })

            # Red "X" in conclusion box if value is in fail region
            if _is_value_in_fail_region(obs, metric, qc_type):
                center_records.append({
                    "x": (X_BOX_LEFT + X_BOX_RIGHT) / 2.0,
                    "y": yb + CENTER_OFFSET,
                    "label": "X",
                    "color": COLOR_RED,
                    "size": 13,
                })

    return (
        pd.DataFrame(bold_records),
        pd.DataFrame(thold_records),
        pd.DataFrame(center_records),
    )


# ---------------------------------------------------------------------------
# Plot assembly
# ---------------------------------------------------------------------------

def make_plot(obs_values: list[float], title: str) -> ggplot:
    """Assemble and return the plotnine ggplot object."""
    rects                        = build_rect_data(obs_values)
    bold_df, thold_df, center_df = build_text_data(obs_values)

    title_parts = title.split("\n", maxsplit=1)
    title_main  = title_parts[0]
    title_sub   = title_parts[1] if len(title_parts) > 1 else ""

    # Column headers placed just below y = 0 (mirror R's bottom axis labels)
    header_df = pd.DataFrame([
        {"x": X_BOOL_CENTER,                    "y": -0.40, "label": "Observed values"},
        {"x": (X_BOX_LEFT + X_BOX_RIGHT) / 2.0, "y": -0.40, "label": "Conclusions"},
    ])

    base_df = pd.DataFrame({"x": [0.0], "y": [0.0]})

    p = (
        ggplot(base_df, aes("x", "y"))

        # ---- row alternating backgrounds (drawn first, behind everything) --
        + geom_rect(
            data        = rects["row_bgs"],
            mapping     = aes(xmin="xmin", xmax="xmax", ymin="ymin", ymax="ymax",
                              fill="fill"),
            inherit_aes = False,
        )

        # ---- rectangles (back-to-front order) ----------------------------

        # Coloured (or white for NoThold) metric bands
        + geom_rect(
            data        = rects["bands"],
            mapping     = aes(xmin="xmin", xmax="xmax", ymin="ymin", ymax="ymax",
                              fill="fill"),
            color       = "#505050",
            size        = 0.25,
            inherit_aes = False,
        )
        # Observation markers (thin white/red notches, no border)
        + geom_rect(
            data        = rects["markers"],
            mapping     = aes(xmin="xmin", xmax="xmax", ymin="ymin", ymax="ymax",
                              fill="fill"),
            inherit_aes = False,
        )
        # Boolean YES/NO boxes (dark border)
        + geom_rect(
            data        = rects["bools"],
            mapping     = aes(xmin="xmin", xmax="xmax", ymin="ymin", ymax="ymax",
                              fill="fill"),
            color       = "#505050",
            size        = 0.5,
            inherit_aes = False,
        )
        # Empty conclusion boxes (white fill, dark border)
        + geom_rect(
            data        = rects["boxes"],
            mapping     = aes(xmin="xmin", xmax="xmax", ymin="ymin", ymax="ymax",
                              fill="fill"),
            color       = "#505050",
            size        = 0.5,
            inherit_aes = False,
        )
        + scale_fill_identity()

        # ---- vertical separator line ------------------------------------
        + geom_vline(xintercept=X_DIVIDER, color="#666666", size=0.5)

        # ---- text labels (back-to-front hierarchy) ----------------------

        # Main row labels (bold, right-aligned)
        + geom_text(
            data        = bold_df,
            mapping     = aes(x="x", y="y", label="label"),
            ha          = "right",
            va          = "center",
            size        = 8,
            fontweight  = "bold",
            inherit_aes = False,
        )
        # Threshold descriptions (muted, right-aligned, smaller)
        + geom_text(
            data        = thold_df,
            mapping     = aes(x="x", y="y", label="label"),
            ha          = "right",
            va          = "center",
            size        = 6,
            color       = "#777777",
            inherit_aes = False,
        )
        # Centred annotation text (YES/NO, min/max range, observed values, and red X for fails)
        + geom_text(
            data        = center_df,
            mapping     = aes(x="x", y="y", label="label", color="color", size="size"),
            ha          = "center",
            va          = "center",
            inherit_aes = False,
        )
        + scale_color_identity()
        + scale_size_identity()
        # Column headers (bold, centred)
        + geom_text(
            data        = header_df,
            mapping     = aes(x="x", y="y", label="label"),
            ha          = "center",
            va          = "center",
            size        = 7.5,
            fontweight  = "bold",
            color       = "#444444",
            inherit_aes = False,
        )

        # ---- title and subtitle (rendered above panel, never clipped) ---
        + labs(title=title_main, subtitle=title_sub)

        # ---- coordinate system ------------------------------------------
        + coord_cartesian(
            xlim=(X_PLOT_LEFT, X_PLOT_RIGHT),
            ylim=(-0.6, float(N_ROWS)),
        )
        + scale_x_continuous(expand=(0, 0))
        + scale_y_continuous(expand=(0, 0))
        + theme_void()
        + theme(
            plot_title         = element_text(size=10, ha=0.5, fontweight="bold"),
            plot_subtitle      = element_text(size=8,  ha="center", color="#666666"),
            plot_margin_top    = 0.02,
            plot_margin_right  = 0.01,
            plot_margin_bottom = 0.02,
            plot_margin_left   = 0.05,
            figure_size        = (9, 7),
        )
    )
    return p


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    args = parse_args()

    use_inpred  = args.inpred_nomenclature == "True"
    sample_type = get_sample_type(args.sample_id, use_inpred)
    title       = f'QC for RNA sample "{args.sample_id}"\n{sample_type}'

    # Derive output PNG path from output PDF (replace .pdf with .png)
    output_png = args.output_pdf.rsplit(".", 1)[0] + ".png"

    raw_values = load_and_extract(args.input_tsv, args.sample_id)
    obs_values = compute_plot_obs_values(raw_values)

    p = make_plot(obs_values, title)

    p.save(args.output_pdf, width=9, height=7, units="in", verbose=False)
    p.save(output_png, width=9, height=7, units="in", dpi=150, verbose=False)

    print(f"Saved PDF : {args.output_pdf}")
    print(f"Saved PNG : {output_png}")


if __name__ == "__main__":
    main()
