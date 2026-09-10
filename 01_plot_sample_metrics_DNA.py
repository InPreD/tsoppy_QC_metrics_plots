#!/usr/bin/env python3
"""
Plot DNA sample QC metrics from a TSO500 LocalApp metrics table.

Replicates the behaviour of 01_plot_sample_metrics_DNA.R using plotnine.

Usage
-----
    python3 01_plot_sample_metrics_DNA.py \\
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
# The 10 variables extracted from the TSV. Index 0 maps to the "pipeline
# completed" boolean row; indices 1-2 are combined into the "contamination"
# boolean row; indices 3-9 map one-to-one to the seven continuous rows.

# Float metrics that should display with 1 decimal place; all others as integers
FLOAT_METRICS = {"DNA_COVERAGE_MAD"}

QC_METRICS: list[dict] = [
    # variable                         lowest    threshold   highest
    {"variable": "COMPLETED_ALL_STEPS",             "lowest": np.nan, "threshold": np.nan,  "highest": np.nan},
    {"variable": "DNA_CONTAMINATION_SCORE",         "lowest": 0.0,    "threshold": 3106.0,  "highest": np.nan},
    {"variable": "DNA_CONTAMINATION_P_VALUE",       "lowest": 0.0,    "threshold": 0.049,   "highest": 1.0},
    {"variable": "DNA_MEDIAN_EXON_COVERAGE",        "lowest": 0.0,    "threshold": 150.0,   "highest": 1000.0},
    {"variable": "DNA_PCT_EXON_50X",                "lowest": 0.0,    "threshold": 90.0,    "highest": 100.0},
    {"variable": "DNA_MEDIAN_INSERT_SIZE",          "lowest": 0.0,    "threshold": 70.0,    "highest": 200.0},
    {"variable": "DNA_PCT_ALIGNED_READS",           "lowest": 0.0,    "threshold": 70.0,    "highest": 100.0},
    {"variable": "DNA_USABLE_MSI_SITES",            "lowest": 0.0,    "threshold": 40.0,    "highest": 130.0},
    {"variable": "DNA_MEDIAN_BIN_COUNT_CNV_TARGET", "lowest": 0.0,    "threshold": 1.0,     "highest": 20.0},
    {"variable": "DNA_COVERAGE_MAD",                "lowest": 0.0,    "threshold": 0.21,    "highest": 0.3},
]

# ---------------------------------------------------------------------------
# Plot row definitions (9 rows: 2 boolean + 7 continuous)
# ---------------------------------------------------------------------------

PLOT_ROWS: list[dict] = [
    {
        "label":           "PIPELINE COMPLETED",
        "threshold_label": "",
        "qc_type":         "TrueFalse",   # 1 = pass (green), 0 = fail (amber)
    },
    {
        "label":           "CONTAMINATION",
        "threshold_label": "",
        "qc_type":         "TrueFalse",   # 1 = contaminated (amber), 0 = clean (green)
    },
    {
        "label":           "MEDIAN EXON COVERAGE",
        "threshold_label": ">= 150   [small variants, TMB]",
        "qc_type":         "FailPass",
    },
    {
        "label":           "EXON BASES WITH COVERAGE >=50 (%)",
        "threshold_label": ">= 90   [small variants, TMB]",
        "qc_type":         "FailPass",
    },
    {
        "label":           "MEDIAN INSERT SIZE",
        "threshold_label": ">= 70   [small variants, TMB]",
        "qc_type":         "FailPass",
    },
    {
        "label":           "ALIGNED READS (%)",
        "threshold_label": ">= 70",
        "qc_type":         "FailPass",
    },
    {
        "label":           "USABLE MSI SITES",
        "threshold_label": ">= 40   [MSI]",
        "qc_type":         "FailPass",
    },
    {
        "label":           "MEDIAN RAW BIN COUNT PER CNV TARGET",
        "threshold_label": ">= 1   [CNV]",
        "qc_type":         "FailPass",
    },
    {
        "label":           "CNV COVERAGE MAD",
        "threshold_label": "<= 0.21   [CNV]",
        "qc_type":         "PassFail",    # colour semantics reversed
    },
]

N_ROWS = len(PLOT_ROWS)

# ---------------------------------------------------------------------------
# Coordinate layout constants (matching the R base-graphics system)
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
X_MIN_LABEL  = -0.1    # left of bar: shows the min-range value
X_MAX_LABEL   = 2.2    # right of bar: shows the max-range value
X_ROW_LABEL   = -0.35  # right-aligned row labels (left margin)

# -- overall plot x limits (left sized to centre content, right has 0.3 unit margin) --
X_PLOT_LEFT  = -2.5
X_PLOT_RIGHT =  3.5

# -- y-offset within each row slot (row height = 1 unit; bar height = 0.5) --
BAR_HEIGHT    = 0.5
LABEL_OFFSET  = 0.30   # y above row_bottom for the main (bold) label
THOLD_OFFSET  = 0.0    # y above row_bottom for the threshold description
CENTER_OFFSET = 0.25   # y above row_bottom for centred text inside bars
OBS_OFFSET    = 0.75   # y above row_bottom for the observed-value label

# -- observation marker --
MARKER_HALF  = 0.01    # half-width of the thin vertical marker rect
MARKER_TOP   = BAR_HEIGHT + 0.1   # extends slightly above the bar top
MARKER_BOT   = -0.1               # extends slightly below the bar bottom


# -- alternating row background colours (zebra-stripe) --
ROW_BG_EVEN = "#F2F2F2"
ROW_BG_ODD  = "#FFFFFF"


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot DNA QC metrics for a single TSO500 sample.",
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
        code = output_id[12]  # 13th character, 0-based
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
    Derive the 9 observed values used by the plot.

    Index 0   : COMPLETED_ALL_STEPS  (1.0 = TRUE, 0.0 = FALSE)
    Index 1   : Contamination flag   (1.0 = contaminated, 0.0 = clean)
    Indices 2-8: Seven continuous DNA metrics (QC_METRICS indices 3-9)
    """
    obs: list[float] = [np.nan] * N_ROWS

    obs[0] = raw["COMPLETED_ALL_STEPS"]

    # Contamination: both score AND p-value must exceed their respective
    # upper thresholds to call a sample contaminated (replicates R logic).
    score     = raw["DNA_CONTAMINATION_SCORE"]
    pval      = raw["DNA_CONTAMINATION_P_VALUE"]
    score_thr = QC_METRICS[1]["threshold"]
    pval_thr  = QC_METRICS[2]["threshold"]
    contaminated = (
        not np.isnan(score)
        and not np.isnan(pval)
        and score > score_thr
        and pval > pval_thr
    )
    obs[1] = 1.0 if contaminated else 0.0

    # Continuous metrics: QC_METRICS[3..9] → PLOT_ROWS[2..8]
    for j, metric in enumerate(QC_METRICS[3:], start=2):
        obs[j] = raw[metric["variable"]]

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
    Build four DataFrames of rectangle geometry keyed by rendering category:

    "boxes"   – empty conclusion boxes (white fill, black border)
    "bools"   – boolean YES/NO coloured boxes (black border)
    "bands"   – left/right coloured bands for continuous metrics (no border)
    "markers" – thin vertical observation marker (no border)

    All DataFrames share the columns: xmin, xmax, ymin, ymax, fill.
    """
    row_bgs: list[dict] = []
    boxes:   list[dict] = []
    bools:   list[dict] = []
    bands:   list[dict] = []
    markers: list[dict] = []

    # Alternating row backgrounds spanning the full content width.
    # Each stripe is centered on the bar (BAR_HEIGHT tall) with equal padding
    # of (1.0 - BAR_HEIGHT) / 2 = 0.25 above and below.
    _half_pad = (1.0 - BAR_HEIGHT) / 2
    for i in range(N_ROWS):
        yb = _row_bottom(i)
        row_bgs.append({
            "xmin": X_BAR_LEFT - 0.05,    "xmax": X_BOX_RIGHT + 0.05,
            "ymin": yb - _half_pad,        "ymax": yb + BAR_HEIGHT + _half_pad,
            "fill": ROW_BG_EVEN if i % 2 == 0 else ROW_BG_ODD,
        })

    for i, row in enumerate(PLOT_ROWS):
        yb      = _row_bottom(i)
        yt      = yb + BAR_HEIGHT
        qc_type = row["qc_type"]
        obs     = obs_values[i]

        # -- empty conclusion box (always present) -------------------------
        boxes.append({
            "xmin": X_BOX_LEFT, "xmax": X_BOX_RIGHT,
            "ymin": yb,         "ymax": yt,
            "fill": COLOR_WHITE,
        })

        if qc_type == "TrueFalse":
            # ---- boolean coloured box ------------------------------------
            # Row 0 (pipeline complete): 1 = PASS (green), 0 = FAIL (amber)
            # Row 1 (contamination):     1 = FAIL (amber), 0 = PASS (green)
            if i == 0:
                box_fill = COLOR_PASS if obs == 1.0 else COLOR_FAIL
            else:
                box_fill = COLOR_FAIL if obs == 1.0 else COLOR_PASS

            bools.append({
                "xmin": X_BOOL_LEFT,  "xmax": X_BOOL_RIGHT,
                "ymin": yb,           "ymax": yt,
                "fill": box_fill,
            })

        else:
            # ---- continuous metric: coloured band + observation marker ---
            # PLOT_ROWS[i] for i >= 2 maps to QC_METRICS[i + 1]
            metric  = QC_METRICS[i + 1]
            highest = metric["highest"]
            thr     = metric["threshold"]

            if np.isnan(highest) or highest == 0.0:
                continue

            thr_scaled = _scale_to_bar(thr, highest)
            obs_x      = _scale_to_bar(obs if not np.isnan(obs) else 0.0, highest)

            # Colour order depends on whether the threshold is a lower (FailPass)
            # or upper (PassFail) bound.
            if qc_type == "FailPass":
                color_left, color_right = COLOR_FAIL, COLOR_PASS
            else:  # "PassFail"
                color_left, color_right = COLOR_PASS, COLOR_FAIL

            # Left region: x in [0, threshold]
            bands.append({
                "xmin": X_BAR_LEFT,  "xmax": thr_scaled,
                "ymin": yb,          "ymax": yt,
                "fill": color_left,
            })
            # Right region: x in [threshold, 2]
            bands.append({
                "xmin": thr_scaled,  "xmax": X_BAR_RIGHT,
                "ymin": yb,          "ymax": yt,
                "fill": color_right,
            })

            # Observation marker: lighter red for visibility against bands and threshold lines.
            # Becomes bright red only when value exceeds the plotted range.
            marker_fill = (
                COLOR_RED if (not np.isnan(obs) and obs > highest) else "#F67979"
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
    center_df  – all other text: YES/NO, min/max, observed values, headers
    """
    bold_records:   list[dict] = []
    thold_records:  list[dict] = []
    center_records: list[dict] = []

    for i, row in enumerate(PLOT_ROWS):
        yb  = _row_bottom(i)
        obs = obs_values[i]

        # -- main row label (bold, right-aligned) --------------------------
        bold_records.append({
            "x": X_ROW_LABEL,
            "y": yb + LABEL_OFFSET,
            "label": row["label"],
        })

        # -- threshold description (smaller, right-aligned) ----------------
        if row["threshold_label"]:
            thold_records.append({
                "x": X_ROW_LABEL,
                "y": yb + THOLD_OFFSET,
                "label": row["threshold_label"],
            })

        qc_type = row["qc_type"]

        if qc_type == "TrueFalse":
            # YES/NO text centred inside the coloured boolean box
            center_records.append({
                "x": X_BOOL_CENTER,
                "y": yb + CENTER_OFFSET,
                "label": "YES" if obs == 1.0 else "NO",
                "color": "#333333",
                "size": 7,
            })

        else:
            # PLOT_ROWS[i] for i >= 2 maps to QC_METRICS[i + 1]
            metric  = QC_METRICS[i + 1]
            highest = metric["highest"]
            lowest  = metric["lowest"]

            if np.isnan(highest) or highest == 0.0:
                continue

            obs_x = _scale_to_bar(obs if not np.isnan(obs) else 0.0, highest)

            # Min label (left of bar, just below centre)
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
                    "size": 18,
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

    # Split the title string into main heading and subtitle (sample type).
    title_parts = title.split("\n", maxsplit=1)
    title_main  = title_parts[0]
    title_sub   = title_parts[1] if len(title_parts) > 1 else ""

    # Column headers placed just below y = 0 (mirror R's bottom axis labels).
    header_df = pd.DataFrame([
        {"x": X_BOOL_CENTER,                    "y": -0.40, "label": "Observed values"},
        {"x": (X_BOX_LEFT + X_BOX_RIGHT) / 2.0, "y": -0.40, "label": "Conclusions"},
    ])

    # Minimal placeholder required by plotnine for the base ggplot call.
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

        # Coloured band backgrounds (subtle border to frame each bar)
        + geom_rect(
            data        = rects["bands"],
            mapping     = aes(xmin="xmin", xmax="xmax", ymin="ymin", ymax="ymax",
                              fill="fill"),
            color       = "#505050",
            size        = 0.25,
            inherit_aes = False,
        )
        # Observation markers (thin white / red notches, no border)
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

        # Main row labels (bold, right-aligned, larger)
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

        # ---- title and subtitle (rendered above the panel, never clipped) --
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
            # ha=0.5 (float) positions the text at 50% of the panel width,
            # equivalent to ggplot2's hjust=0.5, centering both lines.
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
    title       = f'QC for DNA sample "{args.sample_id}"\n{sample_type}'

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
