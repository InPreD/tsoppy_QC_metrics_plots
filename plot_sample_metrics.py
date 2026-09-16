#!/usr/bin/env python3
"""
Plot per-sample QC metrics from a TSO500 LocalApp metrics table.

Supports DNA and RNA sample types via --seq_type.  The seq-type determines
which metrics, thresholds, and QC logic are applied; the rendering pipeline
is shared and generic.

Usage
-----
    python3 01_plot_sample_metrics.py \\
        --input_tsv  <path>         \\
        --sample_id  <id>           \\
        --output_pdf <path>         \\
        --seq_type   [DNA|RNA]      \\
        --inpred_nomenclature [True|False]

Arguments
---------
    --input_tsv            Path to the input TSV metrics table.
    --sample_id            Sample ID to extract (also used in the plot title).
    --output_pdf           Path for the output PDF figure.
    --seq_type             "DNA" or "RNA" — selects metrics and QC thresholds.
    --inpred_nomenclature  "True" | "False" — decode InPreD sample-type codes
                           from sample_id (optional, default: False).

Note: A PNG is automatically saved alongside the PDF (same path, .png extension).
"""

import argparse
import sys

import numpy as np
import polars as pl
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
# Shared visual constants
# ---------------------------------------------------------------------------

VERSION = "1.0.0"

COLOR_FAIL  = "#E69F00"
COLOR_PASS  = "#009E73"
COLOR_RED   = "red"
COLOR_WHITE = "white"

# InPreD nomenclature: 13th character (0-based index 12) of a 19-char ID.
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

# -- metric bar area (x in [0, 2]) --
X_BAR_LEFT  = 0.0
X_BAR_RIGHT = 2.0   # normalisation target

# -- boolean YES/NO box --
X_BOOL_LEFT   = 0.8
X_BOOL_RIGHT  = 1.2
X_BOOL_CENTER = 1.0

# -- vertical separator and conclusion boxes --
X_DIVIDER   = 2.5
X_BOX_LEFT  = 2.8
X_BOX_RIGHT = 3.2

# -- annotation text x positions --
X_MIN_LABEL = -0.1    # range minimum, left of bar
X_MAX_LABEL  = 2.2    # range maximum, right of bar
X_ROW_LABEL  = -0.35  # row label, right-aligned in left margin

# -- overall plot x limits --
X_PLOT_LEFT  = -2.5
X_PLOT_RIGHT =  3.5

# -- y offsets within each 1-unit row slot (bar height = 0.5) --
BAR_HEIGHT    = 0.5
LABEL_OFFSET  = 0.30   # main (bold) row label
THOLD_OFFSET  = 0.0    # threshold description
CENTER_OFFSET = 0.25   # centred text inside bars / conclusion boxes

OBS_LABEL_OFFSET = 0.09
X_BAR_MID        = (X_BAR_LEFT + X_BAR_RIGHT) / 2.0   # = 1.0

# -- observation marker --
MARKER_HALF = 0.01
MARKER_TOP  = BAR_HEIGHT + 0.1   # extends slightly above bar top
MARKER_BOT  = -0.1               # extends slightly below bar bottom

# -- alternating row backgrounds --
ROW_BG_EVEN = "#F2F2F2"
ROW_BG_ODD  = "#FFFFFF"


# ---------------------------------------------------------------------------
# DNA configuration
# ---------------------------------------------------------------------------
# QC_METRICS layout:
#   Index 0:   pipeline completion  → boolean row 0
#   Index 1-2: contamination score + p-value → boolean row 1 (derived)
#   Index 3-9: seven continuous metrics → plot rows 2-8
#
# Each PLOT_ROW carries a direct reference to its QC_METRICS entry ("metric"),
# eliminating index-offset arithmetic in the rendering functions.
# Boolean rows store the value that constitutes PASS in "pass_value".

DNA_FLOAT_METRICS: set[str] = {"DNA_COVERAGE_MAD"}

DNA_QC_METRICS: list[dict] = [
    #                                               lowest    threshold   highest
    {"variable": "COMPLETED_ALL_STEPS",             "lowest": np.nan, "threshold": np.nan, "highest": np.nan},
    {"variable": "DNA_CONTAMINATION_SCORE",         "lowest": 0.0,   "threshold": 3106.0, "highest": np.nan},
    {"variable": "DNA_CONTAMINATION_P_VALUE",       "lowest": 0.0,   "threshold": 0.049,  "highest": 1.0},
    {"variable": "DNA_MEDIAN_EXON_COVERAGE",        "lowest": 0.0,   "threshold": 150.0,  "highest": 1000.0},
    {"variable": "DNA_PCT_EXON_50X",                "lowest": 0.0,   "threshold": 90.0,   "highest": 100.0},
    {"variable": "DNA_MEDIAN_INSERT_SIZE",          "lowest": 0.0,   "threshold": 70.0,   "highest": 200.0},
    {"variable": "DNA_PCT_ALIGNED_READS",           "lowest": 0.0,   "threshold": 70.0,   "highest": 100.0},
    {"variable": "DNA_USABLE_MSI_SITES",            "lowest": 0.0,   "threshold": 40.0,   "highest": 130.0},
    {"variable": "DNA_MEDIAN_BIN_COUNT_CNV_TARGET", "lowest": 0.0,   "threshold": 1.0,    "highest": 20.0},
    {"variable": "DNA_COVERAGE_MAD",                "lowest": 0.0,   "threshold": 0.21,   "highest": 0.3},
]

DNA_PLOT_ROWS: list[dict] = [
    {
        "label":           "PIPELINE COMPLETED",
        "threshold_label": "",
        "qc_type":         "TrueFalse",
        "pass_value":      1.0,    # obs == 1 → pipeline ran → PASS (green)
        "metric":          None,
    },
    {
        "label":           "CONTAMINATION",
        "threshold_label": "",
        "qc_type":         "TrueFalse",
        "pass_value":      0.0,    # obs == 0 → sample clean → PASS (green)
        "metric":          None,
    },
    {
        "label":           "MEDIAN EXON COVERAGE",
        "threshold_label": ">= 150   [small variants, TMB]",
        "qc_type":         "FailPass",
        "pass_value":      None,
        "metric":          DNA_QC_METRICS[3],
    },
    {
        "label":           "EXON BASES WITH COVERAGE >=50 (%)",
        "threshold_label": ">= 90   [small variants, TMB]",
        "qc_type":         "FailPass",
        "pass_value":      None,
        "metric":          DNA_QC_METRICS[4],
    },
    {
        "label":           "MEDIAN INSERT SIZE",
        "threshold_label": ">= 70   [small variants, TMB]",
        "qc_type":         "FailPass",
        "pass_value":      None,
        "metric":          DNA_QC_METRICS[5],
    },
    {
        "label":           "ALIGNED READS (%)",
        "threshold_label": ">= 70",
        "qc_type":         "FailPass",
        "pass_value":      None,
        "metric":          DNA_QC_METRICS[6],
    },
    {
        "label":           "USABLE MSI SITES",
        "threshold_label": ">= 40   [MSI]",
        "qc_type":         "FailPass",
        "pass_value":      None,
        "metric":          DNA_QC_METRICS[7],
    },
    {
        "label":           "MEDIAN RAW BIN COUNT PER CNV TARGET",
        "threshold_label": ">= 1   [CNV]",
        "qc_type":         "FailPass",
        "pass_value":      None,
        "metric":          DNA_QC_METRICS[8],
    },
    {
        "label":           "CNV COVERAGE MAD",
        "threshold_label": "<= 0.21   [CNV]",
        "qc_type":         "PassFail",   # lower is better; right half = fail
        "pass_value":      None,
        "metric":          DNA_QC_METRICS[9],
    },
]


# ---------------------------------------------------------------------------
# RNA configuration
# ---------------------------------------------------------------------------
# QC_METRICS layout:
#   Index 0:   pipeline completion  → boolean row 0
#   Index 1-7: continuous metrics   → plot rows 1-7 (direct 1:1 mapping)
#
# RNA_TOTAL_ON_TARGET_READS (index 2) and RNA_TOTAL_PF_READS (index 7) are
# raw counts in the TSV; compute_plot_obs_values scales them to millions.
# Thresholds are therefore expressed in millions.
#
# qc_type "NoThold": metric is informational only; no pass/fail threshold.

RNA_FLOAT_METRICS: set[str] = {"RNA_PCT_CHIMERIC_READS"}

RNA_QC_METRICS: list[dict] = [
    #                                                   lowest    threshold  highest
    {"variable": "COMPLETED_ALL_STEPS",                 "lowest": np.nan, "threshold": np.nan, "highest": np.nan},
    {"variable": "RNA_MEDIAN_CV_GENE_500X",             "lowest": 0.0,   "threshold": 93.0,   "highest": 130.0},
    {"variable": "RNA_TOTAL_ON_TARGET_READS",           "lowest": 0.0,   "threshold": 9.0,    "highest": 30.0},
    {"variable": "RNA_MEDIAN_INSERT_SIZE",              "lowest": 0.0,   "threshold": 80.0,   "highest": 200.0},
    {"variable": "RNA_PCT_CHIMERIC_READS",              "lowest": 0.0,   "threshold": 0.0,    "highest": 0.3},
    {"variable": "RNA_PCT_ON_TARGET_READS",             "lowest": 0.0,   "threshold": 0.0,    "highest": 100.0},
    {"variable": "RNA_SCALED_MEDIAN_GENE_COVERAGE",     "lowest": 0.0,   "threshold": 0.0,    "highest": 7000.0},
    {"variable": "RNA_TOTAL_PF_READS",                  "lowest": 0.0,   "threshold": 0.0,    "highest": 30.0},
]

RNA_PLOT_ROWS: list[dict] = [
    {
        "label":           "PIPELINE COMPLETED",
        "threshold_label": "",
        "qc_type":         "TrueFalse",
        "pass_value":      1.0,
        "metric":          None,
    },
    {
        "label":           "CV FOR GENES AT MEDIAN COVERAGE >500",
        "threshold_label": "<= 93   [fusions, splice variants]",
        "qc_type":         "PassFail",   # lower CV is better; right half = fail
        "pass_value":      None,
        "metric":          RNA_QC_METRICS[1],
    },
    {
        "label":           "TOTAL ON-TARGET READS (MILLIONS)",
        "threshold_label": ">= 9   [fusions, splice variants]",
        "qc_type":         "FailPass",
        "pass_value":      None,
        "metric":          RNA_QC_METRICS[2],
    },
    {
        "label":           "MEDIAN INSERT SIZE",
        "threshold_label": ">= 80   [fusions, splice variants]",
        "qc_type":         "FailPass",
        "pass_value":      None,
        "metric":          RNA_QC_METRICS[3],
    },
    {
        "label":           "CHIMERIC READS (%)",
        "threshold_label": "No reference value",
        "qc_type":         "NoThold",
        "pass_value":      None,
        "metric":          RNA_QC_METRICS[4],
    },
    {
        "label":           "ON-TARGET READS (%)",
        "threshold_label": "No reference value",
        "qc_type":         "NoThold",
        "pass_value":      None,
        "metric":          RNA_QC_METRICS[5],
    },
    {
        "label":           "SCALED MEDIAN GENE COVERAGE",
        "threshold_label": "No reference value",
        "qc_type":         "NoThold",
        "pass_value":      None,
        "metric":          RNA_QC_METRICS[6],
    },
    {
        "label":           "TOTAL READS PASSING FILTERS (MILLIONS)",
        "threshold_label": "No reference value",
        "qc_type":         "NoThold",
        "pass_value":      None,
        "metric":          RNA_QC_METRICS[7],
    },
]


# ---------------------------------------------------------------------------
# Configuration registry
# ---------------------------------------------------------------------------

SEQ_TYPE_CONFIGS: dict[str, dict] = {
    "DNA": {
        "qc_metrics":    DNA_QC_METRICS,
        "plot_rows":     DNA_PLOT_ROWS,
        "float_metrics": DNA_FLOAT_METRICS,
        "seq_label":     "DNA",
    },
    "RNA": {
        "qc_metrics":    RNA_QC_METRICS,
        "plot_rows":     RNA_PLOT_ROWS,
        "float_metrics": RNA_FLOAT_METRICS,
        "seq_label":     "RNA",
    },
}


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot QC metrics for a single TSO500 DNA or RNA sample.",
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
        help="Sample ID to extract from the table (also used in the plot title).",
    )
    parser.add_argument(
        "--output_pdf",
        required=True,
        help="Output PDF path; a PNG is automatically saved alongside with .png extension.",
    )
    parser.add_argument(
        "--seq_type",
        required=True,
        choices=["DNA", "RNA"],
        help='"DNA" or "RNA" — selects which metrics and thresholds to apply.',
    )
    parser.add_argument(
        "--inpred_nomenclature",
        choices=["True", "False"],
        default="False",
        help="Decode InPreD sample-type codes from sample_id (default: False).",
    )
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Sample type helper
# ---------------------------------------------------------------------------

def get_sample_type(output_id: str, use_inpred: bool) -> str:
    """Return a sample-type descriptor decoded from InPreD nomenclature when requested."""
    if not use_inpred:
        return "Sample type unspecified"
    if len(output_id) == 19:
        code = output_id[12]   # 13th character, 0-based
        return INPRED_SAMPLE_TYPES.get(code, "Sample type unknown")
    return "Sample type unspecified"


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def _to_float(val) -> float:
    """Convert a raw table cell to float; handles TRUE/FALSE strings, NA, and numerics."""
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


def load_and_extract(
    tsv_path: str, sample_id: str, qc_metrics: list[dict]
) -> dict[str, float]:
    """Load the metrics TSV and return {variable: float_value} for the given sample."""
    df       = pl.read_csv(tsv_path, separator="\t", infer_schema=False)
    filtered = df.filter(pl.col("SAMPLE_ID") == sample_id)
    if filtered.is_empty():
        sys.exit(
            f"ERROR: sample ID '{sample_id}' not found in '{tsv_path}'.\n"
            f"Available IDs: {df['SAMPLE_ID'].to_list()}"
        )
    row = filtered.row(0, named=True)
    return {m["variable"]: _to_float(row.get(m["variable"])) for m in qc_metrics}


# ---------------------------------------------------------------------------
# Observed values computation (seq-type-specific)
# ---------------------------------------------------------------------------

def _compute_obs_dna(raw: dict[str, float], qc_metrics: list[dict]) -> list[float]:
    """
    Derive 9 observed values for the DNA plot.

    Row 0   : COMPLETED_ALL_STEPS flag (1 = TRUE, 0 = FALSE)
    Row 1   : contamination flag, derived from score + p-value (1 = contaminated)
    Rows 2-8: continuous metrics at qc_metrics indices 3-9
    """
    obs: list[float] = [np.nan] * 9
    obs[0] = raw["COMPLETED_ALL_STEPS"]

    score, pval = raw["DNA_CONTAMINATION_SCORE"], raw["DNA_CONTAMINATION_P_VALUE"]
    obs[1] = 1.0 if (
        not np.isnan(score) and not np.isnan(pval)
        and score > qc_metrics[1]["threshold"]
        and pval  > qc_metrics[2]["threshold"]
    ) else 0.0

    for j, metric in enumerate(qc_metrics[3:], start=2):
        obs[j] = raw[metric["variable"]]

    return obs


def _compute_obs_rna(raw: dict[str, float], qc_metrics: list[dict]) -> list[float]:
    """
    Derive 8 observed values for the RNA plot.

    Row 0   : COMPLETED_ALL_STEPS flag
    Rows 1-7: continuous metrics at qc_metrics indices 1-7;
              read counts at rows 2 and 7 are scaled from raw counts to millions.
    """
    obs: list[float] = [np.nan] * 8
    obs[0] = raw["COMPLETED_ALL_STEPS"]

    for j, metric in enumerate(qc_metrics[1:], start=1):
        obs[j] = raw[metric["variable"]]

    # Scale raw read counts to millions (thresholds defined in millions above)
    obs[2] = round(obs[2] / 1_000_000, 1) if not np.isnan(obs[2]) else np.nan
    obs[7] = round(obs[7] / 1_000_000, 1) if not np.isnan(obs[7]) else np.nan

    return obs


def compute_plot_obs_values(
    raw: dict[str, float], seq_type: str, qc_metrics: list[dict]
) -> list[float]:
    """Dispatch to the seq-type-specific observed-values computation."""
    if seq_type == "DNA":
        return _compute_obs_dna(raw, qc_metrics)
    return _compute_obs_rna(raw, qc_metrics)


# ---------------------------------------------------------------------------
# Low-level coordinate helpers
# ---------------------------------------------------------------------------

def _row_bottom(row_index: int, n_rows: int) -> float:
    """Y coordinate of the bottom of a row slot (rows are ordered top-down)."""
    return float(n_rows - 1 - row_index)


def _scale_to_bar(value: float, highest: float) -> float:
    """Linearly scale *value* from [0, highest] to [0, X_BAR_RIGHT]; clamps at X_BAR_RIGHT."""
    if np.isnan(highest) or highest == 0.0:
        return 0.0
    return min(value / highest * X_BAR_RIGHT, X_BAR_RIGHT)


def _format_number(value: float, variable: str, float_metrics: set[str]) -> str:
    """
    Format a float for axis annotation.
    Metrics in float_metrics display 1 decimal place; all others are rounded to integer.
    """
    if variable in float_metrics:
        return str(round(value, 1))
    return str(int(round(value)))


def _is_value_in_fail_region(obs: float, metric: dict, qc_type: str) -> bool:
    """Return True when obs falls in the fail (amber) colour band."""
    if np.isnan(obs):
        return False
    threshold = metric["threshold"]
    if np.isnan(threshold):
        return False
    if qc_type == "FailPass":
        return obs < threshold
    if qc_type == "PassFail":
        return obs > threshold
    return False   # NoThold and TrueFalse rows never trigger a red X


# ---------------------------------------------------------------------------
# Per-row-type rendering helpers
# ---------------------------------------------------------------------------

def _build_boolean_row(
    i: int, row: dict, obs: float, n_rows: int, float_metrics: set[str]
) -> dict[str, list[dict]]:
    """
    Build geometry and text records for a TrueFalse (boolean) row.

    The coloured box spans [X_BOOL_LEFT, X_BOOL_RIGHT].  pass_value in the
    row definition declares which observed value constitutes a pass (green);
    any other value is rendered as a fail (amber).
    """
    yb       = _row_bottom(i, n_rows)
    yt       = yb + BAR_HEIGHT
    box_fill = COLOR_PASS if obs == row["pass_value"] else COLOR_FAIL
    return {
        "bools":     [{"xmin": X_BOOL_LEFT, "xmax": X_BOOL_RIGHT, "ymin": yb, "ymax": yt, "fill": box_fill}],
        "bands":     [],
        "markers":   [],
        "center":    [{"x": X_BOOL_CENTER, "y": yb + CENTER_OFFSET,
                        "label": "YES" if obs == 1.0 else "NO", "color": "#333333", "size": 7}],
        "obs_left":  [],
        "obs_right": [],
    }


def _build_continuous_row(
    i: int, row: dict, obs: float, n_rows: int, float_metrics: set[str]
) -> dict[str, list[dict]] | None:
    """
    Build geometry and text records for a continuous metric row.

    Handles FailPass, PassFail, and NoThold qc_types.
    Returns None when the row has no valid display range (highest is nan or 0),
    so the caller can skip it cleanly.
    """
    metric  = row["metric"]
    highest = metric["highest"]
    lowest  = metric["lowest"]
    thr     = metric["threshold"]
    qc_type = row["qc_type"]

    if np.isnan(highest) or highest == 0.0:
        return None

    yb         = _row_bottom(i, n_rows)
    yt         = yb + BAR_HEIGHT
    thr_scaled = _scale_to_bar(thr, highest)
    obs_x      = _scale_to_bar(obs if not np.isnan(obs) else 0.0, highest)

    if qc_type == "FailPass":
        color_left, color_right = COLOR_FAIL, COLOR_PASS
    elif qc_type == "PassFail":
        color_left, color_right = COLOR_PASS, COLOR_FAIL
    else:  # NoThold: informational only, rendered white throughout
        color_left, color_right = COLOR_WHITE, COLOR_WHITE

    # Geometry ---------------------------------------------------------------
    bands: list[dict] = []
    # Left region [0, threshold] — skipped when threshold is at the origin
    if thr_scaled > X_BAR_LEFT:
        bands.append({"xmin": X_BAR_LEFT, "xmax": thr_scaled, "ymin": yb, "ymax": yt, "fill": color_left})
    # Right region [threshold, bar_right]
    bands.append({"xmin": thr_scaled, "xmax": X_BAR_RIGHT, "ymin": yb, "ymax": yt, "fill": color_right})

    # Softer red marker; bright red only when value exceeds the plotted range
    marker_fill = COLOR_RED if (not np.isnan(obs) and obs > highest) else "#D4534A"
    markers = [{"xmin": obs_x - MARKER_HALF, "xmax": obs_x + MARKER_HALF,
                 "ymin": yb + MARKER_BOT,     "ymax": yt + MARKER_TOP - BAR_HEIGHT,
                 "fill": marker_fill}]

    # Text annotations -------------------------------------------------------
    center: list[dict] = []
    if not np.isnan(lowest):
        center.append({"x": X_MIN_LABEL, "y": yb + CENTER_OFFSET,
                        "label": _format_number(lowest, metric["variable"], float_metrics),
                        "color": "#444444", "size": 8.5})
    center.append({"x": X_MAX_LABEL, "y": yb + CENTER_OFFSET,
                   "label": _format_number(highest, metric["variable"], float_metrics),
                   "color": "#444444", "size": 8.5})
    # Red X in the conclusion box for metrics that fall in the fail region
    if _is_value_in_fail_region(obs, metric, qc_type):
        center.append({"x": (X_BOX_LEFT + X_BOX_RIGHT) / 2.0, "y": yb + CENTER_OFFSET,
                        "label": "X", "color": COLOR_RED, "size": 18})

    # Observed value: placed left or right of the marker to avoid overlap
    obs_left:  list[dict] = []
    obs_right: list[dict] = []
    if not np.isnan(obs):
        fmt = _format_number(obs, metric["variable"], float_metrics)
        if obs_x > X_BAR_MID:
            obs_left.append({"x": obs_x - OBS_LABEL_OFFSET, "y": yb + CENTER_OFFSET,
                              "label": fmt, "color": "#444444", "size": 8.5})
        else:
            obs_right.append({"x": obs_x + OBS_LABEL_OFFSET, "y": yb + CENTER_OFFSET,
                               "label": fmt, "color": "#444444", "size": 8.5})

    return {
        "bools":     [],
        "bands":     bands,
        "markers":   markers,
        "center":    center,
        "obs_left":  obs_left,
        "obs_right": obs_right,
    }


# Dispatch table: maps each qc_type to its row-builder function.
# To support a new qc_type, implement a builder with the same signature
# (i, row, obs, n_rows, float_metrics) and register it here.
_ROW_BUILDERS: dict[str, callable] = {
    "TrueFalse": _build_boolean_row,
    "FailPass":  _build_continuous_row,
    "PassFail":  _build_continuous_row,
    "NoThold":   _build_continuous_row,
}


# ---------------------------------------------------------------------------
# Unified layer builder (single pass over plot_rows)
# ---------------------------------------------------------------------------

def build_layer_data(
    obs_values:    list[float],
    plot_rows:     list[dict],
    float_metrics: set[str],
) -> dict[str, pl.DataFrame]:
    """
    Build all geometry and annotation DataFrames in a single pass over plot_rows.

    Row-type-specific rendering is delegated to the functions registered in
    _ROW_BUILDERS.  Adding or removing a metric only requires updating the
    plot_rows configuration list; adding a new qc_type only requires registering
    a new builder in _ROW_BUILDERS.

    Returns a dict keyed by layer name.

    Geometry layers (columns: xmin, xmax, ymin, ymax, fill):
      "row_bgs"   – alternating zebra-stripe backgrounds (full row width)
      "boxes"     – white conclusion boxes with dark border
      "bools"     – coloured boolean YES/NO boxes with dark border
      "bands"     – coloured (or white for NoThold) metric bands with thin border
      "markers"   – thin vertical observation-position notches, no border

    Text layers:
      "bold"      – bold row labels (columns: x, y, label)
      "thold"     – threshold descriptions (columns: x, y, label)
      "center"    – centred annotations: YES/NO, range limits, red X  (columns: x, y, label, color, size)
      "obs_left"  – observed values left of the marker  (columns: x, y, label, color, size)
      "obs_right" – observed values right of the marker (columns: x, y, label, color, size)
    """
    n_rows   = len(plot_rows)
    half_pad = (1.0 - BAR_HEIGHT) / 2

    acc: dict[str, list[dict]] = {
        "row_bgs": [], "boxes": [], "bools": [], "bands": [], "markers": [],
        "bold": [], "thold": [], "center": [], "obs_left": [], "obs_right": [],
    }

    for i, row in enumerate(plot_rows):
        yb  = _row_bottom(i, n_rows)
        obs = obs_values[i]

        # Structural elements common to every row
        acc["row_bgs"].append({
            "xmin": X_BAR_LEFT - 0.05, "xmax": X_BOX_RIGHT + 0.05,
            "ymin": yb - half_pad,     "ymax": yb + BAR_HEIGHT + half_pad,
            "fill": ROW_BG_EVEN if i % 2 == 0 else ROW_BG_ODD,
        })
        acc["boxes"].append({
            "xmin": X_BOX_LEFT, "xmax": X_BOX_RIGHT,
            "ymin": yb,         "ymax": yb + BAR_HEIGHT,
            "fill": COLOR_WHITE,
        })

        # Labels common to every row
        acc["bold"].append({"x": X_ROW_LABEL, "y": yb + LABEL_OFFSET, "label": row["label"]})
        if row["threshold_label"]:
            acc["thold"].append({"x": X_ROW_LABEL, "y": yb + THOLD_OFFSET, "label": row["threshold_label"]})

        # Dispatch row-type-specific geometry and text
        builder = _ROW_BUILDERS[row["qc_type"]]
        result  = builder(i, row, obs, n_rows, float_metrics)
        if result is not None:
            for key in ("bools", "bands", "markers", "center", "obs_left", "obs_right"):
                acc[key].extend(result[key])

    rect_cols  = ["xmin", "xmax", "ymin", "ymax", "fill"]
    text_cols  = ["x", "y", "label"]
    annot_cols = ["x", "y", "label", "color", "size"]

    def _df(records: list[dict], cols: list[str]) -> pl.DataFrame:
        return pl.DataFrame(records) if records else pl.DataFrame({c: [] for c in cols})

    return {
        "row_bgs":  _df(acc["row_bgs"],   rect_cols),
        "boxes":    _df(acc["boxes"],     rect_cols),
        "bools":    _df(acc["bools"],     rect_cols),
        "bands":    _df(acc["bands"],     rect_cols),
        "markers":  _df(acc["markers"],   rect_cols),
        "bold":     _df(acc["bold"],      text_cols),
        "thold":    _df(acc["thold"],     text_cols),
        "center":   _df(acc["center"],    annot_cols),
        "obs_left": _df(acc["obs_left"],  annot_cols),
        "obs_right":_df(acc["obs_right"], annot_cols),
    }


# ---------------------------------------------------------------------------
# Plot assembly
# ---------------------------------------------------------------------------

def make_plot(
    obs_values:    list[float],
    title:         str,
    plot_rows:     list[dict],
    float_metrics: set[str],
) -> ggplot:
    """Assemble and return the complete plotnine ggplot object."""
    n_rows    = len(plot_rows)
    layers    = build_layer_data(obs_values, plot_rows, float_metrics)
    pd_layers = {k: v.to_pandas() for k, v in layers.items()}

    title_parts = title.split("\n", maxsplit=1)
    title_main  = title_parts[0]
    title_sub   = title_parts[1] if len(title_parts) > 1 else ""

    header_df = pl.DataFrame([
        {"x": X_BOOL_CENTER,                     "y": -0.40, "label": "Observed values"},
        {"x": (X_BOX_LEFT + X_BOX_RIGHT) / 2.0,  "y": -0.40, "label": "Conclusions"},
    ]).to_pandas()

    base_df = pl.DataFrame({"x": [0.0], "y": [0.0]}).to_pandas()

    p = (
        ggplot(base_df, aes("x", "y"))

        # Row alternating backgrounds (drawn first, behind everything)
        + geom_rect(
            data        = pd_layers["row_bgs"],
            mapping     = aes(xmin="xmin", xmax="xmax", ymin="ymin", ymax="ymax", fill="fill"),
            inherit_aes = False,
        )
        # Coloured metric bands (thin border frames each bar)
        + geom_rect(
            data        = pd_layers["bands"],
            mapping     = aes(xmin="xmin", xmax="xmax", ymin="ymin", ymax="ymax", fill="fill"),
            color       = "#505050",
            size        = 0.25,
            inherit_aes = False,
        )
        # Observation markers (thin vertical notches, no border)
        + geom_rect(
            data        = pd_layers["markers"],
            mapping     = aes(xmin="xmin", xmax="xmax", ymin="ymin", ymax="ymax", fill="fill"),
            inherit_aes = False,
        )
        # Boolean YES/NO boxes (dark border)
        + geom_rect(
            data        = pd_layers["bools"],
            mapping     = aes(xmin="xmin", xmax="xmax", ymin="ymin", ymax="ymax", fill="fill"),
            color       = "#505050",
            size        = 0.5,
            inherit_aes = False,
        )
        # Conclusion boxes (white fill, dark border)
        + geom_rect(
            data        = pd_layers["boxes"],
            mapping     = aes(xmin="xmin", xmax="xmax", ymin="ymin", ymax="ymax", fill="fill"),
            color       = "#505050",
            size        = 0.5,
            inherit_aes = False,
        )
        + scale_fill_identity()

        # Vertical separator between the metric bar area and the conclusions column
        + geom_vline(xintercept=X_DIVIDER, color="#666666", size=0.5)

        # Text labels (back-to-front rendering order)
        + geom_text(
            data        = pd_layers["bold"],
            mapping     = aes(x="x", y="y", label="label"),
            ha          = "right",
            va          = "center",
            size        = 8,
            fontweight  = "bold",
            inherit_aes = False,
        )
        + geom_text(
            data        = pd_layers["thold"],
            mapping     = aes(x="x", y="y", label="label"),
            ha          = "right",
            va          = "center",
            size        = 6,
            color       = "#777777",
            inherit_aes = False,
        )
        + geom_text(
            data        = pd_layers["center"],
            mapping     = aes(x="x", y="y", label="label", color="color", size="size"),
            ha          = "center",
            va          = "center",
            inherit_aes = False,
        )
        + geom_text(
            data        = pd_layers["obs_left"],
            mapping     = aes(x="x", y="y", label="label", color="color", size="size"),
            ha          = "right",
            va          = "center",
            inherit_aes = False,
        )
        + geom_text(
            data        = pd_layers["obs_right"],
            mapping     = aes(x="x", y="y", label="label", color="color", size="size"),
            ha          = "left",
            va          = "center",
            inherit_aes = False,
        )
        + scale_color_identity()
        + scale_size_identity()
        # Column headers (drawn last so they sit on top)
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

        + labs(title=title_main, subtitle=title_sub)
        + coord_cartesian(
            xlim=(X_PLOT_LEFT, X_PLOT_RIGHT),
            ylim=(-0.6, float(n_rows)),
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
    args   = parse_args()
    config = SEQ_TYPE_CONFIGS[args.seq_type]

    use_inpred  = args.inpred_nomenclature == "True"
    sample_type = get_sample_type(args.sample_id, use_inpred)
    title       = f'QC for {config["seq_label"]} sample "{args.sample_id}"\n{sample_type}'
    output_png  = args.output_pdf.rsplit(".", 1)[0] + ".png"

    raw_values = load_and_extract(args.input_tsv, args.sample_id, config["qc_metrics"])
    obs_values = compute_plot_obs_values(raw_values, args.seq_type, config["qc_metrics"])

    p = make_plot(obs_values, title, config["plot_rows"], config["float_metrics"])

    p.save(args.output_pdf, width=9, height=7, units="in", verbose=False)
    p.save(output_png,      width=9, height=7, units="in", dpi=150, verbose=False)

    print(f"Saved PDF : {args.output_pdf}")
    print(f"Saved PNG : {output_png}")


if __name__ == "__main__":
    main()
