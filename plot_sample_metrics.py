#!/usr/bin/env python3
"""
Plot per-sample QC metrics from a TSO500 LocalApp metrics table.

Supports DNA and RNA sample types via --seq_type.  The seq-type determines
which metrics, thresholds, and QC logic are applied; the rendering pipeline
is shared and generic.

Usage
-----
    python3 01_plot_sample_metrics.py \\
        --input_tsv  <path>           \\
        --sample_id  <id>             \\
        --output_pdf <path>           \\
        --seq_type   [DNA|RNA]        \\
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

from specs.plot_specs_sample_metrics import SAMPLE_METRICS_SPECS

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
# Configuration builders  (derived from SAMPLE_METRICS_SPECS)
# ---------------------------------------------------------------------------

def _get_sorted_specs(seq_type: str, workflow: str) -> list[tuple[str, dict]]:
    """
    Return (key, spec) pairs for all plotted rows of *seq_type* / *workflow*,
    sorted by index.
    """
    wkey  = workflow.lower()
    stype = seq_type.lower()
    return sorted(
        [
            (k, v) for k, v in SAMPLE_METRICS_SPECS.items()
            if wkey in v
            and v[wkey].get("plot", False)
            and stype in v[wkey].get("seq_type", [])
        ],
        key=lambda x: x[1][wkey]["index"],
    )


def _resolve_threshold(
    spec: dict, thresholds: dict[str, dict[str, float]]
) -> float:
    """
    Look up the QC threshold for a row from the runtime-loaded TSV thresholds.

    - FailPass  → LOWER_THRESHOLD row  ("lower" key)
    - PassFail  → UPPER_THRESHOLD row  ("upper" key)
    - NoThold / TrueFalse / helpers → np.nan

    If the variable has no threshold in the TSV (e.g. DNA_PCT_ALIGNED_READS),
    np.nan is returned and the row renders without a threshold line.
    scale_factor is applied when present (e.g. raw read counts → millions).
    """
    source = spec.get("threshold_source")
    if source is None:
        return np.nan

    variable = spec.get("variable")
    if variable is None:
        return np.nan

    raw_thr = thresholds.get(source, {}).get(variable, np.nan)
    if np.isnan(raw_thr):
        raw_thr = spec.get("threshold_fallback", np.nan)
    if np.isnan(raw_thr):
        return np.nan

    factor = spec.get("scale_factor", 1)
    if factor != 1:
        return round(raw_thr / factor, 1)
    return raw_thr


def build_plot_rows(
    seq_type: str, workflow: str, thresholds: dict[str, dict[str, float]]
) -> list[dict]:
    """
    Build the plot_rows list consumed by build_layer_data() for the given
    seq_type and workflow.

    Each entry has:
      label, threshold_label, qc_type, pass_value, metric (or None for TrueFalse rows).
    """
    rows = []
    for _, spec in _get_sorted_specs(seq_type, workflow):
        row: dict = {
            "label":           spec["label"],
            "threshold_label": spec["threshold_label"],
            "qc_type":         spec["qc_type"],
            "pass_value":      spec["pass_value"],
            "metric":          None,
        }
        if spec["qc_type"] != "TrueFalse":
            row["metric"] = {
                "variable":  spec["variable"],
                "lowest":    spec["lowest"]  if spec["lowest"]  is not None else np.nan,
                "threshold": _resolve_threshold(spec, thresholds),
                "highest":   spec["highest"] if spec["highest"] is not None else np.nan,
            }
        rows.append(row)
    return rows


def build_float_metrics(seq_type: str, workflow: str) -> set[str]:
    """
    Return the set of variable names that use 1 d.p. (float) display for the
    given seq_type and workflow.
    """
    wkey  = workflow.lower()
    stype = seq_type.lower()
    return {
        spec["variable"]
        for spec in SAMPLE_METRICS_SPECS.values()
        if spec.get("float_format") and spec["variable"] is not None
        and wkey in spec
        and spec[wkey].get("plot", False)
        and stype in spec[wkey].get("seq_type", [])
    }


def _get_loadable_variables(seq_type: str, workflow: str) -> list[str]:
    """
    Return all TSV column names that must be loaded for the given seq_type and
    workflow.  Includes plotted rows and any helper entries marked with load=True.
    """
    wkey  = workflow.lower()
    stype = seq_type.lower()
    return [
        spec["variable"]
        for spec in SAMPLE_METRICS_SPECS.values()
        if spec["variable"] is not None
        and wkey in spec
        and stype in spec[wkey].get("seq_type", [])
        and (spec[wkey].get("plot", False) or spec[wkey].get("load", False))
    ]


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
    tsv_path: str, sample_id: str, seq_type: str, workflow: str
) -> dict[str, float]:
    """Load the metrics TSV and return {variable: float_value} for the given sample."""
    variables = _get_loadable_variables(seq_type, workflow)
    df        = pl.read_csv(tsv_path, separator="\t", infer_schema=False)
    filtered  = df.filter(pl.col("SAMPLE_ID") == sample_id)
    if filtered.is_empty():
        sys.exit(
            f"ERROR: sample ID '{sample_id}' not found in '{tsv_path}'.\n"
            f"Available IDs: {df['SAMPLE_ID'].to_list()}"
        )
    row = filtered.row(0, named=True)
    return {var: _to_float(row.get(var)) for var in variables}


def load_thresholds(tsv_path: str) -> dict[str, dict[str, float]]:
    """
    Extract QC threshold values from LOWER_THRESHOLD and UPPER_THRESHOLD rows.

    Returns::

        {
            "lower": {"DNA_MEDIAN_EXON_COVERAGE": 150.0, ...},  # LSL_Guideline row
            "upper": {"DNA_CONTAMINATION_SCORE":  3106.0, ...},  # USL_Guideline row
        }

    Only the first occurrence of each RECORD_TYPE is used; thresholds are
    consistent across all runs within the same TSV file.  NA values are omitted.
    """
    df = pl.read_csv(tsv_path, separator="\t", infer_schema=False)
    result: dict[str, dict[str, float]] = {"lower": {}, "upper": {}}

    for record_type, key in [("LOWER_THRESHOLD", "lower"), ("UPPER_THRESHOLD", "upper")]:
        rows = df.filter(pl.col("RECORD_TYPE") == record_type)
        if rows.is_empty():
            continue
        row = rows.row(0, named=True)
        for col, val in row.items():
            fval = _to_float(val)
            if not np.isnan(fval):
                result[key][col] = fval

    return result


# ---------------------------------------------------------------------------
# Observed values computation (seq-type-specific)
# ---------------------------------------------------------------------------

def _compute_obs_dna(
    raw: dict[str, float], workflow: str,
    thresholds: dict[str, dict[str, float]]
) -> list[float]:
    """
    Derive observed values for the DNA plot rows (in index order).

    The contamination boolean is derived from DNA_CONTAMINATION_SCORE and
    DNA_CONTAMINATION_P_VALUE using the UPPER_THRESHOLD values from the TSV.
    All other rows read their value directly from the raw dict.
    """
    sorted_rows = _get_sorted_specs("DNA", workflow)
    obs: list[float] = [np.nan] * len(sorted_rows)

    for i, (key, spec) in enumerate(sorted_rows):
        if key == "DNA_CONTAMINATION":
            score_thr = thresholds["upper"].get("DNA_CONTAMINATION_SCORE", np.nan)
            pval_thr  = thresholds["upper"].get("DNA_CONTAMINATION_P_VALUE", np.nan)
            score     = raw.get("DNA_CONTAMINATION_SCORE", np.nan)
            pval      = raw.get("DNA_CONTAMINATION_P_VALUE", np.nan)
            obs[i]    = 1.0 if (
                not np.isnan(score)     and not np.isnan(pval)
                and not np.isnan(score_thr) and not np.isnan(pval_thr)
                and score > score_thr and pval > pval_thr
            ) else 0.0
        elif spec["variable"] is not None:
            obs[i] = raw.get(spec["variable"], np.nan)

    return obs


def _compute_obs_rna(raw: dict[str, float], workflow: str) -> list[float]:
    """
    Derive observed values for the RNA plot rows (in index order).

    Entries with a scale_factor in the spec are divided and rounded to 1 d.p.
    (e.g. raw read counts → millions).
    """
    sorted_rows = _get_sorted_specs("RNA", workflow)
    obs: list[float] = [np.nan] * len(sorted_rows)

    for i, (_, spec) in enumerate(sorted_rows):
        if spec["variable"] is None:
            continue
        val    = raw.get(spec["variable"], np.nan)
        factor = spec.get("scale_factor", 1)
        if factor != 1 and not np.isnan(val):
            obs[i] = round(val / factor, 1)
        else:
            obs[i] = val

    return obs


def compute_plot_obs_values(
    raw: dict[str, float], seq_type: str, workflow: str,
    thresholds: dict[str, dict[str, float]]
) -> list[float]:
    """Dispatch to the seq-type-specific observed-values computation."""
    if seq_type == "DNA":
        return _compute_obs_dna(raw, workflow, thresholds)
    return _compute_obs_rna(raw, workflow)


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
    if np.isnan(thr_scaled):
        # No threshold in TSV — no fail region exists, so the whole bar is the
        # pass colour (color_right for FailPass = green; color_left for PassFail = green).
        # NoThold rows always enter this branch because color_left == color_right == white.
        pass_color = color_right if qc_type in ("FailPass", "NoThold") else color_left
        bands.append({"xmin": X_BAR_LEFT, "xmax": X_BAR_RIGHT, "ymin": yb, "ymax": yt, "fill": pass_color})
    else:
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
    args     = parse_args()
    seq_type = args.seq_type

    use_inpred  = args.inpred_nomenclature == "True"
    sample_type = get_sample_type(args.sample_id, use_inpred)
    base_path   = args.output_pdf.rsplit(".", 1)[0]
    thresholds  = load_thresholds(args.input_tsv)

    for workflow in ("LocalApp", "Dragen"):
        plot_rows = build_plot_rows(seq_type, workflow, thresholds)
        if not plot_rows:
            continue

        raw_values    = load_and_extract(args.input_tsv, args.sample_id, seq_type, workflow)
        obs_values    = compute_plot_obs_values(raw_values, seq_type, workflow, thresholds)
        float_metrics = build_float_metrics(seq_type, workflow)

        title   = f'QC for {seq_type} sample "{args.sample_id}" [{workflow}]\n{sample_type}'
        pdf_out = f"{base_path}_{workflow}.pdf"
        png_out = f"{base_path}_{workflow}.png"

        p = make_plot(obs_values, title, plot_rows, float_metrics)
        p.save(pdf_out, width=9, height=7, units="in", verbose=False)
        p.save(png_out, width=9, height=7, units="in", dpi=150, verbose=False)

        print(f"Saved PDF : {pdf_out}")
        print(f"Saved PNG : {png_out}")


if __name__ == "__main__":
    main()
