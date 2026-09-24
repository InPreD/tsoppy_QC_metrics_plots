# pyTSOPPI_plots

Python re-implementation of the TSO500 run metrics plotting from [TSOPPI](https://github.com/InPreD/TSOPPI) (`plot_run_metrics.R`). Reads per-run metrics TSV files produced by the TSOPPI post-processing pipeline, merges them across runs, and generates a multi-page QC PDF comparing runs side by side with configurable reference guidelines.

<br>

<br>

## Contents

- [Requirements](#requirements)
- [Input files](#input-files)
- [Supported workflows](#supported-workflows)
- [Project structure](#project-structure)
- [Usage](#usage)
- [Running locally](#running-locally)
- [Adding new plots](#adding-new-plots)
- [Per-sample QC plots](#per-sample-qc-plots)
- [Modifying per-sample metrics](#modifying-per-sample-metrics)
- [Contributors](#contributors)


<br>

<br>

## Requirements

| Tool | Version |
|------|---------|
| Python | >= 3.12 |
| [uv](https://docs.astral.sh/uv/) | latest |

Python dependencies are declared in `pyproject.toml` and locked via `uv.lock`:

<br>

List of python dependencies:

```
polars
plotnine
pyarrow
numpy
```

<br>

Install the environment:

```bash
uv sync --locked
```

<br>

<br>

## Input files

The script expects intermediate files produced by the TSOPPI post-processing pipeline. Input files are organized by workflow type. The script uses a single set of data files per workflow and filters runs via the `--max_runs` CLI argument.

Using the bundled example data:

**LocalApp** (`data/workflow_tables/`):
```
data/workflow_tables/
├── localapp_plot_joint_sequencing_QC_file.tsv   # Run-level metrics (one row per run)
└── localapp_plot_metrics_table.tsv              # All metrics (sample-level, all runs merged)
```

**DRAGEN** (`data/workflow_tables/`):
```
data/workflow_tables/
├── dragen_plot_joint_sequencing_QC_file.tsv     # Run-level metrics (one row per run)
└── dragen_plot_metrics_table.tsv                # All metrics (sample-level, all runs merged)
```

The `--max_runs` argument controls how many of the most recent runs (lowest `RUN_INDEX` values) are included in each plot.

<br>

<br>

## Supported workflows

Both workflows share the same rendering engine (`plot_QC_metrics.py` and `plots.py`), differing only in which metrics are plotted and how thresholds are defined. The active workflow is selected via the `--workflow_type` CLI argument.

| Workflow | CLI value | Typical use |
|----------|-----------|-------------|
| LocalApp | `localapp` | Illumina LocalApp pipeline output |
| DRAGEN | `dragen` | Illumina DRAGEN pipeline output |

The plot specifications are defined in `specs/plot_specs_workflows.py` as the `PLOT_SPECS` dictionary. Each plot entry contains workflow-specific routing keys (`localapp`, `dragen`) that control whether the plot is generated for each workflow and in what position.

<br>

<br>

## Project structure

```
plot_QC_metrics.py              # Entry point: CLI argument parsing, table loading, plot dispatch
plots.py                        # Reusable plotnine rendering functions
plot_sample_metrics.py          # Per-sample QC panel generator (DNA + RNA, LocalApp + Dragen)
plot_run_metrics.R              # Original R implementation (used as reference)
specs/
    __init__.py                 # Marks specs/ as a Python package
    plot_specs_workflows.py     # Active specs for multi-run plots (run-level QC)
    plot_specs_sample_metrics.py # Active specs for per-sample plots (workflow + seq-type routing)
    plot_specs_Dragen.py        # Plot specifications template for DRAGEN workflow (reference)
    plot_specs_LocalApp.py      # Plot specifications template for LocalApp workflow (reference)
example_command.bash            # Local run example with loop over all run counts
pyproject.toml                  # Project metadata and dependencies
.github/workflows/
    run_metrics_plot.yml        # CI workflow (uv + ubuntu-latest)
```

**Design principle:** Rendering logic lives in `plots.py`; what to plot and how to configure each panel is declared in `specs/plot_specs_workflows.py`—a unified specification file with workflow-specific routing. Each plot entry contains workflow-specific sub-dictionaries (`localapp`, `dragen`) that control whether the plot is generated for each workflow and at what position (index).

**Color palette:** Uses [Tableau 10 palette](https://www.tableau.com/blog/colors-upgrade-tableau-10-56782) for consistent visual distinction across up to 10 runs. If you foresee needing more colors, adapt the variable defined in `plots.py`.

<br>

## Usage

```bash
uv run plot_QC_metrics.py --joint_QC <joint_qc_file.tsv>             \
                          --output_pdf <output.pdf>                  \
                          --metrics_table <master_metrics_table.tsv> \
                          --workflow_type <workflow>                 ;
```

<br>

| Argument | Description |
|----------|-------------|
| `--joint_QC` | Path to `joint_sequencing_QC_file.tsv` (required) |
| `--output_pdf` | Output PDF file path (required); parent directory is created if needed |
| `--metrics_table` | Path to merged metrics TSV file (required) |
| `--workflow_type` | Active workflow: `localapp` or `dragen` (required, case-insensitive) |
| `--max_runs` | Maximum number of most recent runs to include (optional, default: 10, range: 1–10). If fewer runs are available, all are included. |

The script reads the merged metrics table, loads the plot specifications, filters plots based on the selected workflow, and writes a multi-page PDF to the specified output path. The most recent run (lowest `RUN_INDEX`) is automatically highlighted in the plots.

<br>

## Running locally

### Quick start: Single run count with both workflows

```bash
RESULTS_DIR=results_python/LocalApp/1R
mkdir -p ${RESULTS_DIR}

# Run with LocalApp workflow
uv run plot_QC_metrics.py --joint_QC data/workflow_tables/localapp_plot_joint_sequencing_QC_file.tsv \
                          --output_pdf ${RESULTS_DIR}/1R_LocalApp_metric_plots.pdf \
                          --metrics_table data/workflow_tables/localapp_plot_metrics_table.tsv \
                          --max_runs 1 \
                          --workflow_type localapp ;

# Run with Dragen workflow
uv run plot_QC_metrics.py --joint_QC data/workflow_tables/dragen_plot_joint_sequencing_QC_file.tsv \
                          --output_pdf ${RESULTS_DIR}/1R_Dragen_metric_plots.pdf \
                          --metrics_table data/workflow_tables/dragen_plot_metrics_table.tsv \
                          --max_runs 1 \
                          --workflow_type dragen ;
```

<br>

### Batch run: Multiple run counts with loop

A fully parameterized example that runs all combinations of run counts (1R, 2R, 5R, 10R) with both workflows is provided in [example_command.bash](example_command.bash):

```bash
# Review the loop structure
cat example_command.bash

# Execute all runs
bash example_command.bash
```

The script uses:
- `RUN_COUNTS="1R 2R 5R 10R"` — Run count options
- `WORKFLOWS="LocalApp Dragen"` — Workflow options
- Nested loops to generate all combinations
- Organized results directory structure: `results_python/TSO500_${RUN_COUNT}/`
- Separate PDF for each workflow: `*_LocalApp_metric_plots.pdf`, `*_Dragen_metric_plots.pdf`

<br>

<br>

## Adding new plots

All plot definitions live in `specs/plot_specs_workflows.py` as entries in the `PLOT_SPECS` dictionary. Each entry maps a metric name to a specification that controls:
- Which workflows display the plot and at what position (workflow routing via `localapp` and `dragen` sub-dicts)
- Which table to query, what columns to use, axis labels, optional guidelines, and overlay annotations

No changes to `plot_QC_metrics.py` or `plots.py` are needed unless a genuinely new plot kind is required.

The three patterns in use are shown below.

<br>

### Pattern 1 — bar plot, no guideline

The simplest case: cast one column to float and plot it as a bar chart.

```python
"CLUSTERS_PASSING_FILTER": {
    "localapp": {
        "plot": True,   # Include this plot for LocalApp
        "index": 1      # First plot in LocalApp PDF
    },
    "dragen": {
        "plot": False,  # Skip this plot for Dragen
        "index": 0      # Index 0 means not relevant
    },
    
    "plot_kind"        : "bar",              # Plot type
    "source"           : "joint_qc_table",   # table key returned by build_tables()
    "requires_samples" : None,               # Do not display DNA/RNA labels
    "x_var"            : "RUN_ID",
    "y_var"            : "CLUSTERS_PASSING_FILTER",
    "fill_var"         : "RUN_ID",
    "x_lab"            : "Run ID",
    "y_lab"            : "Clusters passing filter (percent)",
    "title"            : "Clusters passing filter",
    "value_spec"       : {"operation": "cast", "column": "CLUSTERS_PASSING_FILTER", "dtype": pl.Float64},
    "na_filter_columns": ["CLUSTERS_PASSING_FILTER"],                                                      # Rows where this column is null or "NA" are dropped
    "cart_ylim"        : (0, 104),
    "skip_if_empty"    : True,               # Skip the panel instead of crashing on empty data
},
```
<br>

### Pattern 2 — bar plot with a guideline line

Adds a horizontal reference line read directly from the guideline rows in the input metrics file (`LSL_Guideline` / `USL_Guideline`).

```python
"DNA_MEDIAN_EXON_COVERAGE": {
    "localapp": {
        "plot": True,
        "index": 8
    },
    "dragen": {
        "plot": False,
        "index": 0
    },
    
    "plot_kind"        : "bar",
    "source"           : "dna_data_table",       # Pre-filtered to DNA samples only
    "requires_samples" : "dna",                  # Panel is skipped when no DNA samples are present
    "x_var"            : "SAMPLE_ID",
    "y_var"            : "DNA_MEDIAN_EXON_COVERAGE",
    "fill_var"         : "RUN",
    "x_lab"            : "Sample ID",
    "y_lab"            : "Median exon coverage",
    "title"            : "[LocalApp core DNA metric] Median exon coverage",
    "value_spec"       : {"operation": "cast", "column": "DNA_MEDIAN_EXON_COVERAGE", "dtype": pl.Float64},
    "na_filter_columns": ["DNA_MEDIAN_EXON_COVERAGE"],
    "y_tick_step"      : 100,                    # Fixed distance between y-axis ticks (supports both int and float, e.g. 0.1, 10, 100)
    "guideline"        : {
        "table"        : "dna_guideline_table",  # Table containing LSL/USL rows
        "sample_id"    : "LSL_Guideline",        # Row to read the threshold from
        "value_spec"   : {"operation": "cast", "column": "DNA_MEDIAN_EXON_COVERAGE", "dtype": pl.Int32},
        "python_cast"  : int,
        "label_prefix" : "LSL_Guideline",        # Annotation text displayed on the plot
        "ann_y_offset" : 200,                    # Vertical nudge for the annotation label
    },
},
```

<br>

### Pattern 3 — scatter plot (contamination)

Used for the DNA contamination score vs. p-value panel. The rendering function (`render_contamination_scatter`) is fixed; only the colour mapping and label column vary between the two contamination panels.

```python
"DNA_CONTAMINATION_P_VALUE_BY_RUN": {
    "localapp": {
        "plot": True,
        "index": 6
    },
    "dragen": {
        "plot": False,
        "index": 0
    },
    
    "plot_kind"    : "contamination_scatter",
    "source"       : "dna_data_table",
    "requires_samples": "dna",
    "color_var"    : "RUN",          # Colour points by run ID
    "label_var"    : None,           # No per-point labels
    "title"        : "[LocalApp core DNA metric] Contamination assessment",
    "skip_if_empty": True,
},
```

<br>

### Workflow routing reference

The workflow routing sub-dicts control plot visibility and ordering:

```python
"localapp": {
    "plot": bool,     # True = include in LocalApp PDF, False = skip
    "index": int      # Plot position in PDF (1-based); 0 = not relevant for this workflow
},
"dragen": {
    "plot": bool,     # True = include in Dragen PDF, False = skip
    "index": int      # Plot position in PDF (1-based); 0 = not relevant for this workflow
}
```

**Important:** Index 0 is treated as "not relevant for this workflow" and is ignored during duplicate checking. Multiple plots can have index 0 for the same workflow without error. Indexes >= 1 must be unique within each workflow.

<br>

### Adding a new bar plot

1. Choose a key name matching the metric column.
2. Copy Pattern 1 or 2 as a starting point and adjust `source`, `value_spec`, axis labels, and workflow routing `index` values.
3. Set `requires_samples` to `"dna"` or `"rna"` if the metric only exists for one library type.
4. Add a `guideline` block if a QC threshold should be drawn as a horizontal line.
5. For each workflow, set `plot` to `True` or `False` and assign a unique `index` >= 1 (or 0 if not relevant).

### Adding a new plot kind

If none of the three existing `plot_kind` values fits, add a `render_<kind>` function in `plot_QC_metrics.py` and a corresponding dispatch branch in `render_plot()`. Register the new kind in `_REQUIRED_FIELDS_BY_KIND`. The new kind is then available to any entry in `PLOT_SPECS`.

<br>

<br>

## Per-sample QC plots

In addition to the multi-run overview (`plot_QC_metrics.py`), `plot_sample_metrics.py` generates a per-sample QC panel for sample-level reporting. It produces a self-contained figure with colour-coded pass/fail rows covering all QC thresholds for a single DNA or RNA sample. The script is a unified Python re-implementation of the original R scripts (`01_plot_sample_metrics_DNA.R`, `plot_sample_metrics_RNA.R`).

A single invocation produces **one PDF + PNG pair per selected workflow**. Output files are named `{base}_{LocalApp}.pdf` / `{base}_{Dragen}.pdf` (and matching `.png`), where `{base}` is the value of `--output_pdf` without its `.pdf` extension. Which workflows are produced is controlled by `--workflow` (`L` = LocalApp only, `D` = Dragen only, `B` = both, default). Workflows with no metrics defined for the requested `--seq_type` are silently skipped.

<br>

### `plot_sample_metrics.py`

Generates a per-sample QC panel (9 rows for DNA, 8 rows for RNA) with colour-coded metrics, threshold comparisons, and red "X" fail indicators.

```bash
uv run plot_sample_metrics.py \
    --input_tsv <path> --sample_id <id> --output_pdf <path> --seq_type [DNA|RNA] \
    [--workflow L|D|B] [--inpred_nomenclature True|False]
```

| Argument | Description |
|----------|-------------|
| `--input_tsv` (required) | TSV metrics table with a `SAMPLE_ID` column and all DNA or RNA QC columns (depending on `--seq_type`). |
| `--sample_id` (required) | Sample ID to extract from the input table and display in the plot title. |
| `--output_pdf` (required) | Output PDF base path (e.g. `results/sample_QC.pdf`). The script appends `_LocalApp` / `_Dragen` before the extension and writes one PDF + PNG pair per selected workflow: `sample_QC_LocalApp.pdf`, `sample_QC_Dragen.pdf`, `sample_QC_LocalApp.png`, `sample_QC_Dragen.png`. |
| `--seq_type` (required) | `DNA` or `RNA` — selects which metrics and QC thresholds to apply. |
| `--workflow` (optional) | `L` (LocalApp only), `D` (Dragen only), or `B` (both; default: `B`). Controls which workflow output files are produced. |
| `--inpred_nomenclature` (optional) | `True` or `False` (default: `False`). When `True`, the 13th character of `sample_id` is decoded as an InPreD sample-type code (requires a 19-character ID; e.g. `C` → cell line, `T` → tumour). |

#### DNA metrics (9 rows)

```bash
uv run plot_sample_metrics.py \
    --input_tsv <path> --sample_id <id> --output_pdf <path> --seq_type DNA
```

| Row | Metric | Threshold | Colour logic |
|-----|--------|-----------|---------------|
| 1 | Pipeline completed | — | boolean (green = complete) |
| 2 | Contamination | score > 3106 **and** p-value > 0.049 | boolean (amber = contaminated) |
| 3 | Median exon coverage | ≥ 150 | FailPass |
| 4 | Exon bases with coverage ≥50× (%) | ≥ 90 | FailPass |
| 5 | Median insert size | ≥ 70 | FailPass |
| 6 | Aligned reads (%) | ≥ 70 | FailPass |
| 7 | Usable MSI sites | ≥ 40 | FailPass |
| 8 | Median raw bin count per CNV target | ≥ 1 | FailPass |
| 9 | CNV coverage MAD | ≤ 0.21 | PassFail |

#### RNA metrics (8 rows)

```bash
uv run plot_sample_metrics.py \
    --input_tsv <path> --sample_id <id> --output_pdf <path> --seq_type RNA
```

| Row | Metric | Threshold | Colour logic |
|-----|--------|-----------|---------------|
| 1 | Pipeline completed | — | boolean (green = complete) |
| 2 | CV for genes at median coverage >500× | ≤ 93 | PassFail |
| 3 | Total on-target reads (millions) | ≥ 9 | FailPass |
| 4 | Median insert size | ≥ 80 | FailPass |
| 5 | Chimeric reads (%) | no reference | no threshold bar |
| 6 | On-target reads (%) | no reference | no threshold bar |
| 7 | Scaled median gene coverage | no reference | no threshold bar |
| 8 | Total reads passing filters (millions) | no reference | no threshold bar |

> **Note on RNA read counts:** `RNA_TOTAL_ON_TARGET_READS` and `RNA_TOTAL_PF_READS` are stored as raw counts in the input TSV. The script divides them by 1 000 000 internally, so thresholds (9.0 and 30.0) are expressed in millions.

<br>

### Visual indicators

Each row shows:
- **Colour-coded bar**: Green indicates pass, orange indicates fail.
- **Observed value**: Displayed as a number inside the bar with a thin **light red vertical line** marking its position on the scale.
- **Red "X" marker**: Large red "X" appears in the conclusion box (right side) when a metric fails the QC threshold, making it immediately visible in quick visual scans.

<br>

### Examples

Using the bundled example data:

**Passing DNA sample (LocalApp only):**
```bash
mkdir -p results_python/sample_metrics

uv run plot_sample_metrics.py \
    --input_tsv data/workflow_tables/localapp_plot_metrics_table.tsv \
    --sample_id IPD0195-D01-P91-A25-0001_P \
    --output_pdf results_python/sample_metrics/IPD0195-D01-P91-A25-0001_P_DNA_QC.pdf \
    --seq_type DNA \
    --workflow L
# Produces:
#   results_python/sample_metrics/IPD0195-D01-P91-A25-0001_P_DNA_QC_LocalApp.pdf  (+ .png)
```

**Passing DNA sample (Dragen only):**
```bash
uv run plot_sample_metrics.py \
    --input_tsv data/workflow_tables/dragen_plot_metrics_table.tsv \
    --sample_id TVD0018 \
    --output_pdf results_python/sample_metrics/TVD0018_DNA_QC.pdf \
    --seq_type DNA \
    --workflow D
# Produces:
#   results_python/sample_metrics/TVD0018_DNA_QC_Dragen.pdf  (+ .png)
```

**Passing DNA sample (both workflows):**
```bash
uv run plot_sample_metrics.py \
    --input_tsv data/workflow_tables/localapp_plot_metrics_table.tsv \
    --sample_id IPD0195-D01-P91-A25-0001_P \
    --output_pdf results_python/sample_metrics/IPD0195-D01-P91-A25-0001_P_DNA_QC.pdf \
    --seq_type DNA
# Produces (--workflow defaults to B):
#   results_python/sample_metrics/IPD0195-D01-P91-A25-0001_P_DNA_QC_LocalApp.pdf  (+ .png)
#   results_python/sample_metrics/IPD0195-D01-P91-A25-0001_P_DNA_QC_Dragen.pdf   (+ .png)
```

**Passing RNA sample (LocalApp only):**
```bash
uv run plot_sample_metrics.py \
    --input_tsv data/workflow_tables/localapp_plot_metrics_table.tsv \
    --sample_id IPD0195-R03-P91-A25-0001_P \
    --output_pdf results_python/sample_metrics/IPD0195-R03-P91-A25-0001_P_RNA_QC.pdf \
    --seq_type RNA \
    --workflow L
# Produces:
#   results_python/sample_metrics/IPD0195-R03-P91-A25-0001_P_RNA_QC_LocalApp.pdf  (+ .png)
```

**Failing DNA sample (multiple QC failures with red "X" markers):**
```bash
uv run plot_sample_metrics.py \
    --input_tsv data/workflow_tables/localapp_plot_metrics_table.tsv \
    --sample_id IPD1881-D01-D01-A25_P \
    --output_pdf results_python/sample_metrics/IPD1881-D01-D01-A25_P_DNA_QC.pdf \
    --seq_type DNA \
    --workflow L
# Produces:
#   results_python/sample_metrics/IPD1881-D01-D01-A25_P_DNA_QC_LocalApp.pdf  (+ .png)
```

Both scripts automatically create `.png` files alongside the `.pdf` files.

<br>

<br>

## Modifying per-sample metrics

All per-sample metric row definitions live in `specs/plot_specs_sample_metrics.py` as entries in the `SAMPLE_METRICS_SPECS` dictionary. No changes to `plot_sample_metrics.py` are needed unless a genuinely new row rendering type is required.

Each entry controls:
- Which workflows display the row and at what position (`localapp` / `dragen` sub-dicts with `plot`, `index`, and `seq_type` keys)
- Visual rendering mode (`qc_type`: `TrueFalse`, `FailPass`, `PassFail`, `NoThold`)
- Display range, pass/fail threshold, label text, and optional value scaling

### Workflow routing sub-dict format

```python
"METRIC_NAME": {
    "localapp": {
        "plot":     bool,        # True = include this row in the LocalApp plot
        "index":    int,         # Row position (1-based, top-to-bottom); 0 = not relevant
        "seq_type": list[str],   # ["dna"], ["rna"], or ["dna", "rna"]
    },
    "dragen": {                  # Omit if this metric does not exist in Dragen
        "plot":     bool,
        "index":    int,
        "seq_type": list[str],
    },

    # Row rendering fields
    "qc_type":         str | None,   # "TrueFalse" | "FailPass" | "PassFail" | "NoThold" | None
    "label":           str | None,   # Display label (None for helper entries)
    "threshold_label": str | None,   # Threshold annotation shown below the label
    "pass_value":      float | None, # TrueFalse only: the value that maps to PASS (green)
    "variable":        str | None,   # TSV column name; None for derived boolean rows
    "lowest":          float | None, # Left edge of the display bar
    "threshold":       float | None, # Pass/fail boundary
    "highest":         float | None, # Right edge of the display bar
    "float_format":    bool,         # True = display values with 1 decimal place
    "scale_factor":    int,          # Optional; divide raw value before display (e.g. 1_000_000)
},
```

### Adding a new metric row

1. Add an entry to `SAMPLE_METRICS_SPECS` with the correct workflow sub-dict(s), `qc_type`, label, threshold, and display range.
2. Indexes must be unique per workflow × seq_type combination.
3. No changes to `plot_sample_metrics.py` are required.

### Adding Dragen-specific thresholds or metrics

- To use **different thresholds** for Dragen, update the top-level `threshold` field. If LocalApp and Dragen need different threshold values, the threshold can be moved into the workflow sub-dict in a future extension.
- To **enable an existing entry for Dragen**, add or update the `"dragen"` sub-dict.
- To add a **Dragen-exclusive metric**, create a new entry with only a `"dragen"` sub-dict (no `"localapp"` key needed).

<br>

<br>

## Contributors

- [Jaime A Castro-Mondragon](https://jaimicore.github.io/) (jamondra@uio.no)