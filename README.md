# pyTSOPPI_plots

Python re-implementation of the TSO500 run metrics plotting from [TSOPPI](https://github.com/InPreD/TSOPPI) (`plot_run_metrics.R`). Reads per-run metrics TSV files produced by the TSOPPI post-processing pipeline, merges them across runs, and generates a multi-page QC PDF comparing runs side by side with configurable reference guidelines.

<br>

<br>

## Contents

- [Requirements](#requirements)
- [Input files](#input-files)
- [Project structure](#project-structure)
- [Usage](#usage)
- [Running locally](#running-locally)
- [Adding new plots](#adding-new-plots)
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
```

<br>

Install the environment:

```bash
uv sync --locked
```

<br>

<br>

## Input files

The script expects intermediate files produced by the TSOPPI post-processing pipeline. All inputs reside under one directory per run. Using the bundled example run `251211_NB552811_0117_AHTJM2BGYX`:

```
data/TSO500_postprocessing/
└── 251211_NB552811_0117_AHTJM2BGYX/
    └── 251211_NB552811_0117_AHTJM2BGYX_metrics_plots/
        └── intermediate_metrics_files/
            ├── joint_sequencing_QC_file.tsv                 # One row per run
            ├── master_metrics_table.tsv                     # All metrics join table (written by this script)
            ├── 251211_NB552811_0117_AHTJM2BGYX_metrics.tsv  # Current run
            ├── 251210_NB552811_0117_AHTJM2BGYX_metrics.tsv  # Additional run 1
            └── ...                                          # Additional runs [2 .. N]
```

<br>

`joint_sequencing_QC_file.tsv`

One row per sequencing run. Contains run-level metrics used for the cluster density vs. estimated yield scatter plot.

```
RUN_ID                              PCT_PF_READS  PCT_Q30_R1  PCT_Q30_R2  CLUSTER_DENSITY  ESTIMATED_YIELD  CLUSTERS_PASSING_FILTER  RUN_NUMBER
251211_NB552811_0117_AHTJM2BGYX    83.8          88.8        85.4        241.9            112.3            83.8                     1
251210_NB552811_0117_AHTJM2BGYX    83.8          88.8        85.4        241.9            112.3            83.8                     2
...
```

Additionally, rows with `RUN_ID` equal to `LSL_Guideline`, `USL_Guideline`, or `Internal Guideline` encode QC thresholds and are extracted separately to draw reference lines on the plots.

<br>

`<RUN_ID>_metrics.tsv`

One row per sample. Each file covers both DNA (suffix `D01`) and RNA (suffix `R03`) samples for a single run, along with three special guideline rows at the top.

```
SAMPLE_ID            COMPLETED_ALL_STEPS  DNA_CONTAMINATION_SCORE  DNA_CONTAMINATION_P_VALUE  DNA_MEDIAN_INSERT_SIZE  DNA_MEDIAN_EXON_COVERAGE  ...  RNA_MEDIAN_CV_GENE_500X  RNA_TOTAL_PF_READS  ...
LSL_Guideline        NA                   0                        0                          70                      150                             NA                       NA
USL_Guideline        NA                   3106                     0.049                      NA                      NA                              NA                       NA
AU565-D01-C01-X06    TRUE                 790                      0.999                      110                     548                             NA                       NA
AU565-R03-C01-X06    TRUE                 NA                       NA                         NA                      NA                              65.6                     30002600
BT20-D01-C01-X06     TRUE                 ...
...
```

- **Guideline rows** (`LSL_Guideline`, `USL_Guideline`) define lower and upper specification limits. They are separated from sample data internally and used to draw threshold lines on each plot.

<br>

<br>

## Project structure

```
plot_run_metrics.py   # Entry point: argument parsing, table merging, plot dispatch
plot_specs.py         # Declarative plot definitions (no rendering logic)
plots.py              # Reusable plotnine rendering functions
plot_run_metrics.R    # Original R implementation (reference)
example_command.bash  # Local run example with variable setup
pyproject.toml        # Project metadata and dependencies
.github/workflows/
    run_metrics_plot.yml  # CI workflow (uv + ubuntu-latest)
```

**Design principle:** rendering logic lives in `plots.py`; what to plot and how to configure each panel is declared in `plot_specs.py`. Adding a new plot does not require touching the rendering engine.


<br>

<br>

## Usage

```
uv run plot_run_metrics.py <N>                                     \
    <joint_qc_file.tsv>                                            \
    <output.pdf>                                                   \
    <master_metrics_table.tsv>                                     \
    <highlighted_run_id>                                           \
    <create_plots>                                                 \
    <run1_metrics.tsv> [<run2_metrics.tsv> ... <runN_metrics.tsv>] \
    <run1_id>          [<run2_id>          ... <runN_id>]          ;
```

<br>

| Position | Argument | Description |
|----------|----------|-------------|
| 1 | `N` | Number of runs to include |
| 2 | `joint_qc_file` | Path to `joint_sequencing_QC_file.tsv` |
| 3 | `output_pdf` | Output PDF path |
| 4 | `master_metrics_table` | Output path for the merged metrics TSV |
| 5 | `highlighted_run_id` | Run ID to highlight in all plots |
| 6 | `create_plots` | `True` to produce the PDF; `False` to only write the merged TSV |
| 7 to 6+N | `*_metrics.tsv` files | One per run, in order |
| 7+N to 6+2N | run labels | One per run, matching the metrics file order |

The script must be called **twice** in practice: once with `create_plots=False` to write the merged TSV, and once with `create_plots=True` to produce the PDF. This separation allows inspecting the merged data independently of the plot generation step.

<br>

<br>

## Running locally

Set up variables for the current run and any historical runs to include:

```bash
RUN_ID=251211_NB552811_0117_AHTJM2BGYX
TSOPPI_DIR=data/TSO500_postprocessing
METRIC_PLOTS_DIR=${TSOPPI_DIR}/${RUN_ID}/${RUN_ID}_metrics_plots
RESULTS_DIR=results_python/${RUN_ID}
METRICS=${METRIC_PLOTS_DIR}/intermediate_metrics_files/${RUN_ID}_metrics.tsv

# Add one block per historical run to include
RUN_ID_PREV1=251210_NB552811_0117_AHTJM2BGYX
METRICS_PREV1=${METRIC_PLOTS_DIR}/intermediate_metrics_files/${RUN_ID_PREV1}_metrics.tsv

# Aggregate all TSVs and run IDs into space-separated strings
TSVs="${METRICS} ${METRICS_PREV1}"
RUNIDS="${RUN_ID} ${RUN_ID_PREV1}"

mkdir -p ${RESULTS_DIR}

# Step 1: merge metrics tables (no PDF)
uv run plot_run_metrics.py 2                                                    \
    ${METRIC_PLOTS_DIR}/intermediate_metrics_files/joint_sequencing_QC_file.tsv \
    ${RESULTS_DIR}/TSO500_run_metrics.pdf                                       \
    ${METRIC_PLOTS_DIR}/intermediate_metrics_files/master_metrics_table.tsv     \
    ${RUN_ID} False ${TSVs} ${RUNIDS}                                           ;

# Step 2: generate PDF
uv run plot_run_metrics.py 2                                                    \
    ${METRIC_PLOTS_DIR}/intermediate_metrics_files/joint_sequencing_QC_file.tsv \
    ${RESULTS_DIR}/TSO500_run_metrics.pdf                                       \
    ${METRIC_PLOTS_DIR}/intermediate_metrics_files/master_metrics_table.tsv     \
    ${RUN_ID} True ${TSVs} ${RUNIDS}                                            ;
```

A fully parameterised example with up to 10 historical runs is available in [example_command.bash](example_command.bash).

<br>

<br>

## Adding new plots

All plot definitions live in `plot_specs.py` as entries in the `PLOT_SPECS` dictionary. Each entry maps a metric name to a specification that controls which table to query, what columns to use, axis labels, optional guidelines, and overlay annotations. No changes to `plot_run_metrics.py` or `plots.py` are needed unless a genuinely new plot kind is required.

The three patterns in use are shown below.

<br>

### Pattern 1 — bar plot, no guideline

The simplest case: cast one column to float and plot it as a bar chart.

```python
"CLUSTERS_PASSING_FILTER": {
    "plot_number"      : 1,                  # Position in the output PDF
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
    "plot_number"      : 8,
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
    "y_tick_step"      : 100,                    # Fixed distance between y-axis ticks
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
    "plot_number"  : 6,
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

### Adding a new bar plot

1. Choose a key name matching the metric column.
2. Copy Pattern 1 or 2 as a starting point and adjust `source`, `value_spec`, axis labels, and `plot_number`.
3. Set `requires_samples` to `"dna"` or `"rna"` if the metric only exists for one library type.
4. Add a `guideline` block if a QC threshold should be drawn as a horizontal line.

### Adding a new plot kind

If none of the three existing `plot_kind` values fits, add a `render_<kind>` function in `plot_run_metrics.py` and a corresponding dispatch branch in `render_plot()`. The new kind is then available to any entry in `PLOT_SPECS`.

<br>

## Contributors

- [Jaime A Castro-Mondragon](https://jaimicore.github.io/) (jamondra@uio.no)