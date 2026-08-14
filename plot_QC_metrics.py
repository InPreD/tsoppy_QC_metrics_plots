#!/usr/bin/env python3
"""
TSO500 Run Metrics Plotting Script

This script processes and visualizes quality metrics from TSO500 sequencing runs.
Converts metrics data from multiple runs into consolidated outputs and generates
comprehensive QC plots.

Command-line parameters:
    --joint_QC path/to/joint_sequencing_QC_file.tsv
    --output_folder path/to/output/directory
    --metrics_table path/to/master_metrics_table.tsv
    --run_id RUN_ID
"""

import argparse
import logging
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import polars as pl
from matplotlib.backends.backend_pdf import PdfPages
from plotnine import *

# PLOTS_SPECS is a dictionary containing the specifications for each plot to be generated
# including the type of plot, the data source, and various parameters for rendering.
from plots import AXIS_LINE_COLOR, TABLEAU_20, plot_contamination_scatter, plot_tsoppy_barplot
from specs.plot_specs_workflows import PLOT_SPECS


logger = logging.getLogger(__name__)

VERSION_STRING = "0.3:22-03-10"
ANGLE_X_NAMES = 90


# Minimum required fields per plot_kind, used by validate_plot_specs().
_REQUIRED_COMMON_FIELDS  = {"plot_kind", "source", "title"}
_REQUIRED_FIELDS_BY_KIND = {
    "bar"                      : {"x_var", "y_var", "fill_var", "x_lab", "y_lab", "value_spec"},
    "cluster_density_scatter"  : set(),
    "contamination_scatter"    : {"color_var", "label_var"},
}

SUPPORTED_WORKFLOW = ["Dragen", "LocalApp"]


_REQUIRED_WORKFLOW_KEYS = {"plot", "index"}
_SUPPORTED_WORKFLOW_KEYS = {w.lower() for w in SUPPORTED_WORKFLOW}


def validate_plot_specs(plot_specs: dict) -> None:
    """
    Validate that every entry in PLOT_SPECS contains the minimal required fields.

    Checks:
      - Workflow routing sub-dicts ("localapp", "dragen") are present and well-formed:
          - Both keys must be present.
          - "plot" must be a bool.
          - "index" must be a non-negative integer.
      - Common fields present in every spec (plot_kind, source, title).
      - plot_kind is one of the supported renderer keys.
      - Kind-specific required fields are present for the given plot_kind.
      - No duplicate indexes within each workflow.

    All structural errors are collected and reported together per spec.

    Raises:
        KeyError: For duplicate indexes or other structural validation errors.
    """
    valid_kinds = set(_REQUIRED_FIELDS_BY_KIND)
    structural_errors = []

    # === Phase 1: Check for duplicate indexes per workflow ===
    for wf_key in _SUPPORTED_WORKFLOW_KEYS:
        seen_indexes = {}
        for spec_name, spec in plot_specs.items():
            if wf_key in spec:
                wf_entry = spec[wf_key]
                if isinstance(wf_entry, dict) and "index" in wf_entry:
                    idx = wf_entry["index"]
                    # Only check for duplicates among indexes >= 1 (0 means plot is not relevant)
                    if idx >= 1:
                        if idx in seen_indexes:
                            error_msg = (
                                f"Workflow '{wf_key}': Duplicate index {idx} found in specs "
                                f"'{seen_indexes[idx]}' and '{spec_name}'"
                            )
                            structural_errors.append(error_msg)
                            logger.error(error_msg)
                        else:
                            seen_indexes[idx] = spec_name

    # Raise early if duplicate indexes found
    if structural_errors:
        raise KeyError(f"Found {len(structural_errors)} duplicate index error(s):\n\n" + "\n\n".join(structural_errors))

    # === Phase 2: Existing structural validation ===
    for spec_name, spec in plot_specs.items():
        spec_errors = []

        # --- Validate workflow routing sub-dicts ---
        for wf_key in _SUPPORTED_WORKFLOW_KEYS:
            if wf_key not in spec:
                spec_errors.append(f"  Missing workflow routing key: '{wf_key}'")
                continue

            wf_entry = spec[wf_key]
            missing_wf_fields = _REQUIRED_WORKFLOW_KEYS - set(wf_entry)
            if missing_wf_fields:
                spec_errors.append(
                    f"  Workflow '{wf_key}' is missing fields: {', '.join(sorted(missing_wf_fields))}"
                )
                continue

            if not isinstance(wf_entry["plot"], bool):
                spec_errors.append(
                    f"  Workflow '{wf_key}': 'plot' must be a bool, "
                    f"got {type(wf_entry['plot']).__name__!r}"
                )
            if not isinstance(wf_entry["index"], int) or wf_entry["index"] < 0:
                spec_errors.append(
                    f"  Workflow '{wf_key}': 'index' must be a non-negative integer, "
                    f"got {wf_entry['index']!r}"
                )

        # --- Validate common plot fields ---
        spec_keys = set(spec)
        required_fields = _REQUIRED_COMMON_FIELDS.copy()

        # Check if plot_kind exists and is valid, then add kind-specific fields to requirements
        plot_kind_error = None
        if "plot_kind" in spec_keys:
            plot_kind = spec["plot_kind"]
            if plot_kind in valid_kinds:
                required_fields.update(_REQUIRED_FIELDS_BY_KIND[plot_kind])
            else:
                # Invalid plot_kind: log it but continue to collect all missing fields
                valid_kinds_str = ", ".join(sorted(valid_kinds))
                plot_kind_error = (
                    f"  plot_kind '{plot_kind}' is not recognized. "
                    f"Valid plot types: {valid_kinds_str}"
                )

        missing_fields = required_fields - spec_keys
        if plot_kind_error:
            spec_errors.append(plot_kind_error)
        if missing_fields:
            spec_errors.append(f"  Missing fields: {', '.join(sorted(missing_fields))}")

        if spec_errors:
            error_msg = f"Error in '{spec_name}':\n" + "\n".join(spec_errors)
            structural_errors.append(error_msg)
            logger.error(error_msg)

    if structural_errors:
        error_summary = "\n\n".join(structural_errors)
        raise KeyError(f"Found {len(structural_errors)} validation error(s) in specs:\n\n{error_summary}")


def parse_args(argv):
    """
    Parse the named command-line arguments used for plotting.

    The script expects a precomputed joint QC table and a precomputed
    merged metrics table, both supplied explicitly on the command line.
    The parent directory of output_pdf is created if it does not exist.
    """
    parser = argparse.ArgumentParser(description="Generate TSO500 QC plots from precomputed metrics tables")
    parser.add_argument("--joint_QC",      required=True, dest="joint_qc_file")
    parser.add_argument("--metrics_table", required=True, dest="metrics_table")
    parser.add_argument("--run_id",        required=True, dest="highlighted_run_id")
    parser.add_argument("--output_pdf",    required=True, dest="output_pdf")
    parser.add_argument("--workflow_type", required=True, dest="workflow_type", choices=[w.lower() for w in SUPPORTED_WORKFLOW], type=str.lower, help=f"Workflow type to plot. Supported values: {', '.join(SUPPORTED_WORKFLOW)}",)
    args = parser.parse_args(argv)

    return {
        "joint_qc_file"      : args.joint_qc_file,
        "output_pdf"         : args.output_pdf,
        "metrics_table"      : args.metrics_table,
        "highlighted_run_id" : args.highlighted_run_id,
        "workflow_type"      : args.workflow_type,
    }


def valid_metric_expr(column_name):
    """
    Treat nulls and literal NA strings as missing metric values.
    """
    return pl.col(column_name).is_not_null() & (pl.col(column_name).cast(pl.Utf8) != "NA")


def build_filter_expression(filter_spec):
    """Translate a simple filter description into a Polars expression."""
    column_expr = pl.col(filter_spec["column"])
    if "contains" in filter_spec:
        return column_expr.cast(pl.Utf8).str.contains(filter_spec["contains"])
    if "equals" in filter_spec:
        return column_expr == filter_spec["equals"]
    if "not_equals" in filter_spec:
        return column_expr != filter_spec["not_equals"]
    raise ValueError(f"Unsupported filter specification: {filter_spec}")


def build_value_expression(value_spec, alias_name=None):
    """Create a Polars expression for a raw, scaled, or derived plotting metric."""
    operation = value_spec.get("operation", "cast")

    if operation == "cast":
        expression = pl.col(value_spec["column"])
        if value_spec.get("dtype") is not None:
            expression = expression.cast(value_spec["dtype"])

    elif operation == "divide":
        expression = pl.col(value_spec["column"]).cast(value_spec.get("dtype", pl.Float64)) / value_spec["divisor"]
    
    elif operation == "ratio":
        expression = pl.col(value_spec["numerator"]).cast(value_spec.get("dtype", pl.Float64))
        numerator_divisor = value_spec.get("numerator_divisor")
        if numerator_divisor is not None:
            expression = expression / numerator_divisor
        expression = expression / pl.col(value_spec["denominator"]).cast(value_spec.get("dtype", pl.Float64))
    else:
        raise ValueError(f"Unsupported value operation: {operation}")

    if alias_name is not None:
        return expression.alias(alias_name)
    return expression


def get_guideline_value(tables, guideline_spec):
    """Extract one scalar guideline value and its plotting parameters."""
    table = tables[guideline_spec["table"]]
    if table.height == 0:
        return None

    id_column = guideline_spec.get("id_column", "SAMPLE_ID")
    if id_column not in table.columns:
        return None

    filtered_table = table.filter(pl.col(id_column) == guideline_spec["sample_id"])
    if filtered_table.height == 0:
        return None

    expression = build_value_expression(guideline_spec["value_spec"])
    try:
        value = filtered_table.select(expression).head(1).item()
    except pl.exceptions.ColumnNotFoundError:
        return None

    if value is None:
        return None

    python_cast = guideline_spec.get("python_cast")
    if python_cast is not None:
        value = python_cast(value)

    return {
        "value": value,
        "label": guideline_spec.get("label", f"{guideline_spec['label_prefix']}: {value}"),
        "alpha": guideline_spec.get("alpha"),
        "color": guideline_spec.get("color"),
        "ann_y_offset": guideline_spec.get("ann_y_offset", 0),
    }


def compute_cart_ylim(spec, plot_data):
    """
    Resolve static or data-driven y-axis limits for a bar plot.
     - If "cart_ylim_dynamic" is not specified in the plot specification, use the static limits defined in "cart_ylim".
     - If "cart_ylim_dynamic" is specified with mode "max_plus", calculate the maximum value of the specified column in the plot data and add the defined offset to determine the upper limit
    """
    dynamic_ylim = spec.get("cart_ylim_dynamic")

    # Use static limits if no dynamic y-limit specification is provided
    if dynamic_ylim is None:
        return spec.get("cart_ylim")

    # Handle dynamic y-limit calculation based on the specified mode
    # 'max_plus' mode calculates the upper limit as the maximum value of the specified column plus an offset, 
    # the lower limit can be set to a defined value or default to 0.
    if dynamic_ylim["mode"] == "max_plus":
        max_value = plot_data[dynamic_ylim["column"]].max()
        return (dynamic_ylim.get("lower", 0), max_value + dynamic_ylim.get("offset", 0))

    raise ValueError(f"Unsupported dynamic y-limit mode: {dynamic_ylim['mode']}")


def prepare_bar_plot_data(table, spec):
    """
    Apply per-plot NULL/NA filtering and column derivation before converting to pandas.
    This allows the plot specifications to define which rows and columns are relevant for each plot, and to derive new metrics for plotting without modifying the original data tables.
    """

    filters = []
    
    # Remove rows with null or NA values in the specified columns for this plot
    for column_name in spec.get("na_filter_columns", []):
        filters.append(valid_metric_expr(column_name))
    
    # Apply additional filters defined in the plot specification, which can include conditions like "contains", "equals", or "not_equals" for specific columns.
    for filter_spec in spec.get("filters", []):
        filters.append(build_filter_expression(filter_spec))

    # Apply all filters to the table, combining them with logical AND.
    # If no filters are specified, the original table is returned unchanged.
    if filters:
        combined_filter = filters[0]
        for filter_expression in filters[1:]:
            combined_filter = combined_filter & filter_expression
        table = table.filter(combined_filter)

    table = table.with_columns(build_value_expression(spec["value_spec"], spec["y_var"]))
    return table.to_pandas()


def save_plot(pdf_handle, plot):
    """Render a plotnine plot to the shared PDF and close the matplotlib figure."""
    figure = plot.draw()
    pdf_handle.savefig(figure, bbox_inches="tight")
    plt.close(figure)



# The main plotting function for bar plots, which is used by multiple plot specifications
# All bar plots are created using this function, and the appearance and content of the plots are controlled by the plot specifications defined in plot_specs.py
def render_bar_plot(pdf_handle, spec, tables):
    """
    Render one generic bar plot from its specification entry.
    """

    # Some data requires DNA, RNA or QC
    plot_data = prepare_bar_plot_data(tables[spec["source"]], spec)
    if spec.get("skip_if_empty") and plot_data.empty:
        return

    # Add a horizontal guideline if specified in the plot specification, which is typically used to indicate a QC threshold or target value.
    guideline = None
    if spec.get("guideline") is not None:
        guideline = get_guideline_value(tables, spec["guideline"])

    plot = plot_tsoppy_barplot(
        data          = plot_data,
        x_var         = spec["x_var"],
        y_var         = spec["y_var"],
        fill_var      = spec["fill_var"],
        guide_title   = spec.get("guide_title", "Run"),
        x_lab         = spec["x_lab"],
        y_lab         = spec["y_lab"],
        cart_ylim     = compute_cart_ylim(spec, plot_data),
        title         = spec["title"],
        x_lab_angle   = spec.get("x_lab_angle", ANGLE_X_NAMES),
        fig_size      = spec.get("fig_size", (15, 6)),
        alpha_value   = spec.get("alpha_value", 0.8),
        hline_y       = None if guideline is None else guideline["value"],
        hline_alpha   = None if guideline is None else guideline["alpha"],
        hline_color   = None if guideline is None else guideline["color"],
        hline_label   = "" if guideline is None else guideline["label"],
        ann_y_offset  = 0 if guideline is None else guideline["ann_y_offset"],
        y_tick_step   = spec.get("y_tick_step"),
    )

    # Some metrics (e.g., contamination labels) require overlaying text annotations on the bars
    # This is specified in the plot specification under "overlay_labels".
    overlay_labels = spec.get("overlay_labels")
    if overlay_labels is not None:
        plot = plot + geom_text(
            aes(label=overlay_labels["label_column"]),
            angle=overlay_labels.get("angle", 0),
            va=overlay_labels.get("va"),
            ha=overlay_labels.get("ha"),
            nudge_x=overlay_labels.get("nudge_x", 0),
            nudge_y=overlay_labels.get("nudge_y", 0),
            size=overlay_labels.get("size", 10),
        )

    save_plot(pdf_handle, plot)


# We use it only once for the ESTIMATED_YIELD plot
def render_cluster_density_scatter(pdf_handle, spec, tables):
    """
    Render the sequencer cluster density versus estimated yield scatter plot.
    """
    plot_table = (
        tables[spec["source"]]
        .filter(valid_metric_expr("ESTIMATED_YIELD"))
        .with_columns([
            pl.col("CLUSTER_DENSITY").cast(pl.Float64),
            pl.col("ESTIMATED_YIELD").cast(pl.Float64),
            (pl.col("RUN_NUMBER").cast(pl.Utf8) + " - " + pl.col("RUN_ID")).alias("RUN_LABEL"),
        ])
    )

    if spec.get("skip_if_empty") and plot_table.height == 0:
        return

    plot_data = plot_table.to_pandas()
    plot = (
        ggplot(plot_data, aes(x="CLUSTER_DENSITY", y="ESTIMATED_YIELD"))
        + geom_point(aes(color="RUN_LABEL"), alpha=0.8, size=4.5)
        + scale_color_manual(values=TABLEAU_20)
        + geom_text(aes(label="RUN_NUMBER"), alpha=0.8, size=10, nudge_y=-5)
        + guides(color=guide_legend("Run"))
        + coord_cartesian(xlim=(0, 400), ylim=(0, 150))
        + xlab("Cluster density")
        + ylab("Estimated yield")
        + ggtitle(spec["title"])
        + scale_x_continuous(breaks=[index * 50 for index in range(9)], minor_breaks=[index * 25 for index in range(17)])
        + scale_y_continuous(breaks=[index * 10 for index in range(16)], minor_breaks=[index * 5 for index in range(31)])
        + theme_minimal()
        + theme(
            axis_text_x=element_text(angle=spec.get("x_lab_angle", ANGLE_X_NAMES), vjust=0.5, hjust=0.5, size=7),
            figure_size=spec.get("fig_size", (15, 6)),
            axis_line=element_line(color=AXIS_LINE_COLOR),
            plot_title=element_text(weight="bold", ha="left"),
            panel_grid_major=element_line(color=AXIS_LINE_COLOR, size=0.5, alpha=0.5),
            panel_grid_minor=element_line(color=AXIS_LINE_COLOR, size=0.25, alpha=0.5),
        )
    )
    save_plot(pdf_handle, plot)


def render_contamination_scatter(pdf_handle, spec, tables):
    """
    Render one of the DNA contamination scatter plots.
    """
    plot_table = (
        tables[spec["source"]]
        .filter(valid_metric_expr("DNA_CONTAMINATION_SCORE") & valid_metric_expr("DNA_CONTAMINATION_P_VALUE"))
        .with_columns([
            pl.col("DNA_CONTAMINATION_SCORE").cast(pl.Float64),
            pl.col("DNA_CONTAMINATION_P_VALUE").cast(pl.Float64),
        ])
    )

    if spec.get("skip_if_empty") and plot_table.height == 0:
        return

    # Extract USL contamination score, handling NA values
    score_result = (
        tables["dna_guideline_table"]
        .filter(pl.col("SAMPLE_ID") == "USL_Guideline")
        .with_columns(
            pl.when(pl.col("DNA_CONTAMINATION_SCORE") == "NA")
            .then(None)
            .otherwise(pl.col("DNA_CONTAMINATION_SCORE"))
            .cast(pl.Float64)
            .alias("DNA_CONTAMINATION_SCORE")
        )
        .select(pl.col("DNA_CONTAMINATION_SCORE"))
        .head(1)
        .item()
    )
    usl_contamination_score = float(score_result) if score_result is not None else 5000.0
    
    # Extract USL contamination p-value, handling NA values
    pval_result = (
        tables["dna_guideline_table"]
        .filter(pl.col("SAMPLE_ID") == "USL_Guideline")
        .with_columns(
            pl.when(pl.col("DNA_CONTAMINATION_P_VALUE") == "NA")
            .then(None)
            .otherwise(pl.col("DNA_CONTAMINATION_P_VALUE"))
            .cast(pl.Float64)
            .alias("DNA_CONTAMINATION_P_VALUE")
        )
        .select(pl.col("DNA_CONTAMINATION_P_VALUE"))
        .head(1)
        .item()
    )
    usl_contamination_pval = float(pval_result) if pval_result is not None else 0.05
    max_contamination_score = max(
        5000,
        plot_table.select(pl.col("DNA_CONTAMINATION_SCORE").max()).item(),
    )

    plot = plot_contamination_scatter(
        data                    = plot_table.to_pandas(),
        color_var               = spec["color_var"],
        label_var               = spec["label_var"],
        guide_title             = spec.get("guide_title", "Run"),
        title                   = spec["title"],
        x_lab_angle             = spec.get("x_lab_angle", ANGLE_X_NAMES),
        fig_size                = spec.get("fig_size", (15, 6)),
        max_contamination_score = max_contamination_score,
        usl_contamination_score = usl_contamination_score,
        usl_contamination_pval  = usl_contamination_pval,
        color_values            = spec.get("color_values"),
    )
    save_plot(pdf_handle, plot)


def render_plot(pdf_handle, spec_name, spec, tables):
    """
    Dispatch to the appropriate renderer for one plot specification.
    """

    # Each entry in PLOT_SPECS must have a "plot_kind" field that determines which rendering function is called for that plot.
    plot_kind = spec["plot_kind"]

    if plot_kind == "bar":
        render_bar_plot(pdf_handle, spec, tables)
        return
    
    if plot_kind == "cluster_density_scatter":
        render_cluster_density_scatter(pdf_handle, spec, tables)
        return
    
    if plot_kind == "contamination_scatter":
        render_contamination_scatter(pdf_handle, spec, tables)
        return

    raise ValueError(f"Unsupported plot kind for {spec_name}: {plot_kind}")


def build_tables(joint_qc_table     : pl.DataFrame,
                 merged_tables      : pl.DataFrame,
                 highlighted_run_id : str,
                 workflow           : str = "Dragen") -> dict:
    """
    Build the logical table subsets used by the plot specifications.
    This includes separating guideline and data samples, splitting DNA and RNA metrics, and adding columns for highlighted runs and labels.
    """

    # Validate the workflow argument to ensure it is one of the supported workflows (Dragen or LocalApp).
    workflow = workflow.strip().lower()
    if workflow not in {supported.lower() for supported in SUPPORTED_WORKFLOW}:
        raise ValueError(f"Unsupported workflow: {workflow}. Expected 'Dragen' or 'LocalApp'.")

    
    # Guidelines table
    # The values are used for plotting the guideline lines and annotations in the plots
    joint_qc_guideline_table = joint_qc_table.filter(
        pl.col("RUN_ID").is_in(["LSL_Guideline", "USL_Guideline", "Internal Guideline"])
    )


    # Select workflow rows in merged table
    merged_tables = merged_tables.filter(pl.col("WORKFLOW_TYPE").str.to_lowercase() == workflow)

    guideline_table = merged_tables.filter(
        (pl.col("SAMPLE_ID") == "LSL_Guideline") | (pl.col("SAMPLE_ID") == "USL_Guideline")
    )


    internal_guideline_table = merged_tables.filter(pl.col("SAMPLE_ID") == "Internal Guideline")
    dna_guideline_table = guideline_table.filter(pl.col("DNA_CONTAMINATION_SCORE").is_not_null())
    rna_guideline_table = guideline_table.filter(pl.col("RNA_MEDIAN_CV_GENE_500X").is_not_null())
    
    data_table = merged_tables.filter(
        (pl.col("SAMPLE_ID") != "LSL_Guideline") & (pl.col("SAMPLE_ID") != "USL_Guideline")
    )
    dna_data_table = data_table.filter(pl.col("DNA_TOTAL_PF_READS").is_not_null()).clone()
    rna_data_table = data_table.filter(pl.col("RNA_TOTAL_PF_READS").is_not_null()).clone()

    dna_data_table = dna_data_table.with_columns([
        pl.when(pl.col("RUN") == highlighted_run_id).then(pl.lit("True")).otherwise(pl.lit("False")).alias("highlighted_run"),
        pl.when(pl.col("RUN") == highlighted_run_id).then(pl.col("SAMPLE_ID")).otherwise(pl.lit("")).alias("contamination_label"),
    ])

    return {
        "joint_qc_table"           : joint_qc_table,
        "joint_qc_guideline_table" : joint_qc_guideline_table,
        "merged_tables"            : merged_tables,
        "guideline_table"          : guideline_table,
        "internal_guideline_table" : internal_guideline_table,
        "dna_guideline_table"      : dna_guideline_table,
        "rna_guideline_table"      : rna_guideline_table,
        "data_table"               : data_table,
        "dna_data_table"           : dna_data_table,
        "rna_data_table"           : rna_data_table,
        "dna_sample_count"         : dna_data_table.height,
        "rna_sample_count"         : rna_data_table.height,
    }


def main():
    """Main execution function for metrics processing and plotting."""
    args = sys.argv[1:]
    
    # Configure logging to display error level in output
    logging.basicConfig(
        level=logging.INFO,
        format='[%(levelname)s] %(message)s'
    )

    # ----------------------- #
    # Reading input arguments #
    # ----------------------- #
    print("; Reading input arguments")
    parsed_args = parse_args(args)

    workflow     = parsed_args["workflow_type"]   # already lowercased by argparse
    output_pdf   = Path(parsed_args["output_pdf"])
    output_pdf.parent.mkdir(parents=True, exist_ok=True)


    # ------------------------------ #
    # Read the supplied input tables #
    # ------------------------------ #
    print("; Reading joint QC table")
    joint_qc_table = pl.read_csv(parsed_args["joint_qc_file"], separator="\t")

    print("; Reading master metrics table")
    merged_tables  = pl.read_csv(parsed_args["metrics_table"], separator="\t")

    # Validate the plot specifications to ensure that all required fields are present for each plot kind.
    validate_plot_specs(PLOT_SPECS)

    print(f"; Generating plots for workflow: {workflow}")

    # Prepare the tables in a format used for plotting:
    #   - Filtering for the specified workflow
    #   - Separating DNA and RNA samples
    #   - Adding columns for highlighted runs and labels.
    tables = build_tables(
        joint_qc_table     = joint_qc_table,
        merged_tables      = merged_tables,
        highlighted_run_id = parsed_args["highlighted_run_id"],
        workflow           = workflow,
    )

    print("; Creating the metrics plots")
    if tables["dna_sample_count"] == 0:
        print("; No DNA samples included in the plotted runs")
    if tables["rna_sample_count"] == 0:
        print("; No RNA samples included in the plotted runs")

    # -------------------- #
    # Plot metrics section #
    # -------------------- #
    # Plots are sorted by the workflow-specific index, which replaces the old plot_number field.
    # Only specs where spec[workflow]["plot"] is True are rendered for the requested workflow.
    with PdfPages(output_pdf) as pdf_handle:


        # Select only the plot specifications that are relevant for the requested workflow and have plotting enabled.
        eligible_specs = [
            (spec_name, spec)
            for spec_name, spec in PLOT_SPECS.items()
            if spec.get(workflow, {}).get("plot", False)
        ]


        # Sort the eligible plots by their workflow-specific index
        # Dtermine the order in which they will be rendered in the PDF.
        for spec_name, spec in sorted(eligible_specs, key=lambda item: item[1][workflow]["index"]):

            print(f"; Rendering plot {spec_name} (index {spec[workflow]['index']})")

            # Skip DNA/RNA-specific plots when no samples of that type are present.
            requires_samples = spec.get("requires_samples")
            if requires_samples == "dna" and tables["dna_sample_count"] == 0:
                continue
            if requires_samples == "rna" and tables["rna_sample_count"] == 0:
                continue

            # This is the main call to render each plot
            # It dispatches to the appropriate rendering function based on the "plot_kind" field in the plot specification.
            # Plot kinds:
            #   - "bar" for bar plots
            #   - "cluster_density_scatter" for the 'cluster density vs. estimated yield' scatter plot
            #   - "contamination_scatter" for the 'DNA contamination score vs. p-value' scatter plot
            render_plot(pdf_handle, spec_name, spec, tables)

    print(f"; Plots saved to PDF: {output_pdf}")


if __name__ == "__main__":
    main()