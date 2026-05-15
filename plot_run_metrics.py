#!/usr/bin/env python3
"""
TSO500 Run Metrics Plotting Script

This script processes and visualizes quality metrics from TSO500 sequencing runs.
Converts metrics data from multiple runs into consolidated outputs and generates
comprehensive QC plots.

Command-line parameters:
    1: N (number of runs)
    2: joint run QC file
    3: output PDF file with plotting results
    4: output tsv file with merged metrics data
    5: label of the run that should be highlighted in the plots
    6: create_plots (True/False flag)
    7 to (6+N): transposed metrics matrices for the individual runs
    (7+N) to (6+2N): labels for the individual runs
"""

import sys

import matplotlib.pyplot as plt
import polars as pl
from matplotlib.backends.backend_pdf import PdfPages
from plotnine import *

from plots import AXIS_LINE_COLOR, plot_contamination_scatter, plot_tsoppy_barplot
from plot_specs import PLOT_SPECS


VERSION_STRING = "0.3:22-03-10"
ANGLE_X_NAMES = 35



def parse_args(argv):
    """
    Parse command-line arguments while preserving the existing CLI contract.
    """

    run_count = int(argv[0])
    return {
        "run_count": run_count,
        "joint_qc_file": argv[1],
        "output_pdf": argv[2],
        "merged_metrics_file": argv[3],
        "highlighted_run_id": argv[4],
        "create_plots": argv[5],
        "metric_files": argv[6:6 + run_count],
        "run_labels": argv[6 + run_count:6 + (2 * run_count)],
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
        value = filtered_table.select(expression).item()
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

def merge_metric_tables(parsed_args):
    """
    Read all per-run metric tables, tag them with RUN, and merge them.
    """
    run_tables = []
    print("; Merging input tables")

    for index, (metric_file, run_label) in enumerate(zip(parsed_args["metric_files"], parsed_args["run_labels"]), start=1):
        run_table = (pl.read_csv(metric_file, separator="\t", infer_schema_length=0)
                       .with_columns(pl.lit(run_label).alias("RUN"))
                    )
        if index > 1:
            print(f";   Iteration {index}")
        run_tables.append(run_table)

    return pl.concat(run_tables).unique()


def build_tables(joint_qc_table, merged_tables, highlighted_run_id):
    """
    Build the logical table subsets used by the plot specifications.
    This includes separating guideline and data samples, splitting DNA and RNA metrics, and adding columns for highlighted runs and labels.
    """

    # Guidelines table
    # The values are used for plotting the guideline lines and annotations in the plots
    joint_qc_guideline_table = joint_qc_table.filter(
        pl.col("RUN_ID").is_in(["LSL_Guideline", "USL_Guideline", "Internal Guideline"])
    )


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
        "joint_qc_table": joint_qc_table,
        "joint_qc_guideline_table": joint_qc_guideline_table,
        "merged_tables": merged_tables,
        "guideline_table": guideline_table,
        "internal_guideline_table": internal_guideline_table,
        "dna_guideline_table": dna_guideline_table,
        "rna_guideline_table": rna_guideline_table,
        "data_table": data_table,
        "dna_data_table": dna_data_table,
        "rna_data_table": rna_data_table,
        "dna_sample_count": dna_data_table.height,
        "rna_sample_count": rna_data_table.height,
    }



def compute_cart_ylim(spec, plot_data):
    """
    Resolve static or data-driven y-axis limits for a bar plot.
     - If "cart_ylim_dynamic" is not specified in the plot specification, use the static limits defined in "cart_ylim".
     - If "cart_ylim_dynamic" is specified with mode "max_plus", calculate the maximum value of the specified column in the plot data and add the defined offset to determine the upper limit
    """
    dynamic_ylim = spec.get("cart_ylim_dynamic")
    if dynamic_ylim is None:
        return spec.get("cart_ylim")

    if dynamic_ylim["mode"] == "max_plus":
        max_value = plot_data[dynamic_ylim["column"]].max()
        return (dynamic_ylim.get("lower", 0), max_value + dynamic_ylim.get("offset", 0))

    raise ValueError(f"Unsupported dynamic y-limit mode: {dynamic_ylim['mode']}")


def prepare_bar_plot_data(table, spec):
    """Apply per-plot filtering and column derivation before converting to pandas."""
    filters = []
    for column_name in spec.get("na_filter_columns", []):
        filters.append(valid_metric_expr(column_name))
    for filter_spec in spec.get("filters", []):
        filters.append(build_filter_expression(filter_spec))

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
            axis_text_x=element_text(angle=spec.get("x_lab_angle", ANGLE_X_NAMES), vjust=0.5, hjust=1, size=7),
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

    usl_contamination_score = float(
        tables["dna_guideline_table"]
        .filter(pl.col("SAMPLE_ID") == "USL_Guideline")
        .select(pl.col("DNA_CONTAMINATION_SCORE").cast(pl.Float64))
        .item()
    )
    usl_contamination_pval = float(
        tables["dna_guideline_table"]
        .filter(pl.col("SAMPLE_ID") == "USL_Guideline")
        .select(pl.col("DNA_CONTAMINATION_P_VALUE").cast(pl.Float64))
        .item()
    )
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



def main():
    """Main execution function for metrics processing and plotting."""
    args = sys.argv[1:]

    # ---------------------- #
    # Reading inpu arguments #
    # ---------------------- #
    print("; Reading input arguments")
    parsed_args = parse_args(args)
    if parsed_args["run_count"] < 1:
        sys.exit("Zero runs provided as input? Exiting.")


    # ------------------------- #
    # Join input metrics tables #
    # ------------------------- #
    joint_qc_table = pl.read_csv(parsed_args["joint_qc_file"], separator="\t")
    merged_tables = merge_metric_tables(parsed_args)

    print("; Writing merged input tables into a master metrics file.")
    merged_tables.write_csv(parsed_args["merged_metrics_file"], separator="\t")

    # ------------------------------------------- #
    # Master metrics is created in the first call #
    # ------------------------------------------- #
    if parsed_args["create_plots"] == "False":
        print("  Skipping metrics plotting. Exiting.")
        return


    tables = build_tables(joint_qc_table, merged_tables, parsed_args["highlighted_run_id"])

    print("; Creating the metrics plots")
    if tables["dna_sample_count"] == 0:
        print("; No DNA samples included in the plotted runs")
    if tables["rna_sample_count"] == 0:
        print("; No RNA samples included in the plotted runs")

    # -------------------- #
    # Plot metrics section #
    # -------------------- #
    # The plots are created in the order of the plot specifications, which is determined by the "plot_number" field in the specifications.
    # This allows to control the order of the plots in the output PDF by adjusting the plot numbers in the plot_specs.py file.
    # The plotting functions are designed to be flexible and driven by the specifications, so that new plots can be added by simply adding new entries to the plot_specs.py file without needing to change the plotting code.
    with PdfPages(parsed_args["output_pdf"]) as pdf_handle:

        # spec_name and spec are the key and value of each entry in the PLOT_SPECS dictionary defined in plot_specs.py
        for spec_name, spec in sorted(PLOT_SPECS.items(), key=lambda item: item[1]["plot_number"]):

            # Check if the plot requires DNA or RNA samples and skip if the corresponding sample count is zero
            # This allows to automatically skip DNA-specific plots if no DNA samples are included in the runs, and similarly for RNA-specific plots, without needing to adjust the plot specifications.
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

    print("; Plots saved to PDF")


if __name__ == "__main__":
    main()
