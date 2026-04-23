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
import polars as pl
from plotnine import *
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

from plots import plot_tsoppy_barplot


# This can be moved to a config file
VERSION_STRING = "0.3:22-03-10"


def main():
    """
    Main execution function for metrics processing and plotting.
    """
    
    # ---------------------------- #
    # Parse command-line arguments #
    # ---------------------------- #
    args = sys.argv[1:]  # Skip script name
    print("; Reading input arguments")
    
    run_count            = int(args[0])
    joint_qc_file        = args[1]
    output_pdf           = args[2]
    merged_metrics_file  = args[3]
    highlighted_run_id   = args[4]
    create_plots         = args[5]
    

    # ---------------------------------- #
    # Check input arguments for validity #
    # ---------------------------------- #

    # Exit if zero runs are submitted
    if run_count < 1:
        sys.exit("Zero runs provided as input? Exiting.")
    

    # --------------------------- #
    # Load and process input data #
    # --------------------------- #

    # Load joint QC table
    joint_qc_table = pl.read_csv(joint_qc_file, separator="\t")
    
    # Load first run's metrics table
    merged_tables = pl.read_csv(args[6], separator="\t", infer_schema_length=0)
    merged_tables = merged_tables.with_columns(
        pl.lit(args[6 + run_count]).alias('RUN')
    )
    
    print("; Merging input tables")
    
    # Merge additional runs if present
    if run_count > 1:
        for run_number in range(2, run_count + 1):
            run_table = pl.read_csv(args[5 + run_number], separator="\t", infer_schema_length=0)
            run_table = run_table.with_columns(
                pl.lit(args[5 + run_count + run_number]).alias('RUN')
            )
            print(f";   Iteration {run_number}")
            merged_tables = pl.concat(
                [merged_tables, run_table]
            ).unique()
    
    # Write merged metrics to output file
    print("; Writing merged input tables into a master metrics file.")
    merged_tables.unique().write_csv(merged_metrics_file, separator="\t")
    
    if create_plots == "False":
        print("  Skipping metrics plotting. Exiting.")
        return
    
    # Separate guideline values from actual data
    # Keep DNA guideline values, RNA guideline values,
    # RNA data and DNA data separate from each other
    # (they are to be used differently in the plotting)
    guideline_table = merged_tables.filter((pl.col('SAMPLE_ID') == 'LSL_Guideline') |
                                           (pl.col('SAMPLE_ID') == 'USL_Guideline') 
                                          )
    
    # Both LSL and USL guidelines are specified
    # for the CONTAMINATION_SCORE DNA metric
    dna_guideline_table = guideline_table.filter(
        pl.col('DNA_CONTAMINATION_SCORE').is_not_null()
    )
    
    # Both LSL and USL guidelines are specified
    # for the MEDIAN_CV_GENE_500X RNA metric
    rna_guideline_table = guideline_table.filter(
        pl.col('RNA_MEDIAN_CV_GENE_500X').is_not_null()
    )
    
    # Filter out guideline rows from data table
    data_table = merged_tables.filter(
        (pl.col('SAMPLE_ID') != 'LSL_Guideline') &
        (pl.col('SAMPLE_ID') != 'USL_Guideline')
    )
    
    # Each sample with specified DNA TOTAL_PF_READS count
    # is included in the DNA data table
    dna_data_table = data_table.filter(
        pl.col('DNA_TOTAL_PF_READS').is_not_null()
    ).clone()
    
    # Each sample with specified RNA TOTAL_PF_READS count
    # is included in the RNA data table
    rna_data_table = data_table.filter(
        pl.col('RNA_TOTAL_PF_READS').is_not_null()
    ).clone()
    
    dna_sample_count = dna_data_table.height
    rna_sample_count = rna_data_table.height
    
    # Samples from the highlighted run will be labeled
    # on the contamination scatter plot
    dna_data_table = dna_data_table.with_columns([
        pl.when(pl.col('RUN') == highlighted_run_id)
          .then(pl.lit('True'))
          .otherwise(pl.lit('False'))
          .alias('highlighted_run'),
        pl.when(pl.col('RUN') == highlighted_run_id)
          .then(pl.col('SAMPLE_ID'))
          .otherwise(pl.lit(''))
          .alias('contamination_label')
    ])
    
    print("; Creating the metrics plots")
    
    # Open PDF file for multi-page output
    pdf = PdfPages(output_pdf)


    # ------------- #
    # Plot settings #
    # ------------- #
    angle_x_names = 35


    # --------------------------- #
    # Sequencer-generated metrics #
    # --------------------------- #
    
    # ------------------------------- #
    # Plot 1: Clusters passing filter #
    # ------------------------------- #
    # Filter for non-null CLUSTERS_PASSING_FILTER values and cast to float
    filtered_qc = (joint_qc_table.filter(pl.col('CLUSTERS_PASSING_FILTER').is_not_null())
                                 .with_columns(pl.col('CLUSTERS_PASSING_FILTER').cast(pl.Float64))
                  )
    

    if filtered_qc.height > 0:
        # Convert to pandas and ensure numeric types
        plot_data = filtered_qc.to_pandas()
        

        plot = plot_tsoppy_barplot(data        = plot_data,
                                   x_var       = 'RUN_ID',
                                   y_var       = 'CLUSTERS_PASSING_FILTER',
                                   fill_var    = 'RUN_ID',
                                   guide_title = 'Run',
                                   x_lab       = 'Run ID',
                                   y_lab       = 'Clusters passing filter (percent)',
                                   cart_ylim   = (0, 104),
                                   title       = '[Sequencer flowcell run metric] Clusters passing filter (based on RunCompletionStatus.xml files)\n(only displaying runs with non-NA values)',
                                   x_lab_angle = angle_x_names,
                                   fig_size    = (15, 6),
                                   alpha_value = 0.8)
        
        fig = plot.draw()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)


    # ----------------------------------------------- #
    # Plot 2: Cluster density against estimated yield #
    # ----------------------------------------------- #
    # Filter for non-null ESTIMATED_YIELD values and cast numeric columns
    filtered_qc2 = (joint_qc_table.filter(pl.col('ESTIMATED_YIELD').is_not_null())
                                   .with_columns([
                                       pl.col('CLUSTER_DENSITY').cast(pl.Float64),
                                       pl.col('ESTIMATED_YIELD').cast(pl.Float64),
                                       (pl.col('RUN_NUMBER').cast(pl.Utf8) + ' - ' + pl.col('RUN_ID')).alias('RUN_LABEL')
                                   ])
                   )
    
    if filtered_qc2.height > 0:
        plot_data2 = filtered_qc2.to_pandas()
        
        plot = (
            ggplot(plot_data2, aes(x='CLUSTER_DENSITY', y='ESTIMATED_YIELD')) +
            geom_point(aes(color='RUN_LABEL'), alpha=0.8, size=4.5) +
            geom_text(aes(label='RUN_NUMBER'), alpha=0.8, size=10, nudge_y=-5) +
            guides(color=guide_legend('Run')) +
            coord_cartesian(xlim=(0, 400), ylim=(0, 150)) +
            xlab('Cluster density') +
            ylab('Estimated yield') +
            ggtitle(
                '[Sequencer flowcell run metric] Cluster density against estimated yield '
                '(based on RunCompletionStatus.xml files)\n'
                '(only displaying runs with non-NA values)'
            ) +
            scale_x_continuous(breaks=[i * 50 for i in range(9)], minor_breaks=[i * 25 for i in range(17)]) +
            scale_y_continuous(breaks=[i * 10 for i in range(16)], minor_breaks=[i * 5 for i in range(31)]) +
            theme_minimal() +
            theme(
                axis_text_x=element_text(angle=angle_x_names, vjust=0.5, hjust=1, size=7),
                figure_size=(15, 6),
                axis_line=element_line(color="#4d4d4d"),
                plot_title=element_text(weight='bold', ha = 'left'),
                panel_grid_major=element_line(color='gray', size=0.5),
                panel_grid_minor=element_line(color='lightgray', size=0.25)
            )
        )
        
        fig = plot.draw()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)


    # --------------------------- #
    # LocalApp core run metrics   #
    # --------------------------- #
    
    # ------------------------------ #
    # Plot 3: Percentage of PF reads #
    # ------------------------------ #
    # Cast numeric column to float
    plot_data3 = (joint_qc_table.with_columns(
                      pl.col('PCT_PF_READS').cast(pl.Float64)
                  ).to_pandas())

    plot = plot_tsoppy_barplot(data        = plot_data3,
                               x_var       = 'RUN_ID',
                               y_var       = 'PCT_PF_READS',
                               fill_var    = 'RUN_ID',
                               guide_title = 'Run',
                               x_lab       = 'Run ID',
                               y_lab       = 'Percentage of PF reads',
                               cart_ylim   = (0, 104),
                               title       = '[Sequencer flowcell run metric] Percentage of PF reads (based on RunCompletionStatus.xml files)\n(only displaying runs with non-NA values)',
                               x_lab_angle = angle_x_names,
                               fig_size    = (15, 6),
                               alpha_value = 0.8,
                               hline_y     = 80,
                               hline_alpha = 0.3,
                               hline_color = "red",
                               hline_label = "LSL_Guideline: 80")

    fig = plot.draw()
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)
    

    # ------------------------------------ #
    # Plot 4: Percentage of Q30 R1 reads   #
    # ------------------------------------ #
    # Cast numeric column to float
    plot_data4 = (joint_qc_table.with_columns(
                      pl.col('PCT_Q30_R1').cast(pl.Float64)
                  ).to_pandas())
    
    plot = plot_tsoppy_barplot(data        = plot_data4,
                               x_var       = 'RUN_ID',
                               y_var       = 'PCT_Q30_R1',
                               fill_var    = 'RUN_ID',
                               guide_title = 'Run',
                               x_lab       = 'Run ID',
                               y_lab       = 'Percentage of Q30 R1 reads',
                               cart_ylim   = (0, 104),
                               title       = '[LocalApp core run metric] Percentage of R1 reads with quality >= 30',
                               x_lab_angle = angle_x_names,
                               fig_size    = (15, 6),
                               alpha_value = 0.8,
                               hline_y     = 80,
                               hline_alpha = 0.3,
                               hline_color = "red",
                               hline_label = "LSL_Guideline: 80")

    fig = plot.draw()
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)
    

    # ------------------------------------ #
    # Plot 5: Percentage of Q30 R2 reads   #
    # ------------------------------------ #
    # Cast numeric column to float
    plot_data5 = (joint_qc_table.with_columns(
                      pl.col('PCT_Q30_R2').cast(pl.Float64)
                  ).to_pandas())
    
    plot = plot_tsoppy_barplot(data        = plot_data5,
                               x_var       = 'RUN_ID',
                               y_var       = 'PCT_Q30_R2',
                               fill_var    = 'RUN_ID',
                               guide_title = 'Run',
                               x_lab       = 'Run ID',
                               y_lab       = 'Percentage of Q30 R2 reads',
                               cart_ylim   = (0, 104),
                               title       = '[LocalApp core run metric] Percentage of R2 reads with quality >= 30',
                               x_lab_angle = angle_x_names,
                               fig_size    = (15, 6),
                               alpha_value = 0.8,
                               hline_y     = 80,
                               hline_alpha = 0.3,
                               hline_color = "red",
                               hline_label = "LSL_Guideline: 80")

    fig = plot.draw()
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)
    
    # ---------------------------------------- #
    # Close PDF file after all plots are saved #
    # ---------------------------------------- # 
    pdf.close()
    print("; Plots saved to PDF")


if __name__ == "__main__":
    main()
