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

from plots import plot_tsoppy_barplot, plot_contamination_scatter, AXIS_LINE_COLOR


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
                axis_line=element_line(color = AXIS_LINE_COLOR),
                plot_title=element_text(weight='bold', ha = 'left'),
                panel_grid_major=element_line(color = AXIS_LINE_COLOR, size=0.5, alpha=0.5),
                panel_grid_minor=element_line(color = AXIS_LINE_COLOR, size=0.25, alpha=0.5)
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
    

    # ----------------------------- #
    # LocalApp core DNA metrics     #
    # ----------------------------- #
    
    if dna_sample_count > 0:
        
        # ------------------------------------------------------------- #
        # Plot 6: Contamination assessment by TSO 500 LocalApp (by Run) #
        # ------------------------------------------------------------- #
        
        # Extract guideline values - cast to Float64 first, handling "NA" strings
        usl_contamination_score = float(
            dna_guideline_table.filter(pl.col('SAMPLE_ID') == 'USL_Guideline')
                               .select(pl.col('DNA_CONTAMINATION_SCORE').cast(pl.Float64))
                               .item()
        )
        
        usl_contamination_pval = float(
            dna_guideline_table.filter(pl.col('SAMPLE_ID') == 'USL_Guideline')
                               .select(pl.col('DNA_CONTAMINATION_P_VALUE').cast(pl.Float64))
                               .item()
        )
        
        # Calculate max contamination score - filter out NA values first
        max_contamination_score_from_data = (
            dna_data_table.filter(pl.col('DNA_CONTAMINATION_SCORE') != 'NA')
                         .select(pl.col('DNA_CONTAMINATION_SCORE').cast(pl.Float64))
                         .max()
                         .item()
        )
        max_contamination_score = max(5000, max_contamination_score_from_data)
        
        # Prepare plot data - filter out NA values and cast to Float64
        plot_data_contamination = (dna_data_table
                                  .filter((pl.col('DNA_CONTAMINATION_SCORE') != 'NA') & 
                                          (pl.col('DNA_CONTAMINATION_P_VALUE') != 'NA'))
                                  .with_columns([
                                      pl.col('DNA_CONTAMINATION_SCORE').cast(pl.Float64),
                                      pl.col('DNA_CONTAMINATION_P_VALUE').cast(pl.Float64)
                                  ])
                                  .to_pandas())
        
        plot = plot_contamination_scatter(data                    = plot_data_contamination,
                                          color_var               = 'RUN',
                                          label_var               = None,
                                          guide_title             = 'Run',
                                          title                   = ('[LocalApp core DNA metric] Contamination assessment by TSO 500 LocalApp\n'
                                                                     '(samples located within the red box are deemed contaminated)'),
                                          x_lab_angle             = angle_x_names,
                                          fig_size                = (15, 6),
                                          max_contamination_score = max_contamination_score,
                                          usl_contamination_score = usl_contamination_score,
                                          usl_contamination_pval  = usl_contamination_pval
        )
        
        fig = plot.draw()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)
        
        
        # ---------------------------------------------------------------------- #
        # Plot 7: Contamination assessment by TSO 500 LocalApp (highlighted run) #
        # ---------------------------------------------------------------------- #
        
        plot = plot_contamination_scatter(data                    = plot_data_contamination,
                                          color_var               = 'highlighted_run',
                                          label_var               = 'contamination_label',
                                          guide_title             = 'Latest run',
                                          title                   = ('[LocalApp core DNA metric] Contamination assessment by TSO 500 LocalApp\n'
                                                                     '(samples located within the red box are deemed contaminated)'),
                                          x_lab_angle             = angle_x_names,
                                          fig_size                = (15, 6),
                                          max_contamination_score = max_contamination_score,
                                          usl_contamination_score = usl_contamination_score,
                                          usl_contamination_pval  = usl_contamination_pval
                                         )
        
        fig = plot.draw()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)
        
        
        # ----------------------------------------- #
        # Plot 8: Median exon coverage (DNA)        #
        # ----------------------------------------- #
        
        # Extract guideline value
        lsl_dna_median_exon_coverage = int(
            dna_guideline_table.filter(pl.col('SAMPLE_ID') == 'LSL_Guideline')
                               .select(pl.col('DNA_MEDIAN_EXON_COVERAGE').cast(pl.Int32))
                               .item()
        )
        
        # Prepare plot data
        plot_data8 = (dna_data_table
                     .filter(pl.col('DNA_MEDIAN_EXON_COVERAGE') != 'NA')
                     .with_columns(pl.col('DNA_MEDIAN_EXON_COVERAGE').cast(pl.Float64))
                     .to_pandas())
        
        plot = plot_tsoppy_barplot(data        = plot_data8,
                                   x_var       = 'SAMPLE_ID',
                                   y_var       = 'DNA_MEDIAN_EXON_COVERAGE',
                                   fill_var    = 'RUN',
                                   guide_title = 'Run',
                                   x_lab       = 'Sample ID',
                                   y_lab       = 'Median exon coverage',
                                   title       = '[LocalApp core DNA metric] Median exon coverage',
                                   fig_size    = (15, 6),
                                   alpha_value = 0.8,
                                   hline_y     = lsl_dna_median_exon_coverage,
                                   hline_label = f"LSL_Guideline: {lsl_dna_median_exon_coverage}"
        )
        
        fig = plot.draw()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)


        # --------------------------------------------------------------- #
        # Plot 9: Percentage of exons with coverage >= 50X (DNA)         #
        # --------------------------------------------------------------- #

        # Extract guideline value
        lsl_dna_pct_exon_50x = int(
            dna_guideline_table.filter(pl.col('SAMPLE_ID') == 'LSL_Guideline')
                               .select(pl.col('DNA_PCT_EXON_50X').cast(pl.Float64))
                               .item()
        )

        # Prepare plot data
        plot_data9 = (dna_data_table
                     .filter(pl.col('DNA_PCT_EXON_50X') != 'NA')
                     .with_columns(pl.col('DNA_PCT_EXON_50X').cast(pl.Float64))
                     .to_pandas())

        plot = plot_tsoppy_barplot(data         = plot_data9,
                                   x_var        = 'SAMPLE_ID',
                                   y_var        = 'DNA_PCT_EXON_50X',
                                   fill_var     = 'RUN',
                                   guide_title  = 'Run',
                                   x_lab        = 'Sample ID',
                                   y_lab        = 'Percentage of exons with coverage >= 50X',
                                   cart_ylim    = (0, 105),
                                   title        = '[LocalApp core DNA metric] Percentage of exons with coverage >= 50X',
                                   fig_size     = (15, 6),
                                   alpha_value  = 0.8,
                                   hline_y      = lsl_dna_pct_exon_50x,
                                   hline_label  = f"LSL_Guideline: {lsl_dna_pct_exon_50x}",
                                   ann_y_offset = -5
        )

        fig = plot.draw()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)


        # ----------------------------------------- #
        # Plot 10: Median insert size (DNA)        #
        # ----------------------------------------- #

        # Extract guideline value
        lsl_median_insert_size = int(
            dna_guideline_table.filter(pl.col('SAMPLE_ID') == 'LSL_Guideline')
                               .select(pl.col('DNA_MEDIAN_INSERT_SIZE').cast(pl.Int32))
                               .item()
        )

        # Prepare plot data
        plot_data10 = (dna_data_table
                      .filter(pl.col('DNA_MEDIAN_INSERT_SIZE') != 'NA')
                      .with_columns(pl.col('DNA_MEDIAN_INSERT_SIZE').cast(pl.Float64))
                      .to_pandas())

        plot = plot_tsoppy_barplot(data         = plot_data10,
                                   x_var        = 'SAMPLE_ID',
                                   y_var        = 'DNA_MEDIAN_INSERT_SIZE',
                                   fill_var     = 'RUN',
                                   guide_title  = 'Run',
                                   x_lab        = 'Sample ID',
                                   y_lab        = 'Median insert size',
                                   title        = '[LocalApp core DNA metric] Median insert size',
                                   fig_size     = (15, 6),
                                   alpha_value  = 0.8,
                                   hline_y      = lsl_median_insert_size,
                                   hline_label  = f"LSL_Guideline: {lsl_median_insert_size}",
                                   ann_y_offset = 10
        )

        fig = plot.draw()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)


        # -------------------------------------------------------------- #
        # Plot 11: Coverage MAD (DNA)                                   #
        # -------------------------------------------------------------- #

        # Extract guideline value
        usl_median_coverage_mad = float(
            dna_guideline_table.filter(pl.col('SAMPLE_ID') == 'USL_Guideline')
                               .select(pl.col('DNA_COVERAGE_MAD').cast(pl.Float64))
                               .item()
        )

        # Prepare plot data
        plot_data11 = (dna_data_table
                      .filter(pl.col('DNA_COVERAGE_MAD') != 'NA')
                      .with_columns(pl.col('DNA_COVERAGE_MAD').cast(pl.Float64))
                      .to_pandas())

        plot = plot_tsoppy_barplot(data         = plot_data11,
                                   x_var        = 'SAMPLE_ID',
                                   y_var        = 'DNA_COVERAGE_MAD',
                                   fill_var     = 'RUN',
                                   guide_title  = 'Run',
                                   x_lab        = 'Sample ID',
                                   y_lab        = 'Coverage MAD (Median Absolute Deviation)',
                                   title        = '[LocalApp core DNA metric] Coverage MAD\n(median normalized deviation across all regions used for CNV calling)',
                                   fig_size     = (15, 6),
                                   alpha_value  = 0.8,
                                   hline_y      = usl_median_coverage_mad,
                                   hline_label  = f"USL_Guideline: {usl_median_coverage_mad}",
                                   ann_y_offset = 0
        )

        fig = plot.draw()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)


        # ---------------------------------------------------------------- #
        # Plot 12: Median raw bin count per CNV target (DNA)             #
        # ---------------------------------------------------------------- #

        # Extract guideline value
        lsl_median_bin_count = int(
            dna_guideline_table.filter(pl.col('SAMPLE_ID') == 'LSL_Guideline')
                               .select(pl.col('DNA_MEDIAN_BIN_COUNT_CNV_TARGET').cast(pl.Int32))
                               .item()
        )

        # Prepare plot data
        plot_data12 = (dna_data_table
                      .filter(pl.col('DNA_MEDIAN_BIN_COUNT_CNV_TARGET') != 'NA')
                      .with_columns(pl.col('DNA_MEDIAN_BIN_COUNT_CNV_TARGET').cast(pl.Float64))
                      .to_pandas())

        plot = plot_tsoppy_barplot(data         = plot_data12,
                                   x_var        = 'SAMPLE_ID',
                                   y_var        = 'DNA_MEDIAN_BIN_COUNT_CNV_TARGET',
                                   fill_var     = 'RUN',
                                   guide_title  = 'Run',
                                   x_lab        = 'Sample ID',
                                   y_lab        = 'Median raw bin count per CNV target',
                                   title        = '[LocalApp core DNA metric] Median raw bin count per CNV target',
                                   fig_size     = (15, 6),
                                   alpha_value  = 0.8,
                                   hline_y      = lsl_median_bin_count,
                                   hline_label  = f"LSL_Guideline: {lsl_median_bin_count}",
                                   ann_y_offset = 1
        )

        fig = plot.draw()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)


        # ---------------------------------------------------------------- #
        # Plot 13: Number of usable MSI sites (DNA)                       #
        # ---------------------------------------------------------------- #

        # Extract guideline value
        lsl_usable_msi_sites = int(
            dna_guideline_table.filter(pl.col('SAMPLE_ID') == 'LSL_Guideline')
                               .select(pl.col('DNA_USABLE_MSI_SITES').cast(pl.Int32))
                               .item()
        )

        # Prepare plot data
        plot_data13 = (dna_data_table
                      .filter(pl.col('DNA_USABLE_MSI_SITES') != 'NA')
                      .with_columns(pl.col('DNA_USABLE_MSI_SITES').cast(pl.Float64))
                      .to_pandas())

        plot = plot_tsoppy_barplot(data         = plot_data13,
                                   x_var        = 'SAMPLE_ID',
                                   y_var        = 'DNA_USABLE_MSI_SITES',
                                   fill_var     = 'RUN',
                                   guide_title  = 'Run',
                                   x_lab        = 'Sample ID',
                                   y_lab        = 'Usable MSI sites',
                                   title        = '[LocalApp core DNA metric] Number of usable MSI sites',
                                   fig_size     = (15, 6),
                                   alpha_value  = 0.8,
                                   hline_y      = lsl_usable_msi_sites,
                                   hline_label  = f"LSL_Guideline: {lsl_usable_msi_sites}",
                                   ann_y_offset = 0
        )

        fig = plot.draw()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)
    
    else:
        print("; No DNA samples included in the plotted runs")
    
    # ---------------------------------------- #
    # Close PDF file after all plots are saved #
    # ---------------------------------------- # 
    pdf.close()
    print("; Plots saved to PDF")


if __name__ == "__main__":
    main()
