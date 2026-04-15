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

from plotnine import ggplot


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
    
    print(dna_data_table)
    print("; Creating the metrics plots")
    # Plotting functionality to be implemented
    # (continues from line 92 onwards in the R script)


if __name__ == "__main__":
    main()
