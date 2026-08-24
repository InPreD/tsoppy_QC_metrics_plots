import polars as pl


# PLOT_SPECS structure
# ====================
#
# Each top-level entry defines one plot and is keyed by the plotted metric name.
# Only include fields that are actually consumed by plot_run_metrics.py.
#
# 1. Common entry fields
#    plot_number
#        Integer ordering in the output PDF.
#    plot_kind
#        Renderer selector: "bar", "cluster_density_scatter", or "contamination_scatter".
#    source
#        Key of the prepared table returned by build_tables().
#    requires_samples
#        Optional sample gate: "dna", "rna", or None.
#    title
#        Plot title.
#    skip_if_empty
#        Optional boolean. Skip rendering when the filtered table is empty.
#
# 2. Bar plot fields
#    x_var, y_var, fill_var
#        Column names used for plotting.
#    x_lab, y_lab
#        Axis labels.
#    value_spec
#        Defines how to build the plotted y-value column.
#    na_filter_columns
#        Optional list of columns that must be non-null and not "NA".
#    filters
#        Optional list of extra row filters.
#    cart_ylim
#        Optional fixed y-axis limits tuple.
#    cart_ylim_dynamic
#        Optional dynamic y-axis rule.
#    y_tick_step
#        Optional fixed y-axis tick spacing.
#    guide_title, x_lab_angle, fig_size, alpha_value
#        Optional plot styling overrides.
#    overlay_labels
#        Optional text annotations drawn on the bars.
#    guideline
#        Optional horizontal reference line definition.
#
# 3. guideline fields
#    table
#        Key of the guideline table returned by build_tables().
#    sample_id
#        Identifier value used to select the guideline row.
#    id_column
#        Optional row-selector column name. Defaults to SAMPLE_ID.
#    value_spec
#        Defines how to read or derive the guideline value.
#    python_cast
#        Optional Python cast applied after extraction.
#    label_prefix
#        Prefix used to build the annotation label.
#    label
#        Optional explicit label override.
#    alpha, color, ann_y_offset
#        Optional line and annotation styling.
#
# 4. value_spec fields
#    operation = "cast"
#        Requires: column
#        Optional: dtype
#    operation = "divide"
#        Requires: column, divisor
#        Optional: dtype
#    operation = "ratio"
#        Requires: numerator, denominator
#        Optional: dtype, numerator_divisor
#
# 5. filters entries
#    Each filter entry must define:
#        column
#    And exactly one of:
#        contains, equals, not_equals


PLOT_SPECS = {

    # ========================================= #
    # LocalApp-specific metrics (not in Dragen) #
    # ========================================= #

    "CLUSTERS_PASSING_FILTER": {

        "localapp": {
            "plot": True,
            "index": 1
        },
        "dragen": {
            "plot": False,
            "index": 0
        },

        "plot_kind"         : "bar",
        "source"            : "joint_qc_table",
        "requires_samples"  : None,
        "x_var"             : "RUN_ID",
        "y_var"             : "CLUSTERS_PASSING_FILTER",
        "fill_var"          : "RUN_ID",
        "x_lab"             : "Run ID",
        "y_lab"             : "Clusters passing filter (percent)",
        "title"             : "[Sequencer flowcell run metric] Clusters passing filter (based on RunCompletionStatus.xml files)\n(only displaying runs with non-NA values)",
        "value_spec"        : {"operation": "cast", "column": "CLUSTERS_PASSING_FILTER", "dtype": pl.Float64},
        "na_filter_columns" : ["CLUSTERS_PASSING_FILTER"],
        "cart_ylim"         : (0, 104),
        "skip_if_empty"     : True,
    },


    "ESTIMATED_YIELD": {

        "localapp": {
            "plot": True,
            "index": 2
        },
        "dragen": {
            "plot": False,
            "index": 0
        },

        "plot_kind"         : "cluster_density_scatter",
        "source"            : "joint_qc_table",
        "requires_samples"  : None,
        "title"             : "[Sequencer flowcell run metric] Cluster density against estimated yield (based on RunCompletionStatus.xml files)\n(only displaying runs with non-NA values)",
        "skip_if_empty"     : True,
    },


    "PCT_PF_READS": {

        "localapp": {
            "plot": True,
            "index": 3
        },
        "dragen": {
            "plot": False,
            "index": 0
        },

        "plot_kind"         : "bar",
        "source"            : "joint_qc_table",
        "requires_samples"  : None,
        "x_var"             : "RUN_ID",
        "y_var"             : "PCT_PF_READS",
        "fill_var"          : "RUN_ID",
        "x_lab"             : "Run ID",
        "y_lab"             : "Percentage of PF reads",
        "title"             : "[Sequencer flowcell run metric] Percentage of PF reads (based on RunCompletionStatus.xml files)\n(only displaying runs with non-NA values)",
        "value_spec"        : {"operation": "cast", "column": "PCT_PF_READS", "dtype": pl.Float64},
        "cart_ylim"         : (0, 104),
        "guideline"         : {
            "table"        : "joint_qc_guideline_table",
            "id_column"    : "RUN_ID",
            "sample_id"    : "LSL_Guideline",
            "value_spec"   : {"operation": "cast", "column": "PCT_PF_READS", "dtype": pl.Float64},
            "python_cast"  : float,
            "label_prefix" : "LSL_Guideline",
            "alpha"        : 0.3,
            "color"        : "red",
        },
    },


    "PCT_Q30_R1": {

        "localapp": {
            "plot": True,
            "index": 4
        },
        "dragen": {
            "plot": False,
            "index": 0
        },

        "plot_kind"         : "bar",
        "source"            : "joint_qc_table",
        "requires_samples"  : None,
        "x_var"             : "RUN_ID",
        "y_var"             : "PCT_Q30_R1",
        "fill_var"          : "RUN_ID",
        "x_lab"             : "Run ID",
        "y_lab"             : "Percentage of Q30 R1 reads",
        "title"             : "[LocalApp core run metric] Percentage of R1 reads with quality >= 30",
        "value_spec"        : {"operation": "cast", "column": "PCT_Q30_R1", "dtype": pl.Float64},
        "cart_ylim"         : (0, 104),
        "guideline"         : {
            "table"        : "joint_qc_guideline_table",
            "id_column"    : "RUN_ID",
            "sample_id"    : "LSL_Guideline",
            "value_spec"   : {"operation": "cast", "column": "PCT_Q30_R1", "dtype": pl.Float64},
            "python_cast"  : float,
            "label_prefix" : "LSL_Guideline",
            "alpha"        : 0.3,
            "color"        : "red",
        },
    },


    "PCT_Q30_R2": {

        "localapp": {
            "plot": True,
            "index": 5
        },
        "dragen": {
            "plot": False,
            "index": 0
        },

        "plot_kind"         : "bar",
        "source"            : "joint_qc_table",
        "requires_samples"  : None,
        "x_var"             : "RUN_ID",
        "y_var"             : "PCT_Q30_R2",
        "fill_var"          : "RUN_ID",
        "x_lab"             : "Run ID",
        "y_lab"             : "Percentage of Q30 R2 reads",
        "title"             : "[LocalApp core run metric] Percentage of R2 reads with quality >= 30",
        "value_spec"        : {"operation": "cast", "column": "PCT_Q30_R2", "dtype": pl.Float64},
        "cart_ylim"         : (0, 104),
        "guideline"         : {
            "table"        : "joint_qc_guideline_table",
            "id_column"    : "RUN_ID",
            "sample_id"    : "LSL_Guideline",
            "value_spec"   : {"operation": "cast", "column": "PCT_Q30_R2", "dtype": pl.Float64},
            "python_cast"  : float,
            "label_prefix" : "LSL_Guideline",
            "alpha"        : 0.3,
            "color"        : "red",
        },
    },


    "DNA_CONTAMINATION_P_VALUE_BY_RUN": {


        "localapp": {
            "plot": True,
            "index": 6
        },
        "dragen": {
            "plot": False,
            "index": 0
        },

        "plot_kind"         : "contamination_scatter",
        "source"            : "dna_data_table",
        "requires_samples"  : "dna",
        "color_var"         : "RUN",
        "label_var"         : None,
        "title"             : "[LocalApp core DNA metric] Contamination assessment by TSO 500 LocalApp\n(samples located within the red box are deemed contaminated)",
        "skip_if_empty"     : True,
    },


    "DNA_CONTAMINATION_P_VALUE_HIGHLIGHTED_RUN": {

        "localapp": {
            "plot": True,
            "index": 7
        },
        "dragen": {
            "plot": False,
            "index": 0
        },
        
        "plot_kind"         : "contamination_scatter",
        "source"            : "dna_data_table",
        "requires_samples"  : "dna",
        "color_var"         : "highlighted_run",
        "label_var"         : "contamination_label",
        "guide_title"       : "Latest run",
        "color_values"      : {"False": "#888686", "True": "#C02F2F"},
        "title"             : "[LocalApp core DNA metric] Contamination assessment by TSO 500 LocalApp\n(samples located within the red box are deemed contaminated)",
        "skip_if_empty"     : True,
    },


    "DNA_COVERAGE_MAD": {

        "localapp": {
            "plot": True,
            "index": 8
        },
        "dragen": {
            "plot": False,
            "index": 0
        },
        
        "plot_kind"         : "bar",
        "source"            : "dna_data_table",
        "requires_samples"  : "dna",
        "x_var"             : "SAMPLE_ID",
        "y_var"             : "DNA_COVERAGE_MAD",
        "fill_var"          : "RUN",
        "x_lab"             : "Sample ID",
        "y_lab"             : "Coverage MAD (Median Absolute Deviation)",
        "title"             : "[LocalApp core DNA metric] Coverage MAD\n(median normalized deviation across all regions used for CNV calling)",
        "value_spec"        : {"operation": "cast", "column": "DNA_COVERAGE_MAD", "dtype": pl.Float64},
        "na_filter_columns" : ["DNA_COVERAGE_MAD"],
        "guideline"         : {
            "table"        : "dna_guideline_table",
            "sample_id"    : "USL_Guideline",
            "value_spec"   : {"operation": "cast", "column": "DNA_COVERAGE_MAD", "dtype": pl.Float64},
            "python_cast"  : float,
            "label_prefix" : "USL_Guideline",
            "ann_y_offset" : 0,
        },
    },


   "DNA_PCT_PF_UQ_READS": {

        "localapp": {
            "plot": True,
            "index": 20
        },
        "dragen": {
            "plot": False,
            "index": 0
        },
        
        "plot_kind"         : "bar",
        "source"            : "dna_data_table",
        "requires_samples"  : "dna",
        "x_var"             : "SAMPLE_ID",
        "y_var"             : "DNA_PCT_PF_UQ_READS",
        "fill_var"          : "RUN",
        "x_lab"             : "Sample ID",
        "y_lab"             : "Percentage of UQ reads passing filters",
        "title"             : "[LocalApp DNA metric] Percentage of UQ reads passing filters",
        "value_spec"        : {"operation": "cast", "column": "DNA_PCT_PF_UQ_READS", "dtype": pl.Float64},
        "na_filter_columns" : ["DNA_PCT_PF_UQ_READS"],
        "cart_ylim"         : (0, 104),
        "y_tick_step"       : 10,
    },



    # ========================================= #
    # Dragen-specific metrics (not in LocalApp) #
    # ========================================= #

    "DNA_GENE_SCALED_MAD": {

        "localapp": {
            "plot": False,
            "index": 0
        },
        "dragen": {
            "plot": True,
            "index": 1
        },

        "plot_kind"         : "bar",
        "source"            : "dna_data_table",
        "requires_samples"  : "dna",
        "x_var"             : "SAMPLE_ID",
        "y_var"             : "DNA_GENE_SCALED_MAD",
        "fill_var"          : "RUN",
        "x_lab"             : "Sample ID",
        "y_lab"             : "Gene-scaled MAD",
        "title"             : "[LocalApp core DNA metric] Gene-scaled MAD",
        "value_spec"        : {"operation": "cast", "column": "DNA_GENE_SCALED_MAD", "dtype": pl.Float64},
        "na_filter_columns" : ["DNA_GENE_SCALED_MAD"],
        "guideline"         : {
            "table"        : "dna_guideline_table",
            "sample_id"    : "USL_Guideline",
            "value_spec"   : {"operation": "cast", "column": "DNA_GENE_SCALED_MAD", "dtype": pl.Float64},
            "python_cast"  : float,
            "label_prefix" : "USL_Guideline",
            "ann_y_offset" : -0.04,
        },
    },


    "DNA_MEDIAN_TARGET_HRD_COVERAGE": {

        "localapp": {
            "plot": False,
            "index": 0
        },
        "dragen": {
            "plot": True,
            "index": 8
        },

        "plot_kind"         : "bar",
        "source"            : "dna_data_table",
        "requires_samples"  : "dna",
        "x_var"             : "SAMPLE_ID",
        "y_var"             : "DNA_MEDIAN_TARGET_HRD_COVERAGE",
        "fill_var"          : "RUN",
        "x_lab"             : "Sample ID",
        "y_lab"             : "Median target HRD coverage",
        "title"             : "[Dragen DNA metric] Median target HRD coverage",
        "value_spec"        : {"operation": "cast", "column": "DNA_MEDIAN_TARGET_HRD_COVERAGE", "dtype": pl.Float64},
        "na_filter_columns" : ["DNA_MEDIAN_TARGET_HRD_COVERAGE"],
        "y_tick_step"       : 100,
    },


    "DNA_PCT_Q30_BASES": {

        "localapp": {
            "plot": False,
            "index": 0
        },
        "dragen": {
            "plot": True,
            "index": 14
        },

        "plot_kind"         : "bar",
        "source"            : "dna_data_table",
        "requires_samples"  : "dna",
        "x_var"             : "SAMPLE_ID",
        "y_var"             : "DNA_PCT_Q30_BASES",
        "fill_var"          : "RUN",
        "x_lab"             : "Sample ID",
        "y_lab"             : "Percentage of Q30 bases",
        "title"             : "[Dragen DNA metric] Percentage of Q30 bases",
        "value_spec"        : {"operation": "cast", "column": "DNA_PCT_Q30_BASES", "dtype": pl.Float64},
        "na_filter_columns" : ["DNA_PCT_Q30_BASES"],
        "cart_ylim"         : (0, 104),
        "y_tick_step"       : 10,
    },


    "DNA_PCT_SOFT_CLIPPED_BASES": {

        "localapp": {
            "plot": False,
            "index": 0
        },
        "dragen": {
            "plot": True,
            "index": 16
        },

        "plot_kind"         : "bar",
        "source"            : "dna_data_table",
        "requires_samples"  : "dna",
        "x_var"             : "SAMPLE_ID",
        "y_var"             : "DNA_PCT_SOFT_CLIPPED_BASES",
        "fill_var"          : "RUN",
        "x_lab"             : "Sample ID",
        "y_lab"             : "Percentage of soft clipped bases",
        "title"             : "[Dragen DNA metric] Percentage of soft clipped bases ",
        "value_spec"        : {"operation": "cast", "column": "DNA_PCT_SOFT_CLIPPED_BASES", "dtype": pl.Float64},
        "na_filter_columns" : ["DNA_PCT_SOFT_CLIPPED_BASES"],
        "cart_ylim"         : (0, 100),
        "y_tick_step"       : 10,
        "overlay_labels"    : {
            "label_column" : "DNA_PCT_SOFT_CLIPPED_BASES",
            "angle"        : 90,
            "va"           : "center",
            "ha"           : "right",
            "nudge_y"      : 5,
            "nudge_x"      : 0.5, # Adapted to 10 runs, less runs would need a different nudge_x value
            "size"         : 7.5, # Adapted to 10 runs, less runs would need a different size value
        },
    },


    "DNA_PCT_TARGET_50X": {

        "localapp": {
            "plot": False,
            "index": 0
        },
        "dragen": {
            "plot": True,
            "index": 20
        },

        "plot_kind"         : "bar",
        "source"            : "dna_data_table",
        "requires_samples"  : "dna",
        "x_var"             : "SAMPLE_ID",
        "y_var"             : "DNA_PCT_TARGET_50X",
        "fill_var"          : "RUN",
        "x_lab"             : "Sample ID",
        "y_lab"             : "Percentage of target 50X",
        "title"             : "[Dragen DNA metric] Percentage of target 50X",
        "value_spec"        : {"operation": "cast", "column": "DNA_PCT_TARGET_50X", "dtype": pl.Float64},
        "na_filter_columns" : ["DNA_PCT_TARGET_50X"],
        "y_tick_step"       : 10,
    },


    "DNA_PCT_TARGET_HRD_50X": {

        "localapp": {
            "plot": False,
            "index": 0
        },
        "dragen": {
            "plot": True,
            "index": 21
        },

        "plot_kind"         : "bar",
        "source"            : "dna_data_table",
        "requires_samples"  : "dna",
        "x_var"             : "SAMPLE_ID",
        "y_var"             : "DNA_PCT_TARGET_HRD_50X",
        "fill_var"          : "RUN",
        "x_lab"             : "Sample ID",
        "y_lab"             : "Percentage of target HRD 50X",
        "title"             : "[Dragen DNA metric] Percentage of target HRD 50X",
        "value_spec"        : {"operation": "cast", "column": "DNA_PCT_TARGET_HRD_50X", "dtype": pl.Float64},
        "na_filter_columns" : ["DNA_PCT_TARGET_HRD_50X"],
        "y_tick_step"       : 10,
    },


    "RNA_GENE_ABOVE_MEDIAN_CUTOFF": {

        "localapp": {
            "plot": False,
            "index": 0
        },
        "dragen": {
            "plot": True,
            "index": 25
        },

        "plot_kind"         : "bar",
        "source"            : "rna_data_table",
        "requires_samples"  : "rna",
        "x_var"             : "SAMPLE_ID",
        "y_var"             : "RNA_GENE_ABOVE_MEDIAN_CUTOFF",
        "fill_var"          : "RUN",
        "x_lab"             : "Sample ID",
        "y_lab"             : "Percentage of genes above median cutoff",
        "title"             : "[Dragen RNA metric] Percentage of genes above median cutoff",
        "value_spec"        : {"operation": "cast", "column": "RNA_GENE_ABOVE_MEDIAN_CUTOFF", "dtype": pl.Float64},
        "na_filter_columns" : ["RNA_GENE_ABOVE_MEDIAN_CUTOFF"],
        "y_tick_step"       : 10,
    },


    "RNA_GENE_MEDIAN_COVERAGE": {

        "localapp": {
            "plot": False,
            "index": 0
        },
        "dragen": {
            "plot": True,
            "index": 26
        },

        "plot_kind"         : "bar",
        "source"            : "rna_data_table",
        "requires_samples"  : "rna",
        "x_var"             : "SAMPLE_ID",
        "y_var"             : "RNA_GENE_MEDIAN_COVERAGE",
        "fill_var"          : "RUN",
        "x_lab"             : "Sample ID",
        "y_lab"             : "Median gene coverage",
        "title"             : "[Dragen RNA metric] Median gene coverage",
        "value_spec"        : {"operation": "cast", "column": "RNA_GENE_MEDIAN_COVERAGE", "dtype": pl.Float64},
        "na_filter_columns" : ["RNA_GENE_MEDIAN_COVERAGE"],
        "y_tick_step"       : 100,
    },



    "RNA_MEDIAN_CV_GENE_500X": {


        "localapp": {
            "plot": False,
            "index": 0
        },
        "dragen": {
            "plot": True,
            "index": 27
        },

        "plot_kind"         : "bar",
        "source"            : "rna_data_table",
        "requires_samples"  : "rna",
        "x_var"             : "SAMPLE_ID",
        "y_var"             : "RNA_MEDIAN_CV_GENE_500X",
        "fill_var"          : "RUN",
        "x_lab"             : "Sample ID",
        "y_lab"             : "Median CV for genes with median coverage >500",
        "title"             : "[Dragen core RNA metric] Uniformity of coverage across transcripts",
        "value_spec"        : {"operation": "cast", "column": "RNA_MEDIAN_CV_GENE_500X", "dtype": pl.Float64},
        "na_filter_columns" : ["RNA_MEDIAN_CV_GENE_500X"],
        "cart_ylim"         : (0, 1.04),
        "y_tick_step"       : float(0.1),
        "guideline"         : {
            "table"        : "rna_guideline_table",
            "sample_id"    : "USL_Guideline",
            "value_spec"   : {"operation": "cast", "column": "RNA_MEDIAN_CV_GENE_500X", "dtype": pl.Float64},
            "python_cast"  : float,
            "label_prefix" : "USL_Guideline",
            "ann_y_offset" : -0.25,
        },
    },


    "RNA_MEDIAN_CV_GENE_500X": {

        "localapp": {
            "plot": True,
            "index": 28
        },
        "dragen": {
            "plot": False,
            "index": 0
        },

        "plot_kind"         : "bar",
        "source"            : "rna_data_table",
        "requires_samples"  : "rna",
        "x_var"             : "SAMPLE_ID",
        "y_var"             : "RNA_MEDIAN_CV_GENE_500X",
        "fill_var"          : "RUN",
        "x_lab"             : "Sample ID",
        "y_lab"             : "Median CV for genes with median coverage >500",
        "title"             : "[LocalApp core RNA metric] Uniformity of coverage across transcripts",
        "value_spec"        : {"operation": "cast", "column": "RNA_MEDIAN_CV_GENE_500X", "dtype": pl.Float64},
        "na_filter_columns" : ["RNA_MEDIAN_CV_GENE_500X"],
        "cart_ylim"         : (0, 125),
        "y_tick_step"       : 20,
        "guideline"         : {
            "table"        : "rna_guideline_table",
            "sample_id"    : "USL_Guideline",
            "value_spec"   : {"operation": "cast", "column": "RNA_MEDIAN_CV_GENE_500X", "dtype": pl.Float64},
            "python_cast"  : int,
            "label_prefix" : "USL_Guideline",
            "ann_y_offset" : 0,
        },
    },


    "RNA_PCT_Q30_BASES": {

       "localapp": {
            "plot": False,
            "index": 0
        },
        "dragen": {
            "plot": True,
            "index": 31
        },

        "plot_kind"         : "bar",
        "source"            : "rna_data_table",
        "requires_samples"  : "rna",
        "x_var"             : "SAMPLE_ID",
        "y_var"             : "RNA_PCT_Q30_BASES",
        "fill_var"          : "RUN",
        "x_lab"             : "Sample ID",
        "y_lab"             : "Percentage of Q30 bases",
        "title"             : "[Dragen RNA metric] Percentage of Q30 bases",
        "value_spec"        : {"operation": "cast", "column": "RNA_PCT_Q30_BASES", "dtype": pl.Float64},
        "na_filter_columns" : ["RNA_PCT_Q30_BASES"],
        "cart_ylim"         : (0, 100),
        "y_tick_step"       : 10,
    },


    "RNA_PCT_SOFT_CLIPPED_BASES": {


       "localapp": {
            "plot": False,
            "index": 0
        },
        "dragen": {
            "plot": True,
            "index": 32
        },

        "plot_kind"         : "bar",
        "source"            : "rna_data_table",
        "requires_samples"  : "rna",
        "x_var"             : "SAMPLE_ID",
        "y_var"             : "RNA_PCT_SOFT_CLIPPED_BASES",
        "fill_var"          : "RUN",
        "x_lab"             : "Sample ID",
        "y_lab"             : "Percentage of soft clipped bases",
        "title"             : "[Dragen RNA metric] Percentage of soft clipped bases",
        "value_spec"        : {"operation": "cast", "column": "RNA_PCT_SOFT_CLIPPED_BASES", "dtype": pl.Float64},
        "na_filter_columns" : ["RNA_PCT_SOFT_CLIPPED_BASES"],
        "cart_ylim"         : (0, 100),
        "y_tick_step"       : 10,
        "overlay_labels"    : {
            "label_column" : "RNA_PCT_SOFT_CLIPPED_BASES",
            "angle"        : 90,
            "va"           : "center",
            "ha"           : "right",
            "nudge_y"      : 5,
            "nudge_x"      : 0.5, # Adapted to 10 runs, less runs would need a different nudge_x value
            "size"         : 7.5, # Adapted to 10 runs, less runs would need a different size value
        },
    },



    # ============================================ #
    # Shared metrics (both in Dragen and LocalApp) #
    # ============================================ #

   "DNA_MEAN_FAMILY_SIZE": {

        "localapp": {
            "plot": True,
            "index": 9
        },
        "dragen": {
            "plot": True,
            "index": 2
        },
        
        "plot_kind"         : "bar",
        "source"            : "dna_data_table",
        "requires_samples"  : "dna",
        "x_var"             : "SAMPLE_ID",
        "y_var"             : "DNA_MEAN_FAMILY_SIZE",
        "fill_var"          : "RUN",
        "x_lab"             : "Sample ID",
        "y_lab"             : "Mean family size",
        "title"             : "[LocalApp DNA metric] Mean family size",
        "value_spec"        : {"operation": "cast", "column": "DNA_MEAN_FAMILY_SIZE", "dtype": pl.Float64},
        "na_filter_columns" : ["DNA_MEAN_FAMILY_SIZE"],
        "y_tick_step"       : 2,
        "cart_ylim_dynamic" : {"mode": "max_plus", "column": "DNA_MEAN_FAMILY_SIZE", "lower": 1, "offset": 10},
        "overlay_labels"    : {
            "label_column" : "DNA_MEAN_FAMILY_SIZE",
            "angle"        : 90,
            "va"           : "center",
            "ha"           : "center",
            "nudge_y"      : 0.85,
            "nudge_x"      : 0,
            "size"         : 9,
        },
        "skip_if_empty"     : True,
    },


    "DNA_MEAN_TARGET_COVERAGE": {

        "localapp": {
            "plot": True,
            "index": 10
        },
        "dragen": {
            "plot": True,
            "index": 3
        },

        "plot_kind"         : "bar",
        "source"            : "dna_data_table",
        "requires_samples"  : "dna",
        "x_var"             : "SAMPLE_ID",
        "y_var"             : "DNA_MEAN_TARGET_COVERAGE",
        "fill_var"          : "RUN",
        "x_lab"             : "Sample ID",
        "y_lab"             : "Mean target coverage",
        "title"             : "[InPreD core DNA metric] Mean target coverage",
        "value_spec"        : {"operation": "cast", "column": "DNA_MEAN_TARGET_COVERAGE", "dtype": pl.Float64},
        "na_filter_columns" : ["DNA_MEAN_TARGET_COVERAGE"],
        "y_tick_step"       : 100,
        "guideline"         : {
            "table"        : "internal_guideline_table",
            "sample_id"    : "Internal Guideline",
            "value_spec"   : {"operation": "cast", "column": "DNA_MEAN_TARGET_COVERAGE", "dtype": pl.Float64},
            "python_cast"  : float,
            "label_prefix" : "Internal Guideline",
            "alpha"        : 0.3,
            "color"        : "red",
            "ann_y_offset" : 0,
        },
    },


   "DNA_MEDIAN_BIN_COUNT_CNV_TARGET": {

        "localapp": {
            "plot": True,
            "index": 11
        },
        "dragen": {
            "plot": True,
            "index": 4
        },
        
        "plot_kind"         : "bar",
        "source"            : "dna_data_table",
        "requires_samples"  : "dna",
        "x_var"             : "SAMPLE_ID",
        "y_var"             : "DNA_MEDIAN_BIN_COUNT_CNV_TARGET",
        "fill_var"          : "RUN",
        "x_lab"             : "Sample ID",
        "y_lab"             : "Median raw bin count per CNV target",
        "title"             : "[LocalApp core DNA metric] Median raw bin count per CNV target",
        "value_spec"        : {"operation": "cast", "column": "DNA_MEDIAN_BIN_COUNT_CNV_TARGET", "dtype": pl.Float64},
        "na_filter_columns" : ["DNA_MEDIAN_BIN_COUNT_CNV_TARGET"],
        "guideline"         : {
            "table"        : "dna_guideline_table",
            "sample_id"    : "LSL_Guideline",
            "value_spec"   : {"operation": "cast", "column": "DNA_MEDIAN_BIN_COUNT_CNV_TARGET", "dtype": pl.Float64},
            "python_cast"  : float,
            "label_prefix" : "LSL_Guideline",
            "ann_y_offset" : 1,
        },
    },


    "DNA_MEDIAN_EXON_COVERAGE": {

        "localapp": {
            "plot": True,
            "index": 12
        },
        "dragen": {
            "plot": True,
            "index": 5
        },
        
        "plot_kind"         : "bar",
        "source"            : "dna_data_table",
        "requires_samples"  : "dna",
        "x_var"             : "SAMPLE_ID",
        "y_var"             : "DNA_MEDIAN_EXON_COVERAGE",
        "fill_var"          : "RUN",
        "x_lab"             : "Sample ID",
        "y_lab"             : "Median exon coverage",
        "title"             : "[LocalApp core DNA metric] Median exon coverage",
        "value_spec"        : {"operation": "cast", "column": "DNA_MEDIAN_EXON_COVERAGE", "dtype": pl.Float64},
        "na_filter_columns" : ["DNA_MEDIAN_EXON_COVERAGE"],
        "y_tick_step"       : 100,
        "guideline"         : {
            "table"        : "dna_guideline_table",
            "sample_id"    : "LSL_Guideline",
            "value_spec"   : {"operation": "cast", "column": "DNA_MEDIAN_EXON_COVERAGE", "dtype": pl.Float64},
            "python_cast"  : float,
            "label_prefix" : "LSL_Guideline",
            "ann_y_offset" : 200,
        },
    },


    "DNA_MEDIAN_INSERT_SIZE": {

        "localapp": {
            "plot": True,
            "index": 13
        },
        "dragen": {
            "plot": True,
            "index": 6
        },
        
        "plot_kind"         : "bar",
        "source"            : "dna_data_table",
        "requires_samples"  : "dna",
        "x_var"             : "SAMPLE_ID",
        "y_var"             : "DNA_MEDIAN_INSERT_SIZE",
        "fill_var"          : "RUN",
        "x_lab"             : "Sample ID",
        "y_lab"             : "Median insert size",
        "title"             : "[LocalApp core DNA metric] Median insert size",
        "value_spec"        : {"operation": "cast", "column": "DNA_MEDIAN_INSERT_SIZE", "dtype": pl.Float64},
        "na_filter_columns" : ["DNA_MEDIAN_INSERT_SIZE"],
        "guideline"         : {
            "table"        : "dna_guideline_table",
            "sample_id"    : "LSL_Guideline",
            "value_spec"   : {"operation": "cast", "column": "DNA_MEDIAN_INSERT_SIZE", "dtype": pl.Int32},
            "python_cast"  : int,
            "label_prefix" : "LSL_Guideline",
            "ann_y_offset" : 10,
        },
    },


    "DNA_MEDIAN_TARGET_COVERAGE": {

        "localapp": {
            "plot": True,
            "index": 14
        },
        "dragen": {
            "plot": True,
            "index": 7
        },

        "plot_kind"         : "bar",
        "source"            : "dna_data_table",
        "requires_samples"  : "dna",
        "x_var"             : "SAMPLE_ID",
        "y_var"             : "DNA_MEDIAN_TARGET_COVERAGE",
        "fill_var"          : "RUN",
        "x_lab"             : "Sample ID",
        "y_lab"             : "Median target coverage",
        "title"             : "[LocalApp DNA metric] Median target coverage",
        "value_spec"        : {"operation": "cast", "column": "DNA_MEDIAN_TARGET_COVERAGE", "dtype": pl.Float64},
        "na_filter_columns" : ["DNA_MEDIAN_TARGET_COVERAGE"],
        "y_tick_step"       : 100,
    },


    "DNA_PCT_ALIGNED_READS": {
        
        "localapp": {
            "plot": True,
            "index": 15
        },
        "dragen": {
            "plot": True,
            "index": 9
        },

        "plot_kind"         : "bar",
        "source"            : "dna_data_table",
        "requires_samples"  : "dna",
        "x_var"             : "SAMPLE_ID",
        "y_var"             : "DNA_PCT_ALIGNED_READS",
        "fill_var"          : "RUN",
        "x_lab"             : "Sample ID",
        "y_lab"             : "Percentage of aligned reads",
        "title"             : "[LocalApp DNA metric] Percentage of aligned reads",
        "value_spec"        : {"operation": "cast", "column": "DNA_PCT_ALIGNED_READS", "dtype": pl.Float64},
        "na_filter_columns" : ["DNA_PCT_ALIGNED_READS"],
        "cart_ylim"         : (0, 104),
        "y_tick_step"       : 10,
    },


   "DNA_PCT_CHIMERIC_READS": {

        "localapp": {
            "plot": True,
            "index": 16
        },
        "dragen": {
            "plot": True,
            "index": 10
        },

        "plot_kind"         : "bar",
        "source"            : "dna_data_table",
        "requires_samples"  : "dna",
        "x_var"             : "SAMPLE_ID",
        "y_var"             : "DNA_PCT_CHIMERIC_READS",
        "fill_var"          : "RUN",
        "x_lab"             : "Sample ID",
        "y_lab"             : "Percentage of chimeric reads",
        "title"             : "[LocalApp DNA metric] Percentage of chimeric reads",
        "value_spec"        : {"operation": "cast", "column": "DNA_PCT_CHIMERIC_READS", "dtype": pl.Float64},
        "na_filter_columns" : ["DNA_PCT_CHIMERIC_READS"],
    },


    "DNA_PCT_CONTAMINATION_EST": {

        "localapp": {
            "plot": True,
            "index": 17
        },
        "dragen": {
            "plot": True,
            "index": 11
        },
        
        "plot_kind"         : "bar",
        "source"            : "dna_data_table",
        "requires_samples"  : "dna",
        "x_var"             : "SAMPLE_ID",
        "y_var"             : "DNA_PCT_CONTAMINATION_EST",
        "fill_var"          : "RUN",
        "x_lab"             : "Sample ID",
        "y_lab"             : "Contamination estimate (percent)",
        "title"             : "[LocalApp DNA metric] Contamination estimate",
        "value_spec"        : {"operation": "cast", "column": "DNA_PCT_CONTAMINATION_EST", "dtype": pl.Float64},
        "na_filter_columns" : ["DNA_PCT_CONTAMINATION_EST"],
        "cart_ylim"         : (0, 104),
        "y_tick_step"       : 10,
        "overlay_labels"    : {
            "label_column" : "DNA_PCT_CONTAMINATION_EST",
            "angle"        : 90,
            "va"           : "center",
            "ha"           : "right",
            "nudge_y"      : 10,
            "nudge_x"      : 0.5,
            "size"         : 9,
        },
    },


    "DNA_PCT_EXON_100X": {

        "localapp": {
            "plot": True,
            "index": 18
        },
        "dragen": {
            "plot": True,
            "index": 12
        },
        
        "plot_kind"         : "bar",
        "source"            : "dna_data_table",
        "requires_samples"  : "dna",
        "x_var"             : "SAMPLE_ID",
        "y_var"             : "DNA_PCT_EXON_100X",
        "fill_var"          : "RUN",
        "x_lab"             : "Sample ID",
        "y_lab"             : "Percentage of exons at 100X coverage",
        "title"             : "[LocalApp DNA metric] Percentage of exons at 100X coverage",
        "value_spec"        : {"operation": "cast", "column": "DNA_PCT_EXON_100X", "dtype": pl.Float64},
        "na_filter_columns" : ["DNA_PCT_EXON_100X"],
        "cart_ylim"         : (0, 104),
        "y_tick_step"       : 10,
    },


   "DNA_PCT_EXON_50X": {

        "localapp": {
            "plot": True,
            "index": 19
        },
        "dragen": {
            "plot": True,
            "index": 13
        },
        
        "plot_kind"         : "bar",
        "source"            : "dna_data_table",
        "requires_samples"  : "dna",
        "x_var"             : "SAMPLE_ID",
        "y_var"             : "DNA_PCT_EXON_50X",
        "fill_var"          : "RUN",
        "x_lab"             : "Sample ID",
        "y_lab"             : "Percentage of exons with coverage >= 50X",
        "title"             : "[LocalApp core DNA metric] Percentage of exons with coverage >= 50X",
        "value_spec"        : {"operation": "cast", "column": "DNA_PCT_EXON_50X", "dtype": pl.Float64},
        "na_filter_columns" : ["DNA_PCT_EXON_50X"],
        "cart_ylim"         : (0, 104),
        "y_tick_step"       : 10,
        "guideline"         : {
            "table"        : "dna_guideline_table",
            "sample_id"    : "LSL_Guideline",
            "value_spec"   : {"operation": "cast", "column": "DNA_PCT_EXON_50X", "dtype": pl.Float64},
            "python_cast"  : int,
            "label_prefix" : "LSL_Guideline",
            "ann_y_offset" : -15,
        },
    },


    "DNA_PCT_READ_ENRICHMENT": {

        "localapp": {
            "plot": True,
            "index": 21
        },
        "dragen": {
            "plot": True,
            "index": 15
        },

        "plot_kind"         : "bar",
        "source"            : "dna_data_table",
        "requires_samples"  : "dna",
        "x_var"             : "SAMPLE_ID",
        "y_var"             : "DNA_PCT_READ_ENRICHMENT",
        "fill_var"          : "RUN",
        "x_lab"             : "Sample ID",
        "y_lab"             : "Read enrichment (percent)",
        "title"             : "[LocalApp DNA metric] Read enrichment",
        "value_spec"        : {"operation": "cast", "column": "DNA_PCT_READ_ENRICHMENT", "dtype": pl.Float64},
        "na_filter_columns" : ["DNA_PCT_READ_ENRICHMENT"],
        "cart_ylim"         : (0, 104),
        "y_tick_step"       : 10,
    },


    "DNA_PCT_TARGET_0_4X_MEAN": {

        "localapp": {
            "plot": True,
            "index": 22
        },
        "dragen": {
            "plot": True,
            "index": 17
        },

        "plot_kind"         : "bar",
        "source"            : "dna_data_table",
        "requires_samples"  : "dna",
        "x_var"             : "SAMPLE_ID",
        "y_var"             : "DNA_PCT_TARGET_0_4X_MEAN",
        "fill_var"          : "RUN",
        "x_lab"             : "Sample ID",
        "y_lab"             : "Percentage of target at 0.4 * mean coverage",
        "title"             : "[LocalApp DNA metric] Percentage of target at 0.4 * mean coverage",
        "value_spec"        : {"operation": "cast", "column": "DNA_PCT_TARGET_0_4X_MEAN", "dtype": pl.Float64},
        "na_filter_columns" : ["DNA_PCT_TARGET_0_4X_MEAN"],
        "cart_ylim"         : (0, 104),
        "y_tick_step"       : 10,
    },


   "DNA_PCT_TARGET_100X": {

        "localapp": {
            "plot": True,
            "index": 23
        },
        "dragen": {
            "plot": True,
            "index": 18
        },
        
        "plot_kind"         : "bar",
        "source"            : "dna_data_table",
        "requires_samples"  : "dna",
        "x_var"             : "SAMPLE_ID",
        "y_var"             : "DNA_PCT_TARGET_100X",
        "fill_var"          : "RUN",
        "x_lab"             : "Sample ID",
        "y_lab"             : "Percentage of target at 100X coverage",
        "title"             : "[LocalApp DNA metric] Percentage of target at 100X coverage",
        "value_spec"        : {"operation": "cast", "column": "DNA_PCT_TARGET_100X", "dtype": pl.Float64},
        "na_filter_columns" : ["DNA_PCT_TARGET_100X"],
        "cart_ylim"         : (0, 104),
        "y_tick_step"       : 10,
    },


    "DNA_PCT_TARGET_250X": {

        "localapp": {
            "plot": True,
            "index": 24
        },
        "dragen": {
            "plot": True,
            "index": 19
        },
        
        "plot_kind"         : "bar",
        "source"            : "dna_data_table",
        "requires_samples"  : "dna",
        "x_var"             : "SAMPLE_ID",
        "y_var"             : "DNA_PCT_TARGET_250X",
        "fill_var"          : "RUN",
        "x_lab"             : "Sample ID",
        "y_lab"             : "Percentage of target at 250X coverage",
        "title"             : "[LocalApp DNA metric] Percentage of target at 250X coverage",
        "value_spec"        : {"operation": "cast", "column": "DNA_PCT_TARGET_250X", "dtype": pl.Float64},
        "na_filter_columns" : ["DNA_PCT_TARGET_250X"],
        "cart_ylim"         : (0, 104),
        "y_tick_step"       : 10,
    },


    "DNA_PCT_USABLE_UMI_READS": {

        "localapp": {
            "plot": True,
            "index": 25
        },
        "dragen": {
            "plot": True,
            "index": 22
        },

        "plot_kind"         : "bar",
        "source"            : "dna_data_table",
        "requires_samples"  : "dna",
        "x_var"             : "SAMPLE_ID",
        "y_var"             : "DNA_PCT_USABLE_UMI_READS",
        "fill_var"          : "RUN",
        "x_lab"             : "Sample ID",
        "y_lab"             : "Usable UMI reads (percent)",
        "title"             : "[Dragen DNA metric] Usable UMI reads",
        "value_spec"        : {"operation": "cast", "column": "DNA_PCT_USABLE_UMI_READS", "dtype": pl.Float64},
        "na_filter_columns" : ["DNA_PCT_USABLE_UMI_READS"],
        "cart_ylim"         : (0, 104),
        "y_tick_step"       : 10,
    },


   "DNA_TOTAL_PF_READS_MILLIONS": {

        "localapp": {
            "plot": True,
            "index": 26
        },
        "dragen": {
            "plot": True,
            "index": 23
        },

        "plot_kind"         : "bar",
        "source"            : "dna_data_table",
        "requires_samples"  : "dna",
        "x_var"             : "SAMPLE_ID",
        "y_var"             : "DNA_TOTAL_PF_READS_MILLIONS",
        "fill_var"          : "RUN",
        "x_lab"             : "Sample ID",
        "y_lab"             : "Number of reads (millions)",
        "title"             : "[InPreD core DNA metric] Number of reads in millions",
        "value_spec"        : {
            "operation" : "divide",
            "column"    : "DNA_TOTAL_PF_READS",
            "dtype"     : pl.Float64,
            "divisor"   : 1000000.0,
        },
        "na_filter_columns" : ["DNA_TOTAL_PF_READS"],
        "y_tick_step"       : 10,
        "guideline"         : {
            "table"        : "internal_guideline_table",
            "sample_id"    : "Internal Guideline",
            "value_spec"   : {
                "operation" : "divide",
                "column"    : "DNA_TOTAL_PF_READS",
                "dtype"     : pl.Float64,
                "divisor"   : 1000000.0,
            },
            "python_cast"  : float,
            "label_prefix" : "Internal Guideline",
            "alpha"        : 0.3,
            "color"        : "red",
            "ann_y_offset" : -20,
        },
    },


    "DNA_USABLE_MSI_SITES": {

        "localapp": {
            "plot": True,
            "index": 27
        },
        "dragen": {
            "plot": True,
            "index": 24
        },
        
        "plot_kind"         : "bar",
        "source"            : "dna_data_table",
        "requires_samples"  : "dna",
        "x_var"             : "SAMPLE_ID",
        "y_var"             : "DNA_USABLE_MSI_SITES",
        "fill_var"          : "RUN",
        "x_lab"             : "Sample ID",
        "y_lab"             : "Usable MSI sites",
        "title"             : "[LocalApp core DNA metric] Number of usable MSI sites",
        "value_spec"        : {"operation": "cast", "column": "DNA_USABLE_MSI_SITES", "dtype": pl.Float64},
        "na_filter_columns" : ["DNA_USABLE_MSI_SITES"],
        "guideline"         : {
            "table"        : "dna_guideline_table",
            "sample_id"    : "LSL_Guideline",
            "value_spec"   : {"operation": "cast", "column": "DNA_USABLE_MSI_SITES", "dtype": pl.Int32},
            "python_cast"  : int,
            "label_prefix" : "LSL_Guideline",
            "ann_y_offset" : 0,
        },
    },


    "RNA_MEDIAN_INSERT_SIZE": {

        "localapp": {
            "plot": True,
            "index": 29
        },
        "dragen": {
            "plot": True,
            "index": 28
        },
        
        "plot_kind"         : "bar",
        "source"            : "rna_data_table",
        "requires_samples"  : "rna",
        "x_var"             : "SAMPLE_ID",
        "y_var"             : "RNA_MEDIAN_INSERT_SIZE",
        "fill_var"          : "RUN",
        "x_lab"             : "Sample ID",
        "y_lab"             : "Median insert size",
        "title"             : "[LocalApp core RNA metric] Median insert size",
        "value_spec"        : {"operation": "cast", "column": "RNA_MEDIAN_INSERT_SIZE", "dtype": pl.Float64},
        "na_filter_columns" : ["RNA_MEDIAN_INSERT_SIZE"],
        "y_tick_step"       : 10,
        "guideline"         : {
            "table"        : "rna_guideline_table",
            "sample_id"    : "LSL_Guideline",
            "value_spec"   : {"operation": "cast", "column": "RNA_MEDIAN_INSERT_SIZE", "dtype": pl.Float64},
            "python_cast"  : int,
            "label_prefix" : "LSL_Guideline",
            "ann_y_offset" : 0,
        },
    },


    "RNA_PCT_CHIMERIC_READS": {

        "localapp": {
            "plot": True,
            "index": 30
        },
        "dragen": {
            "plot": True,
            "index": 29
        },

        "plot_kind"         : "bar",
        "source"            : "rna_data_table",
        "requires_samples"  : "rna",
        "x_var"             : "SAMPLE_ID",
        "y_var"             : "RNA_PCT_CHIMERIC_READS",
        "fill_var"          : "RUN",
        "x_lab"             : "Sample ID",
        "y_lab"             : "Percentage of chimeric reads",
        "title"             : "[LocalApp RNA metric] Percentage of chimeric reads",
        "value_spec"        : {"operation": "cast", "column": "RNA_PCT_CHIMERIC_READS", "dtype": pl.Float64},
        "na_filter_columns" : ["RNA_PCT_CHIMERIC_READS"],
    },


    "RNA_PCT_ON_TARGET_READS": {

        "localapp": {
            "plot": True,
            "index": 31
        },
        "dragen": {
            "plot": True,
            "index": 30
        },

        "plot_kind"         : "bar",
        "source"            : "rna_data_table",
        "requires_samples"  : "rna",
        "x_var"             : "SAMPLE_ID",
        "y_var"             : "RNA_PCT_ON_TARGET_READS",
        "fill_var"          : "RUN",
        "x_lab"             : "Sample ID",
        "y_lab"             : "Percentage of on-target reads",
        "title"             : "[LocalApp RNA metric] Percentage of on-target reads",
        "value_spec"        : {"operation": "cast", "column": "RNA_PCT_ON_TARGET_READS", "dtype": pl.Float64},
        "na_filter_columns" : ["RNA_PCT_ON_TARGET_READS"],
        "cart_ylim"         : (0, 100),
        "y_tick_step"       : 10,
    },


    "RNA_SCALED_MEDIAN_GENE_COVERAGE": {

        "localapp": {
            "plot": True,
            "index": 32
        },
        "dragen": {
            "plot": True,
            "index": 33
        },
        
        "plot_kind"         : "bar",
        "source"            : "rna_data_table",
        "requires_samples"  : "rna",
        "x_var"             : "SAMPLE_ID",
        "y_var"             : "RNA_SCALED_MEDIAN_GENE_COVERAGE",
        "fill_var"          : "RUN",
        "x_lab"             : "Sample ID",
        "y_lab"             : "Scaled median gene coverage",
        "title"             : "[LocalApp RNA metric] Scaled median gene coverage",
        "value_spec"        : {"operation": "cast", "column": "RNA_SCALED_MEDIAN_GENE_COVERAGE", "dtype": pl.Float64},
        "na_filter_columns" : ["RNA_SCALED_MEDIAN_GENE_COVERAGE"],
    },


    "RNA_TOTAL_ON_TARGET_READS_MILLIONS": {
        
        "localapp": {
            "plot": True,
            "index": 33
        },
        "dragen": {
            "plot": True,
            "index": 34
        },

        "plot_kind"         : "bar",
        "source"            : "rna_data_table",
        "requires_samples"  : "rna",
        "x_var"             : "SAMPLE_ID",
        "y_var"             : "RNA_TOTAL_ON_TARGET_READS_MILLIONS",
        "fill_var"          : "RUN",
        "x_lab"             : "Sample ID",
        "y_lab"             : "On target reads (millions)",
        "title"             : "[LocalApp core RNA metric] Total on-target reads in millions",
        "value_spec"        : {
            "operation" : "divide",
            "column"    : "RNA_TOTAL_ON_TARGET_READS",
            "dtype"     : pl.Float64,
            "divisor"   : 1000000.0,
        },
        "na_filter_columns" : ["RNA_TOTAL_ON_TARGET_READS"],
        "guideline"         : {
            "table"        : "rna_guideline_table",
            "sample_id"    : "LSL_Guideline",
            "value_spec"   : {
                "operation" : "divide",
                "column"    : "RNA_TOTAL_ON_TARGET_READS",
                "dtype"     : pl.Float64,
                "divisor"   : 1000000.0,
            },
            "python_cast"  : int,
            "label_prefix" : "LSL_Guideline",
            "ann_y_offset" : 0,
        },
    },


    "RNA_TOTAL_PF_READS_MILLIONS": {

        "localapp": {
            "plot": True,
            "index": 34
        },
        "dragen": {
            "plot": True,
            "index": 35
        },
        
        "plot_kind"         : "bar",
        "source"            : "rna_data_table",
        "requires_samples"  : "rna",
        "x_var"             : "SAMPLE_ID",
        "y_var"             : "RNA_TOTAL_PF_READS_MILLIONS",
        "fill_var"          : "RUN",
        "x_lab"             : "Sample ID",
        "y_lab"             : "Total reads that pass quality filters (millions)",
        "title"             : "[LocalApp RNA metric]  Total filter-passing reads in millions\n(the input is down-sampled to 30 million reads [not read pairs])",
        "value_spec"        : {
            "operation" : "divide",
            "column"    : "RNA_TOTAL_PF_READS",
            "dtype"     : pl.Float64,
            "divisor"   : 1000000.0,
        },
        "na_filter_columns" : ["RNA_TOTAL_PF_READS"],
    },


    "DNA_READS_MILLIONS_PER_FAMILY_SIZE": {

       "localapp": {
            "plot": True,
            "index": 35
        },
        "dragen": {
            "plot": True,
            "index": 36
        },

        "plot_kind"         : "bar",
        "source"            : "dna_data_table",
        "requires_samples"  : "dna",
        "x_var"             : "SAMPLE_ID",
        "y_var"             : "DNA_READS_MILLIONS_PER_FAMILY_SIZE",
        "fill_var"          : "RUN",
        "x_lab"             : "Sample ID",
        "y_lab"             : "Number of reads (millions) / mean family size",
        "title"             : "[LocalApp DNA metric] Number of reads divided by mean read family size",
        "value_spec"        : {
            "operation"         : "ratio",
            "numerator"         : "DNA_TOTAL_PF_READS",
            "denominator"       : "DNA_MEAN_FAMILY_SIZE",
            "dtype"             : pl.Float64,
            "numerator_divisor" : 1000000.0,
        },
        "na_filter_columns" : ["DNA_TOTAL_PF_READS", "DNA_MEAN_FAMILY_SIZE"],
        "y_tick_step"       : 10,
    },


    "COMPLETED_ALL_STEPS_DNA": {

       "localapp": {
            "plot": True,
            "index": 36
        },
        "dragen": {
            "plot": True,
            "index": 37
        },

        "plot_kind"         : "bar",
        "source"            : "dna_data_table",
        "requires_samples"  : "dna",
        "x_var"             : "SAMPLE_ID",
        "y_var"             : "COMPLETED_ALL_STEPS",
        "fill_var"          : "RUN",
        "x_lab"             : "Sample ID",
        "y_lab"             : "Pipeline completion",
        "title"             : "[InPreD core DNA metric] Pipeline completion",
        "value_spec"        : {"operation": "cast", "column": "COMPLETED_ALL_STEPS"},
        "na_filter_columns" : ["COMPLETED_ALL_STEPS"],
        "filters"           : [{"column": "SAMPLE_ID", "contains": r"-D0\d-"}],
    },

    "COMPLETED_ALL_STEPS_RNA": {
       "localapp": {
            "plot": True,
            "index": 37
        },
        "dragen": {
            "plot": True,
            "index": 38
        },

        "plot_kind"         : "bar",
        "source"            : "rna_data_table",
        "requires_samples"  : "rna",
        "x_var"             : "SAMPLE_ID",
        "y_var"             : "COMPLETED_ALL_STEPS",
        "fill_var"          : "RUN",
        "x_lab"             : "Sample ID",
        "y_lab"             : "Pipeline completion",
        "title"             : "[InPreD core RNA metric] Pipeline completion",
        "value_spec"        : {"operation": "cast", "column": "COMPLETED_ALL_STEPS"},
        "na_filter_columns" : ["COMPLETED_ALL_STEPS"],
        "filters"           : [{"column": "SAMPLE_ID", "contains": r"-R0\d-"}],
    },

}








