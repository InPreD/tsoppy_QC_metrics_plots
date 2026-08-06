import polars as pl

PLOT_SPECS = {

    # Example of 'cluster_density_scatter' plot
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


    # Example 1 of 'contamination_scatter' plot
    # Each RUN_ID is colored differently, and samples within the red box are deemed contaminated
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


    # Example 2 of 'contamination_scatter' plot
    # This is a variant of the previous plot: here samples from the latest run are highlighted in red, and all other runs are greyed out
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


    # Example 1 of 'bar' plot
    # Barplots with labels on top
    "DNA_MEAN_FAMILY_SIZE": {
        "localapp": {
            "plot": True,
            "index": 21
        },
        "dragen": {
            "plot": True,
            "index": 14
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
            "ha"           : "right",
            "nudge_y"      : 0.85,
            "nudge_x"      : 0.5,
            "size"         : 9, # dragen has 7.5, maybe this can be programmatically handled instead of hardcoding it; not sure what it depends on but something like "if I have more samples smaller, less makes it bigger"
        },
        "skip_if_empty"     : True,
    },


    # Example 2 of 'bar' plot
    # Barplots with guidelines lines 
    "DNA_MEDIAN_EXON_COVERAGE": {
        "localapp": {
            "plot": True,
            "index": 8
        },
        "dragen": {
            "plot": True,
            "index": 1
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
            "value_spec"   : {"operation": "cast", "column": "DNA_MEDIAN_EXON_COVERAGE", "dtype": pl.Int32},
            "python_cast"  : int,
            "label_prefix" : "LSL_Guideline",
            "ann_y_offset" : 200,
        },
    },

    # Example if they differ too much
    "DNA_MEDIAN_EXON_COVERAGE_LOCALAPP": {
        "localapp": {
            "plot": True,
            "index": 9
        },
        "dragen": {
            "plot": False,
            "index": 0
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
            "value_spec"   : {"operation": "cast", "column": "DNA_MEDIAN_EXON_COVERAGE", "dtype": pl.Int32},
            "python_cast"  : int,
            "label_prefix" : "LSL_Guideline",
            "ann_y_offset" : 200,
        },
    },

    "DNA_MEDIAN_EXON_COVERAGE_DRAGEN": {
        "localapp": {
            "plot": False,    
            "index": 0
        },
        "dragen": {
            "plot": True,
            "index": 2
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
        "y_tick_step"       : 200,
        "guideline"         : {
            "table"        : "dna_guideline_table",
            "sample_id"    : "LSL_Guideline",
            "value_spec"   : {"operation": "cast", "column": "DNA_MEDIAN_EXON_COVERAGE", "dtype": pl.Int32},
            "python_cast"  : int,
            "label_prefix" : "LSL_Guideline",
            "ann_y_offset" : 800,
        },
    },

}