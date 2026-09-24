# SAMPLE_METRICS_SPECS structure
# ================================
#
# Each top-level entry defines one QC metric row or a helper variable used in
# per-sample QC plots produced by plot_sample_metrics.py.
# The structure mirrors plot_specs_workflows.py as closely as possible.
#
# 1. Workflow routing sub-dicts  (localapp / dragen)
#    Each workflow has its own sub-dict.  Only include the sub-dict for workflows
#    where the entry is relevant.  A metric that exists in both workflows will
#    have both "localapp" and "dragen" sub-dicts; a workflow-exclusive metric
#    will have only one.
#
#    Keys within each routing sub-dict:
#    plot     : bool       — True to render this entry as a row in the plot.
#    index    : int        — Row position in the plot (1-based, top-to-bottom).
#                            Set to 0 when the entry is not plotted.
#    seq_type : list[str]  — Seq-type(s) for which this row is shown.
#                            Use ["dna"], ["rna"], or ["dna", "rna"].
#    load     : bool       — Optional.  True to load the variable from the TSV
#                            even when plot=False.  Used for helper variables that
#                            feed derived rows.  Defaults to False when omitted.
#
# 2. Row configuration fields
#
#    qc_type          Rendering mode for the row:
#                       "TrueFalse"  — coloured YES/NO boolean box.
#                       "FailPass"   — bar with fail (amber) on left, pass (green) on right.
#                       "PassFail"   — bar with pass (green) on left, fail (amber) on right.
#                       "NoThold"    — informational only; entire bar is white.
#                       None         — helper entry, never rendered.
#    label            Row display label (bold, left margin).  None for helpers.
#    threshold_label  Threshold description shown below the label.  None for helpers.
#    pass_value       TrueFalse rows only: the observed value that means PASS (green).
#    variable         TSV column name to load.  None for derived boolean rows.
#    lowest           Left (minimum) edge of the display bar.  None for TrueFalse / helpers.
#    threshold_source Where to read the QC threshold from at runtime:
#                       "lower"  — LOWER_THRESHOLD row (LSL_Guideline) in the input TSV.
#                                  Used for FailPass rows (value must be >= lower to pass).
#                       "upper"  — UPPER_THRESHOLD row (USL_Guideline) in the input TSV.
#                                  Used for PassFail rows (value must be <= upper to pass)
#                                  and contamination helper entries.
#                       None     — No threshold (TrueFalse, NoThold, or derived rows).
#                     If the TSV has NA for the variable's threshold, np.nan is returned
#                     and the row renders without a threshold line (equivalent to NoThold).
#    highest          Right (maximum) edge of the display bar.  None for TrueFalse / helpers.
#    float_format     True to display axis annotations with 1 decimal place.
#                     Defaults to False when omitted.
#    scale_factor     Divide the raw TSV value by this before display.
#                     Use 1_000_000 to convert raw read counts to millions.
#                     Also applied to the threshold read from the TSV.
#                     Omit (or set to 1) when no scaling is needed.
#
# To add a new metric row:
#   1. Add an entry here with the correct workflow sub-dict(s) and configuration.
#   2. No changes to plot_sample_metrics.py are required.
#
# To remove a metric row:
#   • Either delete the entry, or set plot=False in the relevant workflow sub-dict.
#
# To reorder rows:
#   • Change the index values.  Indexes must be unique per workflow × seq_type.
#
# To add Dragen-specific metrics or different thresholds:
#   • Add or update the "dragen" sub-dict on existing entries.
#   • Add new entries with only a "dragen" sub-dict for Dragen-exclusive metrics.
#   • Threshold values are always read from the workflow's own input TSV at runtime.


SAMPLE_METRICS_SPECS = {

    # ========================= #
    # Shared (both DNA and RNA) #
    # ========================= #

    "COMPLETED_ALL_STEPS": {

        "localapp": { 
            "plot"     : True, 
            "index"    : 1, 
            "seq_type" : ["dna", "rna"]
        },
        "dragen": {
            "plot"     : True,
            "index"    : 1,
            "seq_type" : ["dna", "rna"]
        },

        "qc_type"          : "TrueFalse",
        "label"            : "PIPELINE COMPLETED",
        "threshold_label"  : "",
        "pass_value"       : 1.0,      # obs == 1.0 → pipeline ran → PASS (green)
        "variable"         : "COMPLETED_ALL_STEPS",
        "lowest"           : None,
        "threshold_source" : None,
        "highest"          : None,
        "float_format"     : False,
    },


    # ================ #
    # DNA metrics rows #
    # ================ #

    "DNA_CONTAMINATION": {
        "localapp": {
            "plot"     : True,
            "index"    : 2,
            "seq_type" : ["dna"]
        },
        "dragen": {
            "plot"     : True,
            "index"    : 2,
            "seq_type" : ["dna"]
        },

        "qc_type"          : "TrueFalse",
        "label"            : "CONTAMINATION",
        "threshold_label"  : "",
        "pass_value"       : 0.0,      # obs == 0.0 → sample clean → PASS (green)
        "variable"         : None,     # derived — see compute_obs_values()
        "lowest"           : None,
        "threshold_source" : None,
        "highest"          : None,
        "float_format"     : False,
    },

    # Helper variables: loaded but not plotted; feed the DNA_CONTAMINATION derived row.
    "DNA_CONTAMINATION_SCORE": {
        "localapp": {"plot": False, "load": True, "index": 0, "seq_type": ["dna"]},
        "dragen":   {"plot": False, "load": True, "index": 0, "seq_type": ["dna"]},

        "qc_type":          None,
        "label":            None,
        "threshold_label":  None,
        "pass_value":       None,
        "variable":         "DNA_CONTAMINATION_SCORE",
        "lowest":           None,
        "threshold_source": "upper",
        "highest":          None,
        "float_format":     False,
    },

    "DNA_CONTAMINATION_P_VALUE": {
        "localapp": {"plot": False, "load": True, "index": 0, "seq_type": ["dna"]},
        "dragen":   {"plot": False, "load": True, "index": 0, "seq_type": ["dna"]},

        "qc_type":          None,
        "label":            None,
        "threshold_label":  None,
        "pass_value":       None,
        "variable":         "DNA_CONTAMINATION_P_VALUE",
        "lowest":           None,
        "threshold_source": "upper",
        "highest":          None,
        "float_format":     False,
    },

    "DNA_MEDIAN_EXON_COVERAGE": {
        "localapp": {
            "plot" : True,
            "index"    : 3,
            "seq_type" : ["dna"]
        },
        "dragen": {
            "plot"     : True,
            "index"    : 3,
            "seq_type" : ["dna"]
        },

        "qc_type"          : "FailPass",
        "label"            : "MEDIAN EXON COVERAGE",
        "threshold_label"  : ">= {threshold} [small variants, TMB]",
        "pass_value"       : None,
        "variable"         : "DNA_MEDIAN_EXON_COVERAGE",
        "lowest"           : 0.0,
        "threshold_source" : "lower",  # LOWER_THRESHOLD row
        "highest"          : 1000.0,
        "float_format"     : False,
    },


    "DNA_PCT_EXON_50X": {
        "localapp": {
            "plot"     : True,
            "index"    : 4,
            "seq_type" : ["dna"]},
        "dragen":   {
            "plot"     : True,
            "index"    : 4,
            "seq_type" : ["dna"]},

        "qc_type":          "FailPass",
        "label":            "EXON BASES WITH COVERAGE >=50 (%)",
        "threshold_label":  ">= {threshold} [small variants, TMB]",
        "pass_value":       None,
        "variable":         "DNA_PCT_EXON_50X",
        "lowest":           0.0,
        "threshold_source": "lower",
        "highest":          100.0,
        "float_format":     False,
    },


    "DNA_MEDIAN_INSERT_SIZE": {
        "localapp": {
            "plot"     : True,
            "index"    : 5,
            "seq_type" : ["dna"]},
        "dragen":   {
            "plot"     : True,
            "index"    : 5,
            "seq_type" : ["dna"]},

        "qc_type":          "FailPass",
        "label":            "MEDIAN INSERT SIZE",
        "threshold_label":  ">= {threshold} [small variants, TMB]",
        "pass_value":       None,
        "variable":         "DNA_MEDIAN_INSERT_SIZE",
        "lowest":           0.0,
        "threshold_source": "lower",
        "highest":          200.0,
        "float_format":     False,
    },


    "DNA_PCT_ALIGNED_READS": {
        "localapp": {
            "plot"     : True,
            "index"    : 6,
            "seq_type" : ["dna"]},
        "dragen":   {
            "plot"     : True,
            "index"    : 6,
            "seq_type" : ["dna"]},

        "qc_type":           "FailPass",
        "label":             "ALIGNED READS (%)",
        "threshold_label":   ">= {threshold}",
        "pass_value":        None,
        "variable":          "DNA_PCT_ALIGNED_READS",
        "lowest":            0.0,
        "threshold_source":  "lower",
        "threshold_fallback": 70.0,   # TSV has NA for this metric; use this value instead
        "highest":           100.0,
        "float_format":      False,
    },


    "DNA_USABLE_MSI_SITES": {
        "localapp": {
            "plot"     : True,
            "index"    : 7,
            "seq_type" : ["dna"]},
        "dragen":   {
            "plot"     : True,
            "index"    : 7,
            "seq_type" : ["dna"]},

        "qc_type":          "FailPass",
        "label":            "USABLE MSI SITES",
        "threshold_label":  ">= {threshold} [MSI]",
        "pass_value":       None,
        "variable":         "DNA_USABLE_MSI_SITES",
        "lowest":           0.0,
        "threshold_source": "lower",
        "highest":          140.0,
        "float_format":     False,
    },


    "DNA_MEDIAN_BIN_COUNT_CNV_TARGET": {
        "localapp": {
            "plot"     : True,
            "index"    : 8,
            "seq_type" : ["dna"]},
        "dragen":   {
            "plot"     : True,
            "index"    : 8,
            "seq_type" : ["dna"]},

        "qc_type":          "FailPass",
        "label":            "MEDIAN RAW BIN COUNT PER CNV TARGET",
        "threshold_label":  ">= {threshold}   [CNV]",
        "pass_value":       None,
        "variable":         "DNA_MEDIAN_BIN_COUNT_CNV_TARGET",
        "lowest":           0.0,
        "threshold_source": "lower",
        "highest":          20.0,
        "float_format":     False,
    },


    "DNA_COVERAGE_MAD": {
        "localapp": {
            "plot"     : True,
            "index"    : 9,
            "seq_type" : ["dna"]},
        "dragen":   {
            "plot"     : False,
            "index"    : 0,
            "seq_type" : ["dna"]},

        "qc_type":          "PassFail",  # lower is better; right half = fail
        "label":            "CNV COVERAGE MAD",
        "threshold_label":  "<= {threshold} [CNV]",
        "pass_value":       None,
        "variable":         "DNA_COVERAGE_MAD",
        "lowest":           0.0,
        "threshold_source": "upper",  # UPPER_THRESHOLD row
        "highest":          0.3,
        "float_format":     True,
    },


    "DNA_GENE_SCALED_MAD": {
        "localapp": {
            "plot"     : False,
            "index"    : 0,
            "seq_type" : ["dna"]},
        "dragen":   {
            "plot"     : True,
            "index"    : 9,
            "seq_type" : ["dna"]},

        "qc_type":          "PassFail",  # lower is better; right half = fail
        "label":            "CNV GENE-SCALED MAD",
        "threshold_label":  "<= {threshold} [CNV]",
        "pass_value":       None,
        "variable":         "DNA_GENE_SCALED_MAD",
        "lowest":           0.0,
        "threshold_source": "upper",  # UPPER_THRESHOLD row
        "highest":          0.3,
        "float_format":     True,
    },



    # ===================================== #
    # RNA rows                             #
    # ===================================== #

    "RNA_MEDIAN_CV_GENE_500X": {
        "localapp": {
            "plot"     : True,
            "index"    : 2,
            "seq_type" : ["rna"]},
        "dragen":   {
            "plot"     : True,
            "index"    : 2,
            "seq_type" : ["rna"]},

        "qc_type":          "PassFail",  # lower CV is better; right half = fail
        "label":            "CV FOR GENES AT MEDIAN COVERAGE >500",
        "threshold_label":  "<= {threshold} [fusions, splice variants]",
        "pass_value":       None,
        "variable":         "RNA_MEDIAN_CV_GENE_500X",
        "lowest":           0.0,
        "threshold_source": "upper",
        "highest":          130.0,
        "float_format":     False,
    },


    "RNA_TOTAL_ON_TARGET_READS": {
        "localapp": {
            "plot"     : True,
            "index"    : 3,
            "seq_type" : ["rna"]},
        "dragen":   {
            "plot"     : True,
            "index"    : 3,
            "seq_type" : ["rna"]},

        "qc_type":          "FailPass",
        "label":            "TOTAL ON-TARGET READS (MILLIONS)",
        "threshold_label":  ">= {threshold} [fusions, splice variants]",
        "pass_value":       None,
        "variable":         "RNA_TOTAL_ON_TARGET_READS",
        "lowest":           0.0,
        "threshold_source": "lower",  # raw count from TSV; scale_factor applied
        "highest":          30.0,     # expressed in millions (after scale_factor)
        "float_format":     False,
        "scale_factor":     1_000_000,
    },


    "RNA_MEDIAN_INSERT_SIZE": {
        "localapp": {
            "plot"     : True,
            "index"    : 4,
            "seq_type" : ["rna"]},
        "dragen":   {
            "plot"     : True,
            "index"    : 4,
            "seq_type" : ["rna"]},

        "qc_type":          "FailPass",
        "label":            "MEDIAN INSERT SIZE",
        "threshold_label":  ">= {threshold} [fusions, splice variants]",
        "pass_value":       None,
        "variable":         "RNA_MEDIAN_INSERT_SIZE",
        "lowest":           0.0,
        "threshold_source": "lower",
        "highest":          200.0,
        "float_format":     False,
    },


    "RNA_PCT_CHIMERIC_READS": {
        "localapp": {
            "plot"     : True,
            "index"    : 5,
            "seq_type" : ["rna"]},
        "dragen":   {
            "plot"     : True,
            "index"    : 5,
            "seq_type" : ["rna"]},

        "qc_type":          "NoThold",
        "label":            "CHIMERIC READS (%)",
        "threshold_label":  "No reference value",
        "pass_value":       None,
        "variable":         "RNA_PCT_CHIMERIC_READS",
        "lowest":           0.0,
        "threshold_source": None,
        "highest":          0.3,
        "float_format":     True,
    },


    "RNA_PCT_ON_TARGET_READS": {
        "localapp": {
            "plot"     : True,
            "index"    : 6,
            "seq_type" : ["rna"]},
        "dragen":   {
            "plot"     : True,
            "index"    : 6,
            "seq_type" : ["rna"]},

        "qc_type":          "NoThold",
        "label":            "ON-TARGET READS (%)",
        "threshold_label":  "No reference value",
        "pass_value":       None,
        "variable":         "RNA_PCT_ON_TARGET_READS",
        "lowest":           0.0,
        "threshold_source": None,
        "highest":          100.0,
        "float_format":     False,
    },


    "RNA_SCALED_MEDIAN_GENE_COVERAGE": {
        "localapp": {
            "plot"     : True,
            "index"    : 7,
            "seq_type" : ["rna"]},
        "dragen":   {
            "plot"     : True,
            "index"    : 7,
            "seq_type" : ["rna"]},

        "qc_type":          "NoThold",
        "label":            "SCALED MEDIAN GENE COVERAGE",
        "threshold_label":  "No reference value",
        "pass_value":       None,
        "variable":         "RNA_SCALED_MEDIAN_GENE_COVERAGE",
        "lowest":           0.0,
        "threshold_source": None,
        "highest":          7000.0,
        "float_format":     False,
    },


    "RNA_TOTAL_PF_READS": {
        "localapp": {
            "plot"     : True,
            "index"    : 8,
            "seq_type" : ["rna"]},
        "dragen":   {
            "plot"     : True,
            "index"    : 8,
            "seq_type" : ["rna"]},

        "qc_type":          "NoThold",
        "label":            "TOTAL READS PASSING FILTERS (MILLIONS)",
        "threshold_label":  "No reference value",
        "pass_value":       None,
        "variable":         "RNA_TOTAL_PF_READS",
        "lowest":           0.0,
        "threshold_source": None,
        "highest":          30.0,     # expressed in millions (after scale_factor)
        "float_format":     False,
        "scale_factor":     1_000_000,
    },
}
