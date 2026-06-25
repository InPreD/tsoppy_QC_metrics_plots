#!/usr/bin/env python3
"""
Pytest configuration and shared fixtures for plot_run_metrics tests.

This module provides reusable test fixtures that load real TSO500 data from the
data/ directory. Using actual production data prevents "cheating" by validating
that code handles real-world structures, missing values, and edge cases.

Fixture Strategy:
    - Load real TSV files from data/TSO500_postprocessing/
    - Cache dataframes to avoid repeated I/O during test runs
    - Provide filtered subsets (DNA data, RNA data, guidelines) for specific tests
    - Enable integration testing with authentic schema and data patterns
"""

from pathlib import Path
import pytest
import polars as pl


# ----------------------- #
# Data Path Configuration #
# ----------------------- #

def get_data_root() -> Path:
    """
    Resolve the data directory path relative to the repository root.
    
    Returns:
        Path to data/TSO500_postprocessing/
    """
    current_dir = Path(__file__).parent.parent  # Go up from tests/ to repo root
    data_root = current_dir / "data" / "TSO500_postprocessing" / "251211_NB552811_0117_AHTJM2BGYX"
    metrics_dir = data_root / "251211_NB552811_0117_AHTJM2BGYX_metrics_plots" / "intermediate_metrics_files"
    return metrics_dir


# ------------------------------------- #
# Core Data Fixtures (Real TSO500 Data) #
# ------------------------------------- #

@pytest.fixture(scope="session")
def metric_files_path() -> Path:
    """
    Path to the directory containing real metric TSV files.
    
    These are actual TSO500 metrics from the repository's data/ folder.
    Session scope: loaded once per test session for performance.
    
    Returns:
        Path to intermediate_metrics_files directory
    """
    return get_data_root()


@pytest.fixture(scope="session")
def single_run_metrics(metric_files_path) -> pl.DataFrame:
    """
    Load a single run's metrics file.
    
    Real TSO500 data containing both DNA and RNA metrics for one run.
    Includes guideline rows (LSL_Guideline, USL_Guideline) and sample rows.
    
    Columns (partial list):
        DNA: CONTAMINATION_SCORE, TOTAL_PF_READS, MEDIAN_INSERT_SIZE, etc.
        RNA: MEDIAN_CV_GENE_500X, TOTAL_PF_READS, etc.
        Meta: SAMPLE_ID, COMPLETED_ALL_STEPS
    
    Returns:
        Polars DataFrame with all rows and columns from a single run
    """
    metric_file = metric_files_path / "251211_NB552811_0117_AHTJM2BGYX_metrics.tsv"
    df = pl.read_csv(metric_file, separator="\t", infer_schema_length=0)
    return df


@pytest.fixture(scope="session")
def dna_samples_only(single_run_metrics) -> pl.DataFrame:
    """
    Filter for DNA-only samples (excludes RNA-only, includes samples with DNA data).
    
    DNA samples contain SAMPLE_ID matching -D01-, -D02-, etc. pattern (not guideline rows).
    Real TSO500 pattern: DNA samples labeled as "AU565-D01-C01-X06", "BT20-D01-C01-X06", etc.
    
    Returns:
        Filtered DataFrame with DNA samples only (6 samples)
    """
    df = single_run_metrics.filter(
        (pl.col("SAMPLE_ID") != "LSL_Guideline")
        & (pl.col("SAMPLE_ID") != "USL_Guideline")
        & (pl.col("SAMPLE_ID") != "Internal Guideline")
        & (pl.col("SAMPLE_ID").str.contains("-D\\d+-"))  # DNA samples: -D01-, -D02-, etc.
    )
    return df


@pytest.fixture(scope="session")
def rna_samples_only(single_run_metrics) -> pl.DataFrame:
    """
    Filter for RNA-only samples (excludes DNA-only, includes samples with RNA data).
    
    RNA samples contain SAMPLE_ID matching -R03-, -R04-, etc. pattern (not guideline rows).
    Real TSO500 pattern: RNA samples labeled as "AU565-R03-C01-X06", "BT20-R03-C01-X06", etc.
    
    Returns:
        Filtered DataFrame with RNA samples only (6 samples)
    """
    df = single_run_metrics.filter(
        (pl.col("SAMPLE_ID") != "LSL_Guideline")
        & (pl.col("SAMPLE_ID") != "USL_Guideline")
        & (pl.col("SAMPLE_ID") != "Internal Guideline")
        & (pl.col("SAMPLE_ID").str.contains("-R\\d+-"))  # RNA samples: -R03-, -R04-, etc.
    )
    return df


@pytest.fixture(scope="session")
def guideline_rows(single_run_metrics) -> pl.DataFrame:
    """
    Extract guideline rows (LSL_Guideline, USL_Guideline, Internal Guideline).
    
    Guidelines define QC thresholds and reference values used in plots.
    In real data, these are special rows mixed with sample rows.
    
    Returns:
        DataFrame containing only guideline rows
    """
    df = single_run_metrics.filter(
        (pl.col("SAMPLE_ID") == "LSL_Guideline")
        | (pl.col("SAMPLE_ID") == "USL_Guideline")
        | (pl.col("SAMPLE_ID") == "Internal Guideline")
    )
    return df


@pytest.fixture(scope="session")
def joint_qc_data(metric_files_path) -> pl.DataFrame:
    """
    Load joint sequencing QC data (run-level metrics).
    
    Contains sequencer run statistics: PCT_PF_READS, CLUSTER_DENSITY, ESTIMATED_YIELD, etc.
    One row per run, used for cross-run comparison plots.
    
    Columns:
        RUN_ID, PCT_PF_READS, PCT_Q30_R1, PCT_Q30_R2, CLUSTER_DENSITY,
        ESTIMATED_YIELD, CLUSTERS_PASSING_FILTER, RUN_NUMBER
    
    Returns:
        Polars DataFrame with run-level metrics
    """
    qc_file = metric_files_path / "joint_sequencing_QC_file.tsv"
    df = pl.read_csv(qc_file, separator="\t", infer_schema_length=0)
    return df


@pytest.fixture(scope="session")
def master_metrics_table(metric_files_path) -> pl.DataFrame:
    """
    Load pre-merged metrics table (all runs combined).
    
    This is the output of merge_metric_tables(): concatenated metrics from
    multiple runs with RUN column added for identification.
    Used to test filtering, transformation, and guideline extraction.
    
    Returns:
        Polars DataFrame with merged metrics from multiple runs
    """
    master_file = metric_files_path / "master_metrics_table.tsv"
    df = pl.read_csv(master_file, separator="\t", infer_schema_length=0)
    return df


# ============================================================================
# Derived Fixtures (Transformed Real Data)
# ============================================================================

@pytest.fixture
def metric_file_list(metric_files_path) -> list[Path]:
    """
    List of all individual metric TSV files in the test data directory.
    
    These are the input files to merge_metric_tables(). In production,
    the calling script provides these paths as argv.
    
    Returns:
        List of Path objects pointing to individual metric files
    """
    pattern = "*_metrics.tsv"
    files = sorted(metric_files_path.glob(pattern))
    return files


@pytest.fixture
def metric_file_paths_as_strings(metric_file_list) -> list[str]:
    """
    Metric file paths as strings (matching argv format).
    
    parse_args() receives file paths as strings in sys.argv.
    This fixture converts Path objects to strings for testing.
    
    Returns:
        List of string paths suitable for parse_args() testing
    """
    return [str(f) for f in metric_file_list]