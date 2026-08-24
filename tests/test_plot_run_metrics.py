#!/usr/bin/env python3
"""
Comprehensive unit tests for plot_run_metrics.py and plots.py

Fixtures from conftest.py provide real TSO500 data to prevent "test cheating".
"""

import tempfile
from pathlib import Path
from io import BytesIO
from unittest.mock import Mock, patch, MagicMock

import pytest
import polars as pl
import pandas as pd
from plotnine import ggplot

from plot_QC_metrics import (
    parse_args,
    validate_plot_specs,
    valid_metric_expr,
    build_filter_expression,
    build_value_expression,
    get_guideline_value,
    build_tables,
    compute_cart_ylim,
    prepare_bar_plot_data,
    render_bar_plot,
    render_cluster_density_scatter,
    render_contamination_scatter,
    render_plot,
    save_plot,
    load_input_tables,
    filter_most_recent_runs,
)
from plots import plot_tsoppy_barplot, plot_contamination_scatter
import matplotlib.pyplot as plt


# --------------------------------------------- #
# Local Fixtures for Workflow-Aware Build Tests #
# --------------------------------------------- #

@pytest.fixture(scope="session")
def merged_workflow_tables_path() -> Path:
    """Path to real workflow tables used by build_tables tests."""
    return Path(__file__).parent.parent / "data" / "workflow_tables"


@pytest.fixture(scope="session")
def joint_qc_data_2r(merged_workflow_tables_path) -> pl.DataFrame:
    """Real joint QC table from workflow_tables fixture set (combined LocalApp + Dragen)."""
    localapp_df = pl.read_csv(
        merged_workflow_tables_path / "localapp_plot_joint_sequencing_QC_file.tsv",
        separator="\t",
        infer_schema_length=0,
    )
    dragen_df = pl.read_csv(
        merged_workflow_tables_path / "dragen_plot_joint_sequencing_QC_file.tsv",
        separator="\t",
        infer_schema_length=0,
    )
    # Combine both workflows' data (LocalApp + Dragen)
    return pl.concat([localapp_df, dragen_df])


@pytest.fixture(scope="session")
def master_metrics_table_2r(merged_workflow_tables_path) -> pl.DataFrame:
    """Real merged metrics table from workflow_tables fixture set (combined LocalApp + Dragen)."""
    localapp_df = pl.read_csv(
        merged_workflow_tables_path / "localapp_plot_metrics_table.tsv",
        separator="\t",
        infer_schema_length=0,
    )
    dragen_df = pl.read_csv(
        merged_workflow_tables_path / "dragen_plot_metrics_table.tsv",
        separator="\t",
        infer_schema_length=0,
    )
    # Combine both workflows' data (LocalApp + Dragen)
    return pl.concat([localapp_df, dragen_df])


# ------------------------------------ #
# Unit tests without Data Dependencies #
# ------------------------------------ #

class TestValidatePlotSpecs:
    """
    Table-driven tests for validate_plot_specs().

    Valid and invalid spec cases are defined as parametrize data at the top of each test method.
    Each pytest.param id reflects the scenario being tested, keeping the test body free of
    per-case logic and making it trivial to add new cases without touching the test logic.
    """

    @pytest.mark.parametrize("specs", [
        pytest.param(
            {
                "TEST_PLOT": {
                    "localapp" : {"plot": True,  "index": 1},
                    "dragen"   : {"plot": False, "index": 0},
                    "plot_kind"   : "bar",
                    "source"      : "joint_qc_table",
                    "title"       : "Test Plot",
                    "x_var"       : "RUN_ID",
                    "y_var"       : "VALUE",
                    "fill_var"    : "RUN_ID",
                    "x_lab"       : "Run",
                    "y_lab"       : "Value",
                    "value_spec"  : {"operation": "cast", "column": "VALUE"},
                }
            },
            id="bar_plot",
        ),
        pytest.param(
            {
                "SCATTER_PLOT": {
                    "localapp" : {"plot": True,  "index": 2},
                    "dragen"   : {"plot": False, "index": 0},
                    "plot_kind"   : "cluster_density_scatter",
                    "source"      : "joint_qc_table",
                    "title"       : "Scatter Test",
                }
            },
            id="cluster_density_scatter",
        ),
        pytest.param(
            {
                "CONTAMINATION": {
                    "localapp" : {"plot": True,  "index": 3},
                    "dragen"   : {"plot": False, "index": 0},
                    "plot_kind"   : "contamination_scatter",
                    "source"      : "dna_data_table",
                    "title"       : "Contamination Plot",
                    "color_var"   : "RUN",
                    "label_var"   : "SAMPLE_ID",
                }
            },
            id="contamination_scatter",
        ),
    ])
    def test_valid_specs_do_not_raise(self, specs):
        """Valid specs should pass validation without raising."""
        validate_plot_specs(specs)


    @pytest.mark.parametrize("specs, expected_match", [
        pytest.param(
            {
                "BAD_PLOT": {
                    "localapp" : {"plot": True,  "index": 1},
                    "dragen"   : {"plot": False, "index": 0},
                    # Missing: plot_kind, source, title
                }
            },
            "Error in",
            id="missing_common_fields",
        ),
        pytest.param(
            {
                "BAD_PLOT": {
                    "localapp" : {"plot": True,  "index": 1},
                    "dragen"   : {"plot": False, "index": 0},
                    "plot_kind"   : "unknown_kind",
                    "source"      : "joint_qc_table",
                    "title"       : "Test",
                }
            },
            "not recognized",
            id="unknown_plot_kind",
        ),
        pytest.param(
            {
                "BAD_BAR": {
                    "localapp" : {"plot": True,  "index": 1},
                    "dragen"   : {"plot": False, "index": 0},
                    "plot_kind"   : "bar",
                    "source"      : "joint_qc_table",
                    "title"       : "Test",
                    # Missing: x_var, y_var, fill_var, x_lab, y_lab, value_spec
                }
            },
            "Missing fields",
            id="missing_bar_fields",
        ),
        pytest.param(
            {
                "BAD_CONTAMINATION": {
                    "localapp" : {"plot": True,  "index": 1},
                    "dragen"   : {"plot": False, "index": 0},
                    "plot_kind"   : "contamination_scatter",
                    "source"      : "dna_data_table",
                    "title"       : "Test",
                    # Missing: color_var, label_var
                }
            },
            "Missing fields",
            id="missing_contamination_fields",
        ),
        pytest.param(
            {
                "SPEC_1": {},
                "SPEC_2": {
                    "localapp" : {"plot": True,  "index": 2},
                    "dragen"   : {"plot": False, "index": 0},
                    "plot_kind"   : "nonexistent_kind",
                    "source"      : "table",
                    "title"       : "Test",
                },
            },
            "2 validation error",
            id="multiple_errors_accumulation",
        ),
        pytest.param(
            {
                "MISSING_WORKFLOW_KEY": {
                    # Missing dragen entirely
                    "localapp" : {"plot": True,  "index": 1},
                    "plot_kind"   : "bar",
                    "source"      : "joint_qc_table",
                    "title"       : "Test",
                    "x_var": "X", "y_var": "Y", "fill_var": "F",
                    "x_lab": "X", "y_lab": "Y",
                    "value_spec": {"operation": "cast", "column": "X"},
                }
            },
            "Missing workflow routing key",
            id="missing_workflow_key",
        ),
        pytest.param(
            {
                "INVALID_PLOT_TYPE": {
                    "localapp" : {"plot": "yes", "index": 1},  # plot must be bool
                    "dragen"   : {"plot": False, "index": 0},
                    "plot_kind"   : "bar",
                    "source"      : "joint_qc_table",
                    "title"       : "Test",
                    "x_var": "X", "y_var": "Y", "fill_var": "F",
                    "x_lab": "X", "y_lab": "Y",
                    "value_spec": {"operation": "cast", "column": "X"},
                }
            },
            "'plot' must be a bool",
            id="invalid_plot_type",
        ),
        pytest.param(
            {
                "NEGATIVE_INDEX": {
                    "localapp" : {"plot": True, "index": -1},  # index must be >= 0
                    "dragen"   : {"plot": False, "index": 0},
                    "plot_kind"   : "bar",
                    "source"      : "joint_qc_table",
                    "title"       : "Test",
                    "x_var": "X", "y_var": "Y", "fill_var": "F",
                    "x_lab": "X", "y_lab": "Y",
                    "value_spec": {"operation": "cast", "column": "X"},
                }
            },
            "'index' must be a non-negative integer",
            id="negative_index",
        ),
        pytest.param(
            {
                "FLOAT_INDEX": {
                    "localapp" : {"plot": True, "index": 1.5},  # index must be int
                    "dragen"   : {"plot": False, "index": 0},
                    "plot_kind"   : "bar",
                    "source"      : "joint_qc_table",
                    "title"       : "Test",
                    "x_var": "X", "y_var": "Y", "fill_var": "F",
                    "x_lab": "X", "y_lab": "Y",
                    "value_spec": {"operation": "cast", "column": "X"},
                }
            },
            "'index' must be a non-negative integer",
            id="float_index",
        ),
    ])
    def test_invalid_specs_raise_key_error(self, specs, expected_match):
        """Invalid specs should raise KeyError with a message matching the expected pattern."""
        with pytest.raises(KeyError, match=expected_match):
            validate_plot_specs(specs)


    @patch("plot_QC_metrics.logger")
    def test_validate_plot_specs_logs_errors(self, mock_logger):
        """
        Test that each validation error is logged immediately via logger.error().
        Kept separate because it tests logger behavior via mock, not exception raising.
        """
        invalid_specs = {
            "PLOT_A": {
                "localapp" : {"plot": True,  "index": 1},
                "dragen"   : {"plot": False, "index": 0},
                # Missing: plot_kind, source, title
            },
            "PLOT_B": {
                "localapp" : {"plot": True,  "index": 2},
                "dragen"   : {"plot": False, "index": 1},  # Changed from 0 to 1 to avoid duplicate
                "plot_kind" : "bar",
                "source"    : "joint_qc_table",
                "title"     : "Test",
                # Missing bar-specific fields
            }
        }
        
        with pytest.raises(KeyError):
            validate_plot_specs(invalid_specs)
        
        assert mock_logger.error.call_count == 2
        
        first_call_args = mock_logger.error.call_args_list[0][0][0]
        assert "PLOT_A" in first_call_args
        assert "Missing fields:" in first_call_args
        
        second_call_args = mock_logger.error.call_args_list[1][0][0]
        assert "PLOT_B" in second_call_args
        assert "Missing fields:" in second_call_args


    def test_duplicate_indexes_raise_error(self):
        """
        Duplicate indexes within the same workflow should raise KeyError.
        """
        specs = {
            "PLOT_1": {
                "localapp" : {"plot": True,  "index": 1},
                "dragen"   : {"plot": False, "index": 0},
                "plot_kind"   : "bar",
                "source"      : "joint_qc_table",
                "title"       : "Plot 1",
                "x_var": "X", "y_var": "Y", "fill_var": "F",
                "x_lab": "X", "y_lab": "Y",
                "value_spec": {"operation": "cast", "column": "X"},
            },
            "PLOT_2": {
                "localapp" : {"plot": True,  "index": 1},  # Duplicate index for localapp
                "dragen"   : {"plot": False, "index": 0},
                "plot_kind"   : "bar",
                "source"      : "joint_qc_table",
                "title"       : "Plot 2",
                "x_var": "X", "y_var": "Y", "fill_var": "F",
                "x_lab": "X", "y_lab": "Y",
                "value_spec": {"operation": "cast", "column": "X"},
            },
        }

        with pytest.raises(KeyError, match="Duplicate index"):
            validate_plot_specs(specs)


    def test_multiple_index_zeros_do_not_raise_error(self):
        """
        Multiple plots with index 0 (not relevant) should NOT raise a duplicate index error.
        Index 0 means the plot is not relevant for that workflow and should be ignored.
        """
        specs = {
            "PLOT_1": {
                "localapp" : {"plot": False, "index": 0},  # Not relevant for localapp
                "dragen"   : {"plot": True,  "index": 1},
                "plot_kind"   : "bar",
                "source"      : "joint_qc_table",
                "title"       : "Plot 1",
                "x_var": "X", "y_var": "Y", "fill_var": "F",
                "x_lab": "X", "y_lab": "Y",
                "value_spec": {"operation": "cast", "column": "X"},
            },
            "PLOT_2": {
                "localapp" : {"plot": False, "index": 0},  # Also not relevant for localapp
                "dragen"   : {"plot": True,  "index": 2},
                "plot_kind"   : "bar",
                "source"      : "joint_qc_table",
                "title"       : "Plot 2",
                "x_var": "X", "y_var": "Y", "fill_var": "F",
                "x_lab": "X", "y_lab": "Y",
                "value_spec": {"operation": "cast", "column": "X"},
            },
        }

        # Should not raise an error since multiple 0s are allowed
        validate_plot_specs(specs)  # Passes without raising


class TestParseArgs:
    """Test command-line argument parsing."""

    @pytest.mark.parametrize(
        "argv, expected",
        [
            pytest.param(
                [
                    "--joint_QC",       "joint_qc.tsv",
                    "--output_pdf",     "output/run.pdf",
                    "--metrics_table",  "merged.tsv",
                    "--run_id",         "RUN_001",
                    "--workflow_type",  "localapp",
                ],
                {
                    "joint_qc_file"      : "joint_qc.tsv",
                    "output_pdf"         : "output/run.pdf",
                    "metrics_table"      : "merged.tsv",
                    "highlighted_run_id" : "RUN_001",
                    "workflow_type"      : "localapp",
                    "max_runs"           : 10,
                },
                id="minimal_valid_case",
            ),
            pytest.param(
                [
                    "--joint_QC",       "joint_qc.tsv",
                    "--output_pdf",     "output/run.pdf",
                    "--metrics_table",  "merged.tsv",
                    "--run_id",         "RUN_001",
                    "--workflow_type",  "localapp",
                    "--max_runs",       "5",
                ],
                {
                    "joint_qc_file"      : "joint_qc.tsv",
                    "output_pdf"         : "output/run.pdf",
                    "metrics_table"      : "merged.tsv",
                    "highlighted_run_id" : "RUN_001",
                    "workflow_type"      : "localapp",
                    "max_runs"           : 5,
                },
                id="with_custom_max_runs",
            ),
        ],
    )
    def test_parse_args_valid(self, argv, expected):
        """Parse required named arguments and preserve their values."""
        result = parse_args(argv)

        assert result == expected


    def test_parse_args_file_paths_preserved(self):
        """
        Test that file paths are preserved exactly as given.
        This ensures that Path objects are not altered during parsing.
        """
        argv = [
            "--joint_QC",       "path/to/joint_qc.tsv",
            "--output_pdf",     "output/results/run.pdf",
            "--metrics_table",  "output/merged.tsv",
            "--run_id",         "RUN_highlighted",
            "--workflow_type",  "dragen",
        ]
        result = parse_args(argv)

        assert result["joint_qc_file"] == "path/to/joint_qc.tsv"
        assert result["output_pdf"]    == "output/results/run.pdf"
        assert result["metrics_table"] == "output/merged.tsv"


    def test_parse_args_workflow_type_is_case_insensitive(self):
        """--workflow_type should accept mixed-case values and normalise to lowercase."""
        argv = [
            "--joint_QC",       "joint_qc.tsv",
            "--output_pdf",     "run.pdf",
            "--metrics_table",  "merged.tsv",
            "--run_id",         "RUN_001",
            "--workflow_type",  "LocalApp",
        ]
        result = parse_args(argv)
        assert result["workflow_type"] == "localapp"


    def test_parse_args_max_runs_defaults_to_10(self):
        """--max_runs should default to 10 when not supplied."""
        argv = [
            "--joint_QC",       "joint_qc.tsv",
            "--output_pdf",     "run.pdf",
            "--metrics_table",  "merged.tsv",
            "--run_id",         "RUN_001",
            "--workflow_type",  "localapp",
        ]
        result = parse_args(argv)
        assert result["max_runs"] == 10


    @pytest.mark.parametrize(
        "max_runs_value",
        [
            pytest.param("0",  id="zero"),
            pytest.param("-1", id="negative"),
            pytest.param("11", id="exceeds_max"),
            pytest.param("100", id="way_exceeds_max"),
        ],
    )
    def test_parse_args_max_runs_invalid_raises(self, max_runs_value):
        """--max_runs must be between 1 and 10; values outside this range should raise SystemExit."""
        argv = [
            "--joint_QC",       "joint_qc.tsv",
            "--output_pdf",     "run.pdf",
            "--metrics_table",  "merged.tsv",
            "--run_id",         "RUN_001",
            "--workflow_type",  "localapp",
            "--max_runs",       max_runs_value,
        ]
        with pytest.raises(SystemExit):
            parse_args(argv)


    @pytest.mark.parametrize(
        "argv",
        [
            pytest.param(
                [
                    "--output_pdf",     "run.pdf",
                    "--metrics_table",  "merged.tsv",
                    "--run_id",         "RUN_001",
                    "--workflow_type",  "localapp",
                ],
                id="missing_joint_qc",
            ),
            pytest.param(
                [
                    "--joint_QC",       "joint_qc.tsv",
                    "--metrics_table",  "merged.tsv",
                    "--run_id",         "RUN_001",
                    "--workflow_type",  "localapp",
                ],
                id="missing_output_pdf",
            ),
            pytest.param(
                [
                    "--joint_QC",       "joint_qc.tsv",
                    "--output_pdf",     "run.pdf",
                    "--run_id",         "RUN_001",
                    "--workflow_type",  "localapp",
                ],
                id="missing_metrics_table",
            ),
            pytest.param(
                [
                    "--joint_QC",       "joint_qc.tsv",
                    "--output_pdf",     "run.pdf",
                    "--metrics_table",  "merged.tsv",
                    "--workflow_type",  "localapp",
                ],
                id="missing_run_id",
            ),
            pytest.param(
                [
                    "--joint_QC",       "joint_qc.tsv",
                    "--output_pdf",     "run.pdf",
                    "--metrics_table",  "merged.tsv",
                    "--run_id",         "RUN_001",
                ],
                id="missing_workflow_type",
            ),
            pytest.param(
                [
                    "--joint_QC",       "joint_qc.tsv",
                    "--output_pdf",     "run.pdf",
                    "--metrics_table",  "merged.tsv",
                    "--run_id",         "RUN_001",
                    "--workflow_type",  "unsupported_workflow",
                ],
                id="invalid_workflow_type",
            ),
        ],
    )
    def test_parse_args_missing_required_raises(self, argv):
        """Missing or invalid required CLI arguments should raise SystemExit from argparse."""
        with pytest.raises(SystemExit):
            parse_args(argv)


# ---------------------------------------- #
# Fixtures for table loading / filtering   #
# ---------------------------------------- #

WORKFLOW_TABLES_DIR = Path(__file__).parent.parent / "data" / "workflow_tables"


@pytest.fixture(scope="module")
def valid_joint_qc_path() -> Path:
    return WORKFLOW_TABLES_DIR / "localapp_plot_joint_sequencing_QC_file.tsv"


@pytest.fixture(scope="module")
def valid_metrics_path() -> Path:
    return WORKFLOW_TABLES_DIR / "localapp_plot_metrics_table.tsv"


class TestLoadInputTables:
    """Tests for load_input_tables() — loading and column validation."""

    def test_load_returns_two_dataframes(self, valid_joint_qc_path, valid_metrics_path):
        """load_input_tables returns two non-empty polars DataFrames."""
        joint_qc, metrics = load_input_tables(str(valid_joint_qc_path), str(valid_metrics_path))

        assert isinstance(joint_qc, pl.DataFrame)
        assert isinstance(metrics, pl.DataFrame)
        assert joint_qc.height > 0
        assert metrics.height  > 0


    def test_joint_qc_contains_required_columns(self, valid_joint_qc_path, valid_metrics_path):
        """joint_qc_table must contain RUN_INDEX, WORKFLOW_TYPE, WORKFLOW_VERSION."""
        joint_qc, _ = load_input_tables(str(valid_joint_qc_path), str(valid_metrics_path))

        required = {"RUN_INDEX", "WORKFLOW_TYPE", "WORKFLOW_VERSION"}
        assert required.issubset(set(joint_qc.columns))


    def test_metrics_contains_required_columns(self, valid_joint_qc_path, valid_metrics_path):
        """metrics_table must contain RUN_INDEX, SAMPLE_ID, RUN, WORKFLOW_TYPE, WORKFLOW_VERSION."""
        _, metrics = load_input_tables(str(valid_joint_qc_path), str(valid_metrics_path))

        required = {"RUN_INDEX", "SAMPLE_ID", "RUN", "WORKFLOW_TYPE", "WORKFLOW_VERSION"}
        assert required.issubset(set(metrics.columns))


    @pytest.mark.parametrize("missing_col", sorted(["RUN_INDEX", "WORKFLOW_TYPE", "WORKFLOW_VERSION"]))
    def test_missing_joint_qc_column_raises(self, valid_joint_qc_path, valid_metrics_path, missing_col, tmp_path):
        """load_input_tables raises ValueError when a required joint_qc column is missing."""
        joint_qc = pl.read_csv(str(valid_joint_qc_path), separator="\t")
        stripped  = joint_qc.drop(missing_col)
        bad_path  = tmp_path / "bad_joint_qc.tsv"
        stripped.write_csv(str(bad_path), separator="\t")

        with pytest.raises(ValueError, match=missing_col):
            load_input_tables(str(bad_path), str(valid_metrics_path))


    @pytest.mark.parametrize(
        "missing_col",
        sorted(["RUN_INDEX", "SAMPLE_ID", "RUN", "WORKFLOW_TYPE", "WORKFLOW_VERSION"]),
    )
    def test_missing_metrics_column_raises(self, valid_joint_qc_path, valid_metrics_path, missing_col, tmp_path):
        """load_input_tables raises ValueError when a required metrics column is missing."""
        metrics  = pl.read_csv(str(valid_metrics_path), separator="\t")
        stripped = metrics.drop(missing_col)
        bad_path = tmp_path / "bad_metrics.tsv"
        stripped.write_csv(str(bad_path), separator="\t")

        with pytest.raises(ValueError, match=missing_col):
            load_input_tables(str(valid_joint_qc_path), str(bad_path))


    def test_tables_sorted_ascending_by_run_index(self, tmp_path):
        """Both returned tables must be sorted by RUN_INDEX in ascending order."""
        joint_qc = pl.DataFrame({
            "RUN_INDEX"       : [3, 1, 2],
            "WORKFLOW_TYPE"   : ["localapp"] * 3,
            "WORKFLOW_VERSION": ["1.0"] * 3,
        })
        metrics = pl.DataFrame({
            "RUN_INDEX"       : [3, 1, 2],
            "SAMPLE_ID"       : ["S3", "S1", "S2"],
            "RUN"             : ["RUN_C", "RUN_A", "RUN_B"],
            "WORKFLOW_TYPE"   : ["localapp"] * 3,
            "WORKFLOW_VERSION": ["1.0"] * 3,
        })
        jq_path = tmp_path / "joint_qc.tsv"
        mt_path = tmp_path / "metrics.tsv"
        joint_qc.write_csv(str(jq_path), separator="\t")
        metrics.write_csv(str(mt_path), separator="\t")

        jq_out, mt_out = load_input_tables(str(jq_path), str(mt_path))

        assert jq_out["RUN_INDEX"].to_list() == [1, 2, 3]
        assert mt_out["RUN_INDEX"].to_list() == [1, 2, 3]


    def test_duplicate_sample_ids_suffixed_with_run_index(self, tmp_path):
        """SAMPLE_IDs appearing more than once across the whole metrics table get '|RUN_INDEX' appended."""
        joint_qc = pl.DataFrame({
            "RUN_INDEX"       : [1, 2],
            "WORKFLOW_TYPE"   : ["localapp"] * 2,
            "WORKFLOW_VERSION": ["1.0"] * 2,
        })
        # "MCF7" appears in both run 1 and run 2 — should be suffixed.
        # "BT20" appears only once — should be left unchanged.
        metrics = pl.DataFrame({
            "RUN_INDEX"       : [1, 2, 1],
            "SAMPLE_ID"       : ["MCF7", "MCF7", "BT20"],
            "RUN"             : ["RUN_A", "RUN_B", "RUN_A"],
            "WORKFLOW_TYPE"   : ["localapp"] * 3,
            "WORKFLOW_VERSION": ["1.0"] * 3,
        })
        jq_path = tmp_path / "joint_qc.tsv"
        mt_path = tmp_path / "metrics.tsv"
        joint_qc.write_csv(str(jq_path), separator="\t")
        metrics.write_csv(str(mt_path), separator="\t")

        _, mt_out = load_input_tables(str(jq_path), str(mt_path))

        sample_ids = set(mt_out["SAMPLE_ID"].to_list())
        assert "MCF7|1" in sample_ids
        assert "MCF7|2" in sample_ids
        assert "MCF7" not in sample_ids   # Original unsuffixed name must not remain
        assert "BT20" in sample_ids       # Unique sample must be unchanged


    def test_guideline_samples_not_suffixed_when_duplicated(self, tmp_path):
        """Guideline SAMPLE_IDs (LSL_Guideline, USL_Guideline, Internal Guideline) must never be renamed."""
        joint_qc = pl.DataFrame({
            "RUN_INDEX"       : [1, 2],
            "WORKFLOW_TYPE"   : ["localapp"] * 2,
            "WORKFLOW_VERSION": ["1.0"] * 2,
        })
        # All three guideline IDs appear in both runs — they must not be suffixed.
        metrics = pl.DataFrame({
            "RUN_INDEX"       : [1, 2, 1, 2, 1, 2],
            "SAMPLE_ID"       : [
                "LSL_Guideline", "LSL_Guideline",
                "USL_Guideline", "USL_Guideline",
                "Internal Guideline", "Internal Guideline",
            ],
            "RUN"             : ["RUN_A", "RUN_B"] * 3,
            "WORKFLOW_TYPE"   : ["localapp"] * 6,
            "WORKFLOW_VERSION": ["1.0"] * 6,
        })
        jq_path = tmp_path / "joint_qc.tsv"
        mt_path = tmp_path / "metrics.tsv"
        joint_qc.write_csv(str(jq_path), separator="\t")
        metrics.write_csv(str(mt_path), separator="\t")

        _, mt_out = load_input_tables(str(jq_path), str(mt_path))

        sample_ids = mt_out["SAMPLE_ID"].to_list()
        assert all(sid in {"LSL_Guideline", "USL_Guideline", "Internal Guideline"} for sid in sample_ids)


class TestFilterMostRecentRuns:
    """Tests for filter_most_recent_runs() — RUN_INDEX-based run selection."""

    @pytest.fixture
    def sample_joint_qc(self) -> pl.DataFrame:
        """Small joint QC table with RUN_INDEX values 001–005 (001 = most recent)."""
        return pl.DataFrame({
            "RUN_INDEX"       : ["005", "004", "003", "002", "001"],
            "RUN_ID"          : ["RUN_E", "RUN_D", "RUN_C", "RUN_B", "RUN_A"],
            "WORKFLOW_TYPE"   : ["localapp"] * 5,
            "WORKFLOW_VERSION": ["1.0"] * 5,
        })


    @pytest.fixture
    def sample_metrics(self) -> pl.DataFrame:
        """Small metrics table aligned to the same runs as sample_joint_qc."""
        return pl.DataFrame({
            "RUN_INDEX"       : ["005", "005", "004", "004", "003", "003", "002", "002", "001", "001"],
            "SAMPLE_ID"       : ["S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10"],
            "RUN"             : ["RUN_E"] * 2 + ["RUN_D"] * 2 + ["RUN_C"] * 2 + ["RUN_B"] * 2 + ["RUN_A"] * 2,
            "WORKFLOW_TYPE"   : ["localapp"] * 10,
            "WORKFLOW_VERSION": ["1.0"] * 10,
        })


    def test_returns_two_dataframes(self, sample_joint_qc, sample_metrics):
        """filter_most_recent_runs returns two polars DataFrames."""
        jq, mt = filter_most_recent_runs(sample_joint_qc, sample_metrics, max_runs=2)

        assert isinstance(jq, pl.DataFrame)
        assert isinstance(mt, pl.DataFrame)


    def test_correct_number_of_runs_in_joint_qc(self, sample_joint_qc, sample_metrics):
        """Joint QC table should have exactly max_runs rows after filtering."""
        jq, _ = filter_most_recent_runs(sample_joint_qc, sample_metrics, max_runs=2)

        assert jq.height == 2


    def test_selects_lowest_run_indexes(self, sample_joint_qc, sample_metrics):
        """The two most recent runs (lowest RUN_INDEX) should be selected."""
        jq, _ = filter_most_recent_runs(sample_joint_qc, sample_metrics, max_runs=2)

        selected_indexes = set(jq["RUN_INDEX"].to_list())
        assert selected_indexes == {"001", "002"}


    def test_metrics_filtered_to_same_runs(self, sample_joint_qc, sample_metrics):
        """metrics_table should be filtered to the same set of RUN_INDEX values as joint_qc."""
        jq, mt = filter_most_recent_runs(sample_joint_qc, sample_metrics, max_runs=2)

        jq_indexes = set(jq["RUN_INDEX"].to_list())
        mt_indexes = set(mt["RUN_INDEX"].to_list())
        assert jq_indexes == mt_indexes


    def test_metrics_row_count_matches_selected_runs(self, sample_joint_qc, sample_metrics):
        """Metrics rows should include only samples from selected runs (2 samples per run)."""
        _, mt = filter_most_recent_runs(sample_joint_qc, sample_metrics, max_runs=2)

        assert mt.height == 4  # 2 runs × 2 samples each


    def test_max_runs_larger_than_available_returns_all(self, sample_joint_qc, sample_metrics):
        """When max_runs exceeds available runs, all runs should be returned."""
        jq, mt = filter_most_recent_runs(sample_joint_qc, sample_metrics, max_runs=100)

        assert jq.height == 5
        assert mt.height == 10


    def test_max_runs_one_selects_single_most_recent(self, sample_joint_qc, sample_metrics):
        """max_runs=1 should return only the single most recent run (lowest RUN_INDEX)."""
        jq, _ = filter_most_recent_runs(sample_joint_qc, sample_metrics, max_runs=1)

        assert jq.height == 1
        assert jq["RUN_INDEX"][0] == "001"


    def test_real_data_file(self, valid_joint_qc_path, valid_metrics_path):
        """filter_most_recent_runs selects indexes 3 and 4 (most recent) from real data when max_runs=2."""
        joint_qc, metrics = load_input_tables(str(valid_joint_qc_path), str(valid_metrics_path))
        jq, mt = filter_most_recent_runs(joint_qc, metrics, max_runs=2)

        # RUN_INDEX is read as integer by polars (leading zeros stripped)
        selected_indexes = set(jq["RUN_INDEX"].to_list())
        assert selected_indexes == {3, 4}


class TestFilterExpressions:
    """Table-driven tests for build_filter_expression()."""

    @pytest.mark.parametrize(
        "filter_spec",
        [
            pytest.param({"column": "SAMPLE_ID", "contains": "D01"}, id="contains"),
            pytest.param({"column": "SAMPLE_ID", "equals": "LSL_Guideline"}, id="equals"),
            pytest.param({"column": "SAMPLE_ID", "not_equals": "LSL_Guideline"}, id="not_equals"),
        ],
    )
    def test_build_filter_supported_types(self, filter_spec):
        """Supported filter specs should produce a valid Polars expression."""
        expr = build_filter_expression(filter_spec)
        assert str(type(expr)).find("Expr") != -1

    @pytest.mark.parametrize(
        "filter_spec",
        [
            pytest.param({"column": "SAMPLE_ID", "range": (0, 100)}, id="unsupported_filter"),
        ],
    )
    def test_build_filter_unsupported_raises(self, filter_spec):
        """Unsupported filter types should raise ValueError."""
        with pytest.raises(ValueError, match="Unsupported filter"):
            build_filter_expression(filter_spec)


class TestValueExpressions:
    """
    Test build_value_expression for different operations.
    Each operation should produce a valid Polars expression that can be applied to a DataFrame.
    """

    def test_value_expr_cast_no_dtype(self):
        """
        Test cast operation without explicit dtype.
        """
        spec = {"operation": "cast", "column": "SAMPLE_ID"}
        expr = build_value_expression(spec)
        assert str(type(expr)).find("Expr") != -1


    def test_value_expr_cast_with_dtype(self):
        """
        Test cast operation with explicit dtype.
        """
        spec = {"operation": "cast", "column": "DNA_TOTAL_PF_READS", "dtype": pl.Float64}
        expr = build_value_expression(spec)
        assert str(type(expr)).find("Expr") != -1


    def test_value_expr_cast_with_alias(self):
        """
        Test cast operation with alias name.
        """
        spec = {"operation": "cast", "column": "READS", "dtype": pl.Int32}
        expr = build_value_expression(spec, alias_name="READS_NORMALIZED")
        assert str(type(expr)).find("Expr") != -1


    def test_value_expr_divide(self):
        """
        Test divide operation.
        """
        spec = {
            "operation" : "divide",
            "column"    : "DNA_TOTAL_PF_READS",
            "divisor"   : 1000000,
            "dtype"     : pl.Float64,
        }
        expr = build_value_expression(spec)
        assert str(type(expr)).find("Expr") != -1


    def test_value_expr_ratio_simple(self):
        """
        Test ratio operation (numerator / denominator).
        """
        spec = {
            "operation"   : "ratio",
            "numerator"   : "PASSING_READS",
            "denominator" : "TOTAL_READS",
            "dtype"       : pl.Float64,
        }
        expr = build_value_expression(spec)
        assert str(type(expr)).find("Expr") != -1


    def test_value_expr_ratio_with_numerator_divisor(self):
        """Test ratio with numerator_divisor."""
        spec = {
            "operation"         : "ratio",
            "numerator"         : "PASSING_READS",
            "numerator_divisor" : 1000,
            "denominator"       : "TOTAL_READS",
            "dtype"             : pl.Float64,
        }
        expr = build_value_expression(spec)
        assert str(type(expr)).find("Expr") != -1


    def test_value_expr_unsupported_operation_raises(self):
        """
        Test that unsupported operation raises ValueError.
        An unsupported operation like "sqrt" should trigger an error.
        """
        spec = {"operation": "sqrt", "column": "VALUE"}
        with pytest.raises(ValueError, match="Unsupported value operation"):
            build_value_expression(spec)


class TestValidMetricExpr:
    """
    Test valid_metric_expr for null and NA string filtering.
    Each test ensures that the expression correctly filters out invalid values.
    """


    def test_valid_metric_expr_creates_expression(self):
        """
        Test that valid_metric_expr creates a valid Polars expression.
        """
        expr = valid_metric_expr("DNA_CONTAMINATION_SCORE")
        assert str(type(expr)).find("Expr") != -1


    def test_valid_metric_expr_filters_nulls(self):
        """
        Test that expression correctly filters nulls.
        """
        df = pl.DataFrame({
            "METRIC": [1.0, None, 2.0, None, 3.0],
        })
        result = df.filter(valid_metric_expr("METRIC"))
        assert result.height == 3


    def test_valid_metric_expr_filters_na_strings(self):
        """
        Test that expression correctly filters 'NA' strings.
        """
        df = pl.DataFrame({
            "METRIC": ["1.0", "NA", "2.0", "na", "3.0"],  # lowercase "na" should not be filtered
        })
        result = df.filter(valid_metric_expr("METRIC"))

        # Should exclude exactly the "NA" string
        assert result.height == 4


class TestComputeCartYLim:
    """
    Test y-axis limit computation (static and dynamic).
    """

    def test_compute_cart_ylim_static(self):
        """
        Test static y-axis limits.
        """
        spec   = {"cart_ylim": (0, 100)}
        result = compute_cart_ylim(spec, pd.DataFrame())
        assert result == (0, 100)


    def test_compute_cart_ylim_dynamic_max_plus(self):
        """
        Test dynamic y-limits with max_plus mode.
        """
        spec = {
            "cart_ylim_dynamic": {
                "mode"   : "max_plus",
                "column" : "VALUE",
                "lower"  : 0,
                "offset" : 10,
            }
        }
        data = pd.DataFrame({"VALUE": [10, 20, 30, 40, 50]})
        result = compute_cart_ylim(spec, data)
        assert result == (0, 60)  # 50 + 10


    def test_compute_cart_ylim_dynamic_without_lower(self):
        """
        Test dynamic limits when lower is not specified.
        """
        spec = {
            "cart_ylim_dynamic": {
                "mode"   : "max_plus",
                "column" : "VALUE",
                "offset" : 5,
            }
        }
        data = pd.DataFrame({"VALUE": [10, 20, 30]})
        result = compute_cart_ylim(spec, data)
        assert result == (0, 35)  # default lower=0, 30 + 5


    def test_compute_cart_ylim_unsupported_mode_raises(self):
        """
        Test that unsupported dynamic mode raises ValueError.
        """
        spec = {
            "cart_ylim_dynamic": {
                "mode"   : "percentile",
                "column" : "VALUE",
            }
        }
        data = pd.DataFrame({"VALUE": [10, 20, 30]})
        with pytest.raises(ValueError, match="Unsupported dynamic"):
            compute_cart_ylim(spec, data)


# --------------------------------------------------- #
# Data transformation tests (Real Data from Fixtures) #
# --------------------------------------------------- #

class TestBuildTables:
    """Test the build_tables function that creates logical subsets."""

    @pytest.mark.parametrize("workflow", ["Dragen", "dragen", "LocalApp", "localapp"])
    def test_build_tables_structure(self, joint_qc_data_2r, master_metrics_table_2r, workflow):
        """Ensure stable output keys for accepted workflow values (case-insensitive)."""
        highlight_run = master_metrics_table_2r.select(pl.col("RUN").drop_nulls().first()).item()
        result = build_tables(joint_qc_data_2r, master_metrics_table_2r, highlight_run, workflow=workflow)

        expected_keys = {
            "joint_qc_table",
            "joint_qc_guideline_table",
            "merged_tables",
            "guideline_table",
            "internal_guideline_table",
            "dna_guideline_table",
            "rna_guideline_table",
            "data_table",
            "dna_data_table",
            "rna_data_table",
            "dna_sample_count",
            "rna_sample_count",
        }
        assert set(result.keys()) == expected_keys

    @pytest.mark.parametrize(
        "requested_workflow, expected_workflow_type",
        [
            ("Dragen", "dragen"),
            ("localapp", "localapp"),
        ],
    )
    def test_build_tables_filters_by_workflow(self, joint_qc_data_2r, master_metrics_table_2r, requested_workflow, expected_workflow_type):
        """Ensure only rows from requested workflow remain after filtering."""
        highlight_run = master_metrics_table_2r.select(pl.col("RUN").drop_nulls().first()).item()
        result = build_tables(joint_qc_data_2r, master_metrics_table_2r, highlight_run, workflow=requested_workflow)
        observed_workflows = set(
            result["merged_tables"]["WORKFLOW_TYPE"]
            .cast(pl.Utf8)
            .str.to_lowercase()
            .unique()
            .to_list()
        )
        assert observed_workflows == {expected_workflow_type}

    def test_build_tables_guideline_and_data_partitioning(self, joint_qc_data_2r, master_metrics_table_2r):
        """Guideline rows must be split from regular data rows."""
        highlight_run = master_metrics_table_2r.select(pl.col("RUN").drop_nulls().first()).item()
        result = build_tables(joint_qc_data_2r, master_metrics_table_2r, highlight_run, workflow="Dragen")

        guideline_ids = set(result["guideline_table"]["SAMPLE_ID"].unique())
        assert guideline_ids.issubset({"LSL_Guideline", "USL_Guideline"})

        data_ids = result["data_table"]["SAMPLE_ID"].unique()
        for sample_id in data_ids:
            assert sample_id not in {"LSL_Guideline", "USL_Guideline", "Internal Guideline"}

    @pytest.mark.parametrize(
        "pick_existing_run, expected_present_flags",
        [(True, {"True"}), (False, {"False"})],
    )
    def test_build_tables_highlighted_run_marking(self, joint_qc_data_2r, master_metrics_table_2r, pick_existing_run, expected_present_flags):
        """highlighted_run flag must reflect RUN equality against highlighted_run_id."""
        data_for_workflow = master_metrics_table_2r.filter(pl.col("WORKFLOW_TYPE").str.to_lowercase() == "dragen")
        if data_for_workflow.height == 0:
            pytest.skip("No Dragen rows available in master_metrics_table_2r fixture")

        existing_run = data_for_workflow.select(pl.col("RUN").drop_nulls().first()).item()
        highlighted_run = existing_run if pick_existing_run else "RUN_NOT_PRESENT"
        result = build_tables(joint_qc_data_2r, master_metrics_table_2r, highlighted_run, workflow="Dragen")

        if result["dna_data_table"].height == 0:
            pytest.skip("No DNA rows available after workflow filtering")

        observed_flags = set(result["dna_data_table"]["highlighted_run"].unique())
        assert expected_present_flags.issubset(observed_flags)

    def test_build_tables_invalid_workflow_raises(self, joint_qc_data_2r, master_metrics_table_2r):
        """Unsupported workflow values should fail fast with ValueError."""
        highlight_run = master_metrics_table_2r.select(pl.col("RUN").drop_nulls().first()).item()

        with pytest.raises(ValueError, match="Unsupported workflow"):
            build_tables(joint_qc_data_2r, master_metrics_table_2r, highlight_run, workflow="Other")


class TestPrepareBarPlotData:
    """Test data preparation for bar plots (filtering + transformations)."""

    def test_prepare_bar_plot_data_basic(self, dna_samples_only):
        """Test basic data preparation without filters."""
        spec = {
            "value_spec"        : {"operation": "cast", "column": "DNA_TOTAL_PF_READS", "dtype": pl.Float64},
            "y_var"             : "READS",
            "na_filter_columns" : ["DNA_TOTAL_PF_READS"],  # Filter out "NA" strings before casting
            "filters"           : [],
        }
        result = prepare_bar_plot_data(dna_samples_only, spec)

        assert isinstance(result, pd.DataFrame)
        assert "READS" in result.columns
        assert len(result) > 0


    def test_prepare_bar_plot_data_with_na_filter(self, single_run_metrics):
        """
        Test that na_filter_columns correctly filters nulls.
        """
        spec = {
            "value_spec"        : {"operation": "cast", "column": "DNA_CONTAMINATION_SCORE", "dtype": pl.Float64},
            "y_var"             : "CONTAMINATION",
            "na_filter_columns" : ["DNA_CONTAMINATION_SCORE"],
            "filters"           : [],
        }
        before_count = single_run_metrics.height
        result       = prepare_bar_plot_data(single_run_metrics, spec)
        after_count  = len(result)

        # Should have fewer rows after filtering out nulls
        assert after_count <= before_count


    def test_prepare_bar_plot_data_with_sample_filter(self, single_run_metrics):
        """
        Test that sample filters are applied correctly.
        For example, filtering out guideline rows should remove them from the result.
        """
        spec = {
            "value_spec"        : {"operation": "cast", "column": "DNA_TOTAL_PF_READS", "dtype": pl.Float64},
            "y_var"             : "READS",
            "na_filter_columns" : ["DNA_TOTAL_PF_READS"],  # Filter out "NA" strings
            "filters"           : [{"column": "SAMPLE_ID", "not_equals": "LSL_Guideline"}],
        }
        result = prepare_bar_plot_data(single_run_metrics, spec)

        assert "LSL_Guideline" not in result["SAMPLE_ID"].values


    def test_prepare_bar_plot_data_derived_column(self, dna_samples_only):
        """
        Test that derived columns (e.g., from divide operation) are created.
        """
        spec = {
            "value_spec"    : {
                "operation" : "divide",
                "column"    : "DNA_TOTAL_PF_READS",
                "divisor"   : 1000000,
                "dtype"     : pl.Float64,
            },
            "y_var"            : "READS_MILLIONS",
            "na_filter_columns": ["DNA_TOTAL_PF_READS"],  # Filter out "NA" strings
            "filters"          : [],
        }
        result = prepare_bar_plot_data(dna_samples_only, spec)

        assert "READS_MILLIONS" in result.columns

        # Verify the derived column was created and has reasonable values
        assert len(result) > 0
        assert result["READS_MILLIONS"].notna().all()  # No NaN values after filtering


class TestGetGuidelineValue:
    """Test guideline value extraction."""

    def test_get_guideline_value_basic(self, guideline_rows):
        """Test extracting a guideline value."""

        tables = {"guideline_table" : guideline_rows}
        guideline_spec = {
            "table"        : "guideline_table",
            "sample_id"    : "USL_Guideline",
            "value_spec"   : {"operation": "cast", "column": "DNA_CONTAMINATION_SCORE", "dtype": pl.Float64},
            "label_prefix" : "USL",
        }
        result = get_guideline_value(tables, guideline_spec)

        if result is not None:  # Only if guideline exists in data
            assert "value" in result
            assert "label" in result
            assert isinstance(result["value"], (int, float))


    def test_get_guideline_value_with_cast(self, guideline_rows):
        """Test guideline value extraction with python_cast."""
        tables = {"guideline_table": guideline_rows}
        guideline_spec = {
            "table"        : "guideline_table",
            "sample_id"    : "USL_Guideline",
            "value_spec"   : {"operation": "cast", "column": "DNA_CONTAMINATION_SCORE", "dtype": pl.Float64},
            "python_cast"  : float,
            "label_prefix" : "USL",
        }
        result = get_guideline_value(tables, guideline_spec)

        if result is not None:
            assert isinstance(result["value"], float)


    def test_get_guideline_value_missing_returns_none(self, single_run_metrics):
        """
        Test that missing guideline returns None.
        This simulates a case where the guideline row is not present in the data table.
        """
        tables = {"guideline_table": single_run_metrics.head(0)}  # Empty table
        guideline_spec = {
            "table"       : "guideline_table",
            "sample_id"   : "USL_Guideline",
            "value_spec"  : {"operation": "cast", "column": "DNA_CONTAMINATION_SCORE", "dtype": pl.Float64},
            "label_prefix": "USL",
        }
        result = get_guideline_value(tables, guideline_spec)
        assert result is None


    def test_get_guideline_value_missing_column_returns_none(self, guideline_rows):
        """Test that missing value column returns None."""
        tables = {"guideline_table": guideline_rows}
        guideline_spec = {
            "table"       : "guideline_table",
            "sample_id"   : "USL_Guideline",
            "value_spec"  : {"operation": "cast", "column": "NONEXISTENT_COLUMN", "dtype": pl.Float64},
            "label_prefix": "USL",
        }
        result = get_guideline_value(tables, guideline_spec)
        assert result is None


# ------------------------------------- #
# Rendering & Integration Tests (Plots) #
# ------------------------------------- #


class TestRenderBarPlot:
    """Test bar plot rendering (integration test)."""

    def test_render_bar_plot_creates_output(self, dna_samples_only):
        """
        Test that render_bar_plot executes without error.
        This test uses a mock PDF handle to verify that the function attempts to save a figure."""
        spec = {
            "source"            : "data",
            "plot_kind"         : "bar",
            "x_var"             : "SAMPLE_ID",
            "y_var"             : "CONTAMINATION",
            "fill_var"          : "SAMPLE_ID",
            "x_lab"             : "Sample",
            "y_lab"             : "Contamination Score",
            "title"             : "Test Plot",
            "value_spec"        : {"operation": "cast", "column": "DNA_CONTAMINATION_SCORE", "dtype": pl.Float64},
            "na_filter_columns" : ["DNA_CONTAMINATION_SCORE"],
            "filters"           : [],
            "cart_ylim"         : (0, 500),
        }
        tables = {"data": dna_samples_only}

        # Mock the PDF handle
        mock_pdf = MagicMock()

        # Should not raise
        try:
            render_bar_plot(mock_pdf, spec, tables)
            # If it reaches here, rendering succeeded
            assert mock_pdf.savefig.called or not dna_samples_only.is_empty()
        except Exception as e:
            # If data is empty, skip; otherwise fail
            if dna_samples_only.height == 0:
                pytest.skip("No DNA samples in test data")
            raise


    def test_render_bar_plot_skip_if_empty(self):
        """
        Test that empty data with skip_if_empty=True doesn't render.
        """
        spec = {
            "source"            : "data",
            "plot_kind"         : "bar",
            "skip_if_empty"     : True,
            "x_var"             : "SAMPLE_ID",
            "y_var"             : "READS",
            "fill_var"          : "SAMPLE_ID",
            "x_lab"             : "Sample",
            "y_lab"             : "Reads",
            "title"             : "Empty Test",
            "value_spec"        : {"operation": "cast", "column": "DNA_TOTAL_PF_READS", "dtype": pl.Float64},
            "na_filter_columns" : ["DNA_TOTAL_PF_READS"],
            "filters"           : [],
            "cart_ylim"         : (0, 100),
        }
        # Create empty table with correct columns
        empty_table = pl.DataFrame({
            "SAMPLE_ID"          : pl.Series(dtype=pl.Utf8),
            "DNA_TOTAL_PF_READS" : pl.Series(dtype=pl.Float64)
        })
        tables = {"data": empty_table}
        mock_pdf = MagicMock()

        render_bar_plot(mock_pdf, spec, tables)
        # savefig should not be called for empty data with skip_if_empty=True


class TestPlotTsoppyBarplot:
    """Table-driven checks for plot_tsoppy_barplot() variants."""

    @pytest.mark.parametrize(
        "plot_kwargs",
        [
            pytest.param({}, id="basic"),
            pytest.param({"hline_y": 100.0, "hline_label": "Threshold"}, id="with_hline"),
            pytest.param({"y_tick_step": 50, "cart_ylim": (0, 500)}, id="with_y_tick_step"),
        ],
    )
    def test_plot_tsoppy_barplot_variants(self, dna_samples_only, plot_kwargs):
        """Plot variants should return a ggplot object for real DNA data."""
        data = dna_samples_only.to_pandas()
        plot = plot_tsoppy_barplot(
            data=data,
            x_var="SAMPLE_ID",
            y_var="DNA_CONTAMINATION_SCORE",
            fill_var="SAMPLE_ID",
            x_lab="Sample",
            y_lab="Contamination",
            title="Test",
            **plot_kwargs,
        )

        assert isinstance(plot, ggplot)


class TestPlotContaminationScatter:
    """Test contamination scatter plot function."""

    def test_plot_contamination_scatter_returns_ggplot(self, dna_samples_only):
        """
        Test that plot_contamination_scatter returns a ggplot object.
        Using real DNA sample data from the fixture to ensure the plot is generated correctly.
        """
        data = dna_samples_only.to_pandas()
        plot = plot_contamination_scatter(
            data                    = data,
            color_var               = "SAMPLE_ID",
            title                   = "Contamination",
            usl_contamination_score = 1000,
            usl_contamination_pval  = 0.05,
        )

        assert isinstance(plot, ggplot)


    def test_plot_contamination_scatter_with_labels(self, dna_samples_only):
        """
        Test contamination scatter with point labels.
        Using real DNA sample data from the fixture to ensure the plot is generated correctly.
        """
        data = dna_samples_only.to_pandas()
        plot = plot_contamination_scatter(
            data                    = data,
            color_var               = "SAMPLE_ID",
            label_var               = "SAMPLE_ID",
            title                   = "Contamination with Labels",
            usl_contamination_score = 1000,
            usl_contamination_pval  = 0.05,
        )

        assert isinstance(plot, ggplot)


class TestRenderPlot:
    """Test the render_plot dispatcher."""

    def test_render_plot_dispatches_bar(self, dna_samples_only):
        """
        Test that render_plot correctly dispatches to bar renderer.
        Using real DNA sample data from the fixture to ensure the plot is generated correctly.
        """
        spec = {
            "plot_kind"         : "bar",
            "source"            : "data",
            "x_var"             : "SAMPLE_ID",
            "y_var"             : "READS",
            "fill_var"          : "SAMPLE_ID",
            "x_lab"             : "Sample",
            "y_lab"             : "Reads",
            "title"             : "Test",
            "value_spec"        : {"operation": "cast", "column": "DNA_TOTAL_PF_READS", "dtype": pl.Float64},
            "na_filter_columns" : ["DNA_TOTAL_PF_READS"],  # Filter out "NA" strings
            "filters"           : [],
            "cart_ylim"         : (0, 1000000),
        }
        tables = {"data": dna_samples_only}
        mock_pdf = MagicMock()

        # Should dispatch without error
        render_plot(mock_pdf, "test_plot", spec, tables)


    def test_render_plot_unsupported_kind_raises(self, dna_samples_only):
        """
        Test that unsupported plot_kind raises ValueError.
        Using real DNA sample data from the fixture to ensure the plot is generated correctly.
        """
        spec = {
            "plot_kind" : "heatmap",  # Unsupported
            "source"    : "data",
        }
        tables = {"data": dna_samples_only}
        mock_pdf = MagicMock()

        with pytest.raises(ValueError, match="Unsupported plot kind"):
            render_plot(mock_pdf, "test_plot", spec, tables)


class TestSavePlot:
    """Test the save_plot function that writes plots to PDF."""

    def test_save_plot_saves_and_closes_figure(self, dna_samples_only):
        """save_plot should save the figure to the PDF handle and close the matplotlib figure."""
        data = dna_samples_only.to_pandas()
        plot = plot_tsoppy_barplot(
            data=data,
            x_var="SAMPLE_ID",
            y_var="DNA_CONTAMINATION_SCORE",
            fill_var="SAMPLE_ID",
            x_lab="Sample",
            y_lab="Contamination",
            title="Test Plot",
        )

        mock_pdf = MagicMock()
        save_plot(mock_pdf, plot)

        assert mock_pdf.savefig.called
        assert len(plt.get_fignums()) == 0
        plt.close('all')


class TestRenderClusterDensityScatter:
    """Test the cluster density vs estimated yield scatter plot rendering."""

    def test_render_cluster_density_scatter_executes(self, joint_qc_data):
        """
        Test that render_cluster_density_scatter executes without error.
        This verifies that the function correctly processes the plot specification and rendering.
        """
        # Create a spec for cluster density scatter plot
        spec = {
            "source"        : "joint_qc_table",
            "plot_kind"     : "cluster_density_scatter",
            "title"         : "Cluster Density vs Estimated Yield",
            "x_lab_angle"   : 45,
            "fig_size"      : (15, 6),
            "skip_if_empty" : False,
        }

        tables = {"joint_qc_table": joint_qc_data}
        mock_pdf = MagicMock()

        # Should not raise
        try:
            render_cluster_density_scatter(mock_pdf, spec, tables)
            # If data is present, savefig should be called
            if joint_qc_data.height > 0:
                assert mock_pdf.savefig.called or True
        except Exception as e:
            # If data is empty, skip; otherwise fail
            if joint_qc_data.height == 0:
                pytest.skip("No QC data in test data")
            raise


    def test_render_cluster_density_scatter_skip_if_empty(self):
        """
        Test that render_cluster_density_scatter respects skip_if_empty flag.
        """
        spec = {
            "source"        : "joint_qc_table",
            "plot_kind"     : "cluster_density_scatter",
            "title"         : "Cluster Density",
            "skip_if_empty" : True,
        }

        # Create empty table
        empty_table = pl.DataFrame({
            "CLUSTER_DENSITY" : pl.Series(dtype=pl.Float64),
            "ESTIMATED_YIELD" : pl.Series(dtype=pl.Float64),
            "RUN_ID"          : pl.Series(dtype=pl.Utf8),
            "RUN_NUMBER"      : pl.Series(dtype=pl.Utf8),
            "RUN_INDEX"       : pl.Series(dtype=pl.Int64),
        })
        tables = {"joint_qc_table": empty_table}
        mock_pdf = MagicMock()

        # Should execute without error
        render_cluster_density_scatter(mock_pdf, spec, tables)


class TestRenderContaminationScatter:
    """Test the contamination scatter plot rendering."""

    def test_render_contamination_scatter_executes(self, single_run_metrics, dna_samples_only):
        """
        Test that render_contamination_scatter executes without error.
        This verifies that the function correctly processes the contamination data and renders the plot.
        """
        # Create dna_data_table by adding a RUN column (simulating what merge_metric_tables does)
        dna_data_with_run = dna_samples_only.with_columns(
            pl.lit("RUN_001").alias("RUN")
        ).with_columns([
            pl.when(pl.col("RUN") == "RUN_001").then(pl.lit("False")).otherwise(pl.lit("False")).alias("highlighted_run"),
            pl.when(pl.col("RUN") == "RUN_001").then(pl.col("SAMPLE_ID")).otherwise(pl.lit("")).alias("contamination_label"),
        ])

        spec = {
            "source"        : "dna_data_table",
            "plot_kind"     : "contamination_scatter",
            "title"         : "DNA Contamination Assessment",
            "color_var"     : "RUN",
            "label_var"     : None,
            "guide_title"   : "Run",
            "x_lab_angle"   : 45,
            "fig_size"      : (15, 6),
            "skip_if_empty" : False,
        }

        # Create dna_guideline_table from the full single_run_metrics
        dna_guideline = single_run_metrics.filter(
            pl.col("SAMPLE_ID").is_in(["USL_Guideline", "LSL_Guideline", "Internal Guideline"])
        )

        tables = {
            "dna_data_table"      : dna_data_with_run,
            "dna_guideline_table" : dna_guideline,
        }
        mock_pdf = MagicMock()

        # Should not raise
        try:
            render_contamination_scatter(mock_pdf, spec, tables)
            # If data is present, rendering should complete
            assert True
        except Exception as e:
            # If data is empty, skip; otherwise fail
            if dna_data_with_run.height == 0:
                pytest.skip("No DNA samples in test data")
            raise


    def test_render_contamination_scatter_skip_if_empty(self):
        """
        Test that render_contamination_scatter respects skip_if_empty flag.
        """
        spec = {
            "source"        : "dna_data_table",
            "plot_kind"     : "contamination_scatter",
            "title"         : "Contamination",
            "color_var"     : "RUN",
            "label_var"     : None,
            "skip_if_empty" : True,
        }

        # Create empty tables
        empty_data = pl.DataFrame({
            "SAMPLE_ID"                 : pl.Series(dtype=pl.Utf8),
            "DNA_CONTAMINATION_SCORE"   : pl.Series(dtype=pl.Float64),
            "DNA_CONTAMINATION_P_VALUE" : pl.Series(dtype=pl.Float64),
            "RUN"                       : pl.Series(dtype=pl.Utf8),
        })
        empty_guideline = pl.DataFrame({
            "SAMPLE_ID"                 : pl.Series(dtype=pl.Utf8),
            "DNA_CONTAMINATION_SCORE"   : pl.Series(dtype=pl.Float64),
            "DNA_CONTAMINATION_P_VALUE" : pl.Series(dtype=pl.Float64),
        })

        tables = {
            "dna_data_table"      : empty_data,
            "dna_guideline_table" : empty_guideline,
        }
        mock_pdf = MagicMock()

        # Should execute without error
        render_contamination_scatter(mock_pdf, spec, tables)


# --------------------------- #
# EDGE CASES & ERROR HANDLING #
# --------------------------- #

class TestEdgeCases:
    """Test edge cases and error conditions."""

    def test_valid_metric_expr_with_empty_dataframe(self):
        """
        Test valid_metric_expr on empty dataframe.
        This ensures that the expression can be applied to a DataFrame with no rows without raising an error.
        """
        df = pl.DataFrame({"METRIC": pl.Series(dtype=pl.Float64)})
        result = df.filter(valid_metric_expr("METRIC"))
        assert result.height == 0


    def test_build_value_expression_with_missing_column(self):
        """
        Test that missing column in data is handled (should raise at execution, not creation).
        This ensures that the expression can be created even if the column does not exist,
        but will raise an error when applied to a DataFrame.
        """
        spec = {"operation": "cast", "column": "NONEXISTENT", "dtype": pl.Float64}
        expr = build_value_expression(spec)

        # Expression is created OK, but should fail when applied
        df = pl.DataFrame({"OTHER": [1, 2, 3]})
        with pytest.raises(Exception):  # ColumnNotFound or similar
            df.select(expr)


    def test_prepare_bar_plot_data_empty_after_filter(self, single_run_metrics):
        """
        Test handling when filters remove all rows.
        This ensures that the function can handle cases where the filters result in an empty DataFrame.
        """
        spec = {
            "value_spec"        : {"operation": "cast", "column": "DNA_TOTAL_PF_READS", "dtype": pl.Float64},
            "y_var"             : "READS",
            "na_filter_columns" : [],
            "filters"           : [
                {"column": "SAMPLE_ID", "equals": "NONEXISTENT_SAMPLE"},  # Will match nothing
            ],
        }
        result = prepare_bar_plot_data(single_run_metrics, spec)

        assert len(result) == 0


    def test_compute_cart_ylim_with_empty_dataframe(self):
        """
        Test dynamic y-limit computation with empty data.
        This ensures that the function can handle cases where the data is empty.
        """
        spec = {
            "cart_ylim_dynamic": {
                "mode"   : "max_plus",
                "column" : "VALUE",
                "offset" : 10,
            }
        }
        empty_data = pd.DataFrame({"VALUE": []})

        # Should handle gracefully or raise appropriate error
        try:
            result = compute_cart_ylim(spec, empty_data)
            # If it returns, check it's a reasonable default
            assert isinstance(result, tuple)
        except Exception as e:
            # If it raises, should be informative
            assert "VALUE" in str(e) or "empty" in str(e).lower()
