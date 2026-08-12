# ========================================================= #
# python script to create LocalApp and Dragen metrics plots #
# ========================================================= #
RUN_ID=250822_A02134_0101_BHYFHLDRX5


# ================= #
# Example 1 : 1 run #
# ================= #
RESULTS_DIR_python=results_python/TSO500_1R
mkdir -p ${RESULTS_DIR_python}
echo "; Running plot_QC_metrics: 1 run"
uv run plot_QC_metrics.py --joint_QC data/merged_workflow_tables/joint_sequencing_QC_file_1R.tsv  \
                          --output ${RESULTS_DIR_python}                                          \
                          --metrics_table data/merged_workflow_tables/master_metrics_table_1R.tsv \
                          --run_id ${RUN_ID}     


# ================== #
# Example 2 : 2 runs #
# ================== #
mkdir -p ${RESULTS_DIR_python}
echo "; Running plot_QC_metrics: 2 runs"
uv run plot_QC_metrics.py --joint_QC data/merged_workflow_tables/joint_sequencing_QC_file_2R.tsv  \
                          --output ${RESULTS_DIR_python}                                          \
                          --metrics_table data/merged_workflow_tables/master_metrics_table_2R.tsv \
                          --run_id ${RUN_ID}                                                      ;


# ================== #
# Example 3 : 5 runs #
# ================== #
RESULTS_DIR_python=results_python/TSO500_5R
mkdir -p ${RESULTS_DIR_python}
echo "; Running plot_QC_metrics: 5 runs"
uv run plot_QC_metrics.py --joint_QC data/merged_workflow_tables/joint_sequencing_QC_file_5R.tsv  \
                          --output ${RESULTS_DIR_python}                                          \
                          --metrics_table data/merged_workflow_tables/master_metrics_table_5R.tsv \
                          --run_id ${RUN_ID}                                                      ;


# =================== #
# Example 4 : 10 runs #
# =================== #
RESULTS_DIR_python=results_python/TSO500_10R
mkdir -p ${RESULTS_DIR_python}
echo "; Running plot_QC_metrics: 10 runs"
uv run plot_QC_metrics.py --joint_QC data/merged_workflow_tables/joint_sequencing_QC_file_10R.tsv  \
                          --output ${RESULTS_DIR_python}                                          \
                          --metrics_table data/merged_workflow_tables/master_metrics_table_10R.tsv \
                          --run_id ${RUN_ID}       



# ======================================= #
# Restore the original plot_specs.py file #
# ======================================= #
cp specs/plot_specs_minimal.py specs/plot_specs.py

# TSVs="${METRICS} ${METRICS_PREV1}"
# RUNIDS="${RUN_ID} ${RUN_ID_PREV1}"

# uv run plot_run_metrics.py D 2 ${METRIC_PLOTS_DIR}/intermediate_metrics_files/joint_sequencing_QC_file.tsv ${RESULTS_DIR_python}/TSO500_run_metrics.pdf ${METRIC_PLOTS_DIR}/intermediate_metrics_files/master_metrics_table.tsv ${RUN_ID} False ${TSVs} ${RUNIDS}
# uv run plot_run_metrics.py D 2 ${METRIC_PLOTS_DIR}/intermediate_metrics_files/joint_sequencing_QC_file.tsv ${RESULTS_DIR_python}/TSO500_run_metrics.pdf ${METRIC_PLOTS_DIR}/intermediate_metrics_files/master_metrics_table.tsv ${RUN_ID} True ${TSVs} ${RUNIDS}