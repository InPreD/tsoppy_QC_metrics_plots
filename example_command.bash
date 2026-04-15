# Example command
# Rscript plot_run_metrics.R    \
#     1                         \ 
#     /inpred/data//sample_data/results/metrics_plots_HUS_all/intermediate_metrics_files/joint_sequencing_QC_file.tsv \
#     /inpred/data//sample_data/results/metrics_plots_HUS_all/TSO500_run_metrics.pdf                                  \
#     /inpred/data//sample_data/results/metrics_plots_HUS_all/intermediate_metrics_files/master_metrics_table.tsv     \
#     260312_NDX551102_RUO_0042_AH3HTFBGY2                                                                            \
#     False                                                                                                           \
#     /inpred/data//sample_data/results/metrics_plots_HUS_all/intermediate_metrics_files/260312_NDX551102_RUO_0042_AH3HTFBGY2_metrics.tsv    \
#     260312_NDX551102_RUO_0042_AH3HTFBGY2                                                                             ;



RUN_ID=251211_NB552811_0117_AHTJM2BGYX
TSOPPI_DIR=data/TSO500_postprocessing
METRIC_PLOTS_DIR=${TSOPPI_DIR}/${RUN_ID}/${RUN_ID}_metrics_plots
RESULTS_DIR_R=results_R/${RUN_ID}
RESULTS_DIR_python=results_python/${RUN_ID}

# ------------------------------------ #
# R script to create the metrics plots #
# ------------------------------------ #
echo "; Running R script to create the metrics plots"
mkdir -p ${RESULTS_DIR_R}

Rscript plot_run_metrics.R 1 ${METRIC_PLOTS_DIR}/intermediate_metrics_files/joint_sequencing_QC_file.tsv ${RESULTS_DIR_R}/TSO500_run_metrics.pdf ${METRIC_PLOTS_DIR}/intermediate_metrics_files/master_metrics_table.tsv ${RUN_ID} False ${METRIC_PLOTS_DIR}/intermediate_metrics_files/${RUN_ID}_metrics.tsv ${RUN_ID}
Rscript plot_run_metrics.R 1 ${METRIC_PLOTS_DIR}/intermediate_metrics_files/joint_sequencing_QC_file.tsv ${RESULTS_DIR_R}/TSO500_run_metrics.pdf ${METRIC_PLOTS_DIR}/intermediate_metrics_files/master_metrics_table.tsv ${RUN_ID} True ${METRIC_PLOTS_DIR}/intermediate_metrics_files/${RUN_ID}_metrics.tsv ${RUN_ID}


# ----------------------------------------- #
# python script to create the metrics plots #
# ----------------------------------------- #
echo "; Running python script to create the metrics plots"
mkdir -p ${RESULTS_DIR_python}

uv run plot_run_metrics.py 1 ${METRIC_PLOTS_DIR}/intermediate_metrics_files/joint_sequencing_QC_file.tsv ${RESULTS_DIR_R}/TSO500_run_metrics.pdf ${METRIC_PLOTS_DIR}/intermediate_metrics_files/master_metrics_table.tsv ${RUN_ID} False ${METRIC_PLOTS_DIR}/intermediate_metrics_files/${RUN_ID}_metrics.tsv ${RUN_ID}
uv run plot_run_metrics.py 1 ${METRIC_PLOTS_DIR}/intermediate_metrics_files/joint_sequencing_QC_file.tsv ${RESULTS_DIR_R}/TSO500_run_metrics.pdf ${METRIC_PLOTS_DIR}/intermediate_metrics_files/master_metrics_table.tsv ${RUN_ID} True ${METRIC_PLOTS_DIR}/intermediate_metrics_files/${RUN_ID}_metrics.tsv ${RUN_ID}