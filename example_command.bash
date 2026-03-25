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
RESULTS_DIR=results/${RUN_ID}

mkdir -p ${RESULTS_DIR}

Rscript plot_run_metrics.R 1 ${METRIC_PLOTS_DIR}/intermediate_metrics_files/joint_sequencing_QC_file.tsv ${RESULTS_DIR}/TSO500_run_metrics.pdf ${METRIC_PLOTS_DIR}/intermediate_metrics_files/master_metrics_table.tsv ${RUN_ID} False ${METRIC_PLOTS_DIR}/intermediate_metrics_files/${RUN_ID}_metrics.tsv ${RUN_ID}
Rscript plot_run_metrics.R 1 ${METRIC_PLOTS_DIR}/intermediate_metrics_files/joint_sequencing_QC_file.tsv ${RESULTS_DIR}/TSO500_run_metrics.pdf ${METRIC_PLOTS_DIR}/intermediate_metrics_files/master_metrics_table.tsv ${RUN_ID} True ${METRIC_PLOTS_DIR}/intermediate_metrics_files/${RUN_ID}_metrics.tsv ${RUN_ID}