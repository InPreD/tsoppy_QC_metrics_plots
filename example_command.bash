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


# ------------------------------------- #
# RUN ID to be highlighted in the plots #
# ------------------------------------- #
RUN_ID=251211_NB552811_0117_AHTJM2BGYX
TSOPPI_DIR=data/TSO500_postprocessing
METRIC_PLOTS_DIR=${TSOPPI_DIR}/${RUN_ID}/${RUN_ID}_metrics_plots
RESULTS_DIR_R=results_R/${RUN_ID}
RESULTS_DIR_python=results_python/${RUN_ID}
METRICS=${METRIC_PLOTS_DIR}/intermediate_metrics_files/${RUN_ID}_metrics.tsv


# ----------------------------------- #
# Previous RUN_ID, add more as needed #
# ----------------------------------- #
RUN_ID_PREV1=260521_NB552811_0117_AHTJM2BGYX
METRICS_PREV1=${METRIC_PLOTS_DIR}/intermediate_metrics_files/${RUN_ID_PREV1}_metrics.tsv

RUN_ID_PREV2=260520_NB552811_0117_AHTJM2BGYX
METRICS_PREV2=${METRIC_PLOTS_DIR}/intermediate_metrics_files/${RUN_ID_PREV2}_metrics.tsv

RUN_ID_PREV3=260519_NB552811_0117_AHTJM2BGYX
METRICS_PREV3=${METRIC_PLOTS_DIR}/intermediate_metrics_files/${RUN_ID_PREV3}_metrics.tsv

RUN_ID_PREV4=260518_NB552811_0117_AHTJM2BGYX
METRICS_PREV4=${METRIC_PLOTS_DIR}/intermediate_metrics_files/${RUN_ID_PREV4}_metrics.tsv

RUN_ID_PREV5=260517_NB552811_0117_AHTJM2BGYX
METRICS_PREV5=${METRIC_PLOTS_DIR}/intermediate_metrics_files/${RUN_ID_PREV5}_metrics.tsv

RUN_ID_PREV6=260516_NB552811_0117_AHTJM2BGYX
METRICS_PREV6=${METRIC_PLOTS_DIR}/intermediate_metrics_files/${RUN_ID_PREV6}_metrics.tsv

RUN_ID_PREV7=260515_NB552811_0117_AHTJM2BGYX
METRICS_PREV7=${METRIC_PLOTS_DIR}/intermediate_metrics_files/${RUN_ID_PREV7}_metrics.tsv

RUN_ID_PREV8=260514_NB552811_0117_AHTJM2BGYX
METRICS_PREV8=${METRIC_PLOTS_DIR}/intermediate_metrics_files/${RUN_ID_PREV8}_metrics.tsv

RUN_ID_PREV9=260513_NB552811_0117_AHTJM2BGYX
METRICS_PREV9=${METRIC_PLOTS_DIR}/intermediate_metrics_files/${RUN_ID_PREV9}_metrics.tsv


# --------------------------------- #
# Strings with all TSVs and RUN_IDs #
# --------------------------------- #
# To be passed to the plotting scripts
# Add more as needed
TSVs="${METRICS} ${METRICS_PREV1} ${METRICS_PREV2} ${METRICS_PREV3} ${METRICS_PREV4} ${METRICS_PREV5} ${METRICS_PREV6} ${METRICS_PREV7} ${METRICS_PREV8} ${METRICS_PREV9}"
RUNIDS="${RUN_ID} ${RUN_ID_PREV1} ${RUN_ID_PREV2} ${RUN_ID_PREV3} ${RUN_ID_PREV4} ${RUN_ID_PREV5} ${RUN_ID_PREV6} ${RUN_ID_PREV7} ${RUN_ID_PREV8} ${RUN_ID_PREV9}"

# # Example with 1 run only
# TSVs="${METRICS}"
# RUNIDS="${RUN_ID}"

# # Example with 2 runs
# TSVs="${METRICS} ${METRICS_PREV1}"
# RUNIDS="${RUN_ID} ${RUN_ID_PREV1}"


# ------------------------------------ #
# R script to create the metrics plots #
# ------------------------------------ #
echo "; Running R script to create the metrics plots"
mkdir -p ${RESULTS_DIR_R}

Rscript plot_run_metrics.R 10 ${METRIC_PLOTS_DIR}/intermediate_metrics_files/joint_sequencing_QC_file.tsv ${RESULTS_DIR_R}/TSO500_run_metrics.pdf ${METRIC_PLOTS_DIR}/intermediate_metrics_files/master_metrics_table.tsv ${RUN_ID} False ${TSVs} ${RUNIDS}
Rscript plot_run_metrics.R 10 ${METRIC_PLOTS_DIR}/intermediate_metrics_files/joint_sequencing_QC_file.tsv ${RESULTS_DIR_R}/TSO500_run_metrics.pdf ${METRIC_PLOTS_DIR}/intermediate_metrics_files/master_metrics_table.tsv ${RUN_ID} True ${TSVs} ${RUNIDS}


# ----------------------------------------- #
# python script to create the metrics plots #
# ----------------------------------------- #
echo "; Running python script to create the metrics plots"
mkdir -p ${RESULTS_DIR_python}

uv run plot_run_metrics.py 10 ${METRIC_PLOTS_DIR}/intermediate_metrics_files/joint_sequencing_QC_file.tsv ${RESULTS_DIR_python}/TSO500_run_metrics.pdf ${METRIC_PLOTS_DIR}/intermediate_metrics_files/master_metrics_table.tsv ${RUN_ID} False ${TSVs} ${RUNIDS}
uv run plot_run_metrics.py 10 ${METRIC_PLOTS_DIR}/intermediate_metrics_files/joint_sequencing_QC_file.tsv ${RESULTS_DIR_python}/TSO500_run_metrics.pdf ${METRIC_PLOTS_DIR}/intermediate_metrics_files/master_metrics_table.tsv ${RUN_ID} True ${TSVs} ${RUNIDS}

#Rscript    /inpred/user_scripts/libs/plot_run_metrics.R    1    /inpred/data//sample_data/results/metrics_plots_HUS_latest/intermediate_metrics_files/joint_sequencing_QC_file.tsv    /inpred/data//sample_data/results/metrics_plots_HUS_latest/TSO500_run_metrics.pdf    /inpred/data//sample_data/results/metrics_plots_HUS_latest/intermediate_metrics_files/master_metrics_table.tsv    260312_NDX551102_RUO_0042_AH3HTFBGY2    True    /inpred/data//sample_data/results/metrics_plots_HUS_latest/intermediate_metrics_files/260312_NDX551102_RUO_0042_AH3HTFBGY2_metrics.tsv    260312_NDX551102_RUO_0042_AH3HTFBGY2