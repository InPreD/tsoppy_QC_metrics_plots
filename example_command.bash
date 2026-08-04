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
TSOPPI_DIR=data/LocalApp/TSO500_postprocessing
METRIC_PLOTS_DIR=${TSOPPI_DIR}/${RUN_ID}/${RUN_ID}_metrics_plots
RESULTS_DIR_R=results_R/LocalApp/${RUN_ID}
RESULTS_DIR_python=results_python/LocalApp/${RUN_ID}
METRICS=${METRIC_PLOTS_DIR}/intermediate_metrics_files/${RUN_ID}_metrics.tsv


# ----------------------------------- #
# Previous RUN_ID, add more as needed #
# ----------------------------------- #
RUN_ID_PREV1=251210_NB552811_0117_AHTJM2BGYX
METRICS_PREV1=${METRIC_PLOTS_DIR}/intermediate_metrics_files/${RUN_ID_PREV1}_metrics.tsv

RUN_ID_PREV2=251209_NB552811_0117_AHTJM2BGYX
METRICS_PREV2=${METRIC_PLOTS_DIR}/intermediate_metrics_files/${RUN_ID_PREV2}_metrics.tsv

RUN_ID_PREV3=251208_NB552811_0117_AHTJM2BGYX
METRICS_PREV3=${METRIC_PLOTS_DIR}/intermediate_metrics_files/${RUN_ID_PREV3}_metrics.tsv

RUN_ID_PREV4=251207_NB552811_0117_AHTJM2BGYX
METRICS_PREV4=${METRIC_PLOTS_DIR}/intermediate_metrics_files/${RUN_ID_PREV4}_metrics.tsv

RUN_ID_PREV5=251206_NB552811_0117_AHTJM2BGYX
METRICS_PREV5=${METRIC_PLOTS_DIR}/intermediate_metrics_files/${RUN_ID_PREV5}_metrics.tsv

RUN_ID_PREV6=251205_NB552811_0117_AHTJM2BGYX
METRICS_PREV6=${METRIC_PLOTS_DIR}/intermediate_metrics_files/${RUN_ID_PREV6}_metrics.tsv

RUN_ID_PREV7=251204_NB552811_0117_AHTJM2BGYX
METRICS_PREV7=${METRIC_PLOTS_DIR}/intermediate_metrics_files/${RUN_ID_PREV7}_metrics.tsv

RUN_ID_PREV8=251203_NB552811_0117_AHTJM2BGYX
METRICS_PREV8=${METRIC_PLOTS_DIR}/intermediate_metrics_files/${RUN_ID_PREV8}_metrics.tsv

RUN_ID_PREV9=251202_NB552811_0117_AHTJM2BGYX
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


# --------------------------------------------- #
# R script to create the LocalApp metrics plots #
# --------------------------------------------- #
echo "; Running R script to create the LocalApp metrics plots"
mkdir -p ${RESULTS_DIR_R}

Rscript plot_run_metrics.R 10 ${METRIC_PLOTS_DIR}/intermediate_metrics_files/joint_sequencing_QC_file.tsv ${RESULTS_DIR_R}/TSO500_run_metrics.pdf ${METRIC_PLOTS_DIR}/intermediate_metrics_files/master_metrics_table.tsv ${RUN_ID} False ${TSVs} ${RUNIDS}
Rscript plot_run_metrics.R 10 ${METRIC_PLOTS_DIR}/intermediate_metrics_files/joint_sequencing_QC_file.tsv ${RESULTS_DIR_R}/TSO500_run_metrics.pdf ${METRIC_PLOTS_DIR}/intermediate_metrics_files/master_metrics_table.tsv ${RUN_ID} True ${TSVs} ${RUNIDS}



# ------------------------------------------------------------------- #
# python script to create the plots using a minimal set of parameters #
# ------------------------------------------------------------------- #
RESULTS_DIR_python=results_python/Minimal/${RUN_ID}
echo "; Running python script to create the metrics plots using minimal set of parameters"
mkdir -p ${RESULTS_DIR_python}

# The first time the script is run, it will create the joint_sequencing_QC_file.tsv file, which will be used in subsequent runs
# The second time it will use the joint_sequencing_QC_file.tsv file to create the plots
cp specs/plot_specs_minimal.py specs/plot_specs.py
uv run plot_run_metrics.py 10 ${METRIC_PLOTS_DIR}/intermediate_metrics_files/joint_sequencing_QC_file.tsv ${RESULTS_DIR_python}/TSO500_run_metrics.pdf ${METRIC_PLOTS_DIR}/intermediate_metrics_files/master_metrics_table.tsv ${RUN_ID} False ${TSVs} ${RUNIDS}
uv run plot_run_metrics.py 10 ${METRIC_PLOTS_DIR}/intermediate_metrics_files/joint_sequencing_QC_file.tsv ${RESULTS_DIR_python}/TSO500_run_metrics.pdf ${METRIC_PLOTS_DIR}/intermediate_metrics_files/master_metrics_table.tsv ${RUN_ID} True ${TSVs} ${RUNIDS}



# -------------------------------------------------- #
# python script to create the LocalApp metrics plots #
# -------------------------------------------------- #
echo "; Running python script to create the LocalApp metrics plots"
mkdir -p ${RESULTS_DIR_python}

# The first time the script is run, it will create the joint_sequencing_QC_file.tsv file, which will be used in subsequent runs
# The second time it will use the joint_sequencing_QC_file.tsv file to create the plots
cp specs/plot_specs_LocalApp.py specs/plot_specs.py
uv run plot_run_metrics.py 10 ${METRIC_PLOTS_DIR}/intermediate_metrics_files/joint_sequencing_QC_file.tsv ${RESULTS_DIR_python}/TSO500_run_metrics.pdf ${METRIC_PLOTS_DIR}/intermediate_metrics_files/master_metrics_table.tsv ${RUN_ID} False ${TSVs} ${RUNIDS}
uv run plot_run_metrics.py 10 ${METRIC_PLOTS_DIR}/intermediate_metrics_files/joint_sequencing_QC_file.tsv ${RESULTS_DIR_python}/TSO500_run_metrics.pdf ${METRIC_PLOTS_DIR}/intermediate_metrics_files/master_metrics_table.tsv ${RUN_ID} True ${TSVs} ${RUNIDS}



# ------------------------------------- #
# RUN ID to be highlighted in the plots #
# ------------------------------------- #
RUN_ID=250822_A02134_0101_BHYFHLDRX5
TSOPPI_DIR=data/DRAGEN/TSO500_postprocessing
METRIC_PLOTS_DIR=${TSOPPI_DIR}/${RUN_ID}/${RUN_ID}_metrics_plots
RESULTS_DIR_python=results_python/DRAGEN/${RUN_ID}
METRICS=${METRIC_PLOTS_DIR}/intermediate_metrics_files/${RUN_ID}_metrics.tsv

# ----------------------------------- #
# Previous RUN_ID, add more as needed #
# ----------------------------------- #
RUN_ID_PREV1=250627_A02134_0093_BHYHNYDRX5
METRICS_PREV1=${METRIC_PLOTS_DIR}/intermediate_metrics_files/${RUN_ID_PREV1}_metrics.tsv

RUN_ID_PREV2=250613_A02134_0091_BHYFHJDRX5
METRICS_PREV2=${METRIC_PLOTS_DIR}/intermediate_metrics_files/${RUN_ID_PREV2}_metrics.tsv

RUN_ID_PREV3=250530_A02134_0089_BHYGC3DRX5
METRICS_PREV3=${METRIC_PLOTS_DIR}/intermediate_metrics_files/${RUN_ID_PREV3}_metrics.tsv

RUN_ID_PREV4=250516_A02134_0087_BHLMFHDRX5
METRICS_PREV4=${METRIC_PLOTS_DIR}/intermediate_metrics_files/${RUN_ID_PREV4}_metrics.tsv

RUN_ID_PREV5=250425_A02134_0083_BHYHMWDRX5
METRICS_PREV5=${METRIC_PLOTS_DIR}/intermediate_metrics_files/${RUN_ID_PREV5}_metrics.tsv

RUN_ID_PREV6=250404_A02134_0079_BHTLC7DRX5
METRICS_PREV6=${METRIC_PLOTS_DIR}/intermediate_metrics_files/${RUN_ID_PREV6}_metrics.tsv

RUN_ID_PREV7=250321_A02134_0073_BHTMFHDRX5
METRICS_PREV7=${METRIC_PLOTS_DIR}/intermediate_metrics_files/${RUN_ID_PREV7}_metrics.tsv

RUN_ID_PREV8=241004_A02134_0026_BHJ2CFDRX5
METRICS_PREV8=${METRIC_PLOTS_DIR}/intermediate_metrics_files/${RUN_ID_PREV8}_metrics.tsv

RUN_ID_PREV9=240927_A02134_0024_BHHGG7DRX5
METRICS_PREV9=${METRIC_PLOTS_DIR}/intermediate_metrics_files/${RUN_ID_PREV9}_metrics.tsv

TSVs="${METRICS} ${METRICS_PREV1} ${METRICS_PREV2} ${METRICS_PREV3} ${METRICS_PREV4} ${METRICS_PREV5} ${METRICS_PREV6} ${METRICS_PREV7} ${METRICS_PREV8} ${METRICS_PREV9}"
RUNIDS="${RUN_ID} ${RUN_ID_PREV1} ${RUN_ID_PREV2} ${RUN_ID_PREV3} ${RUN_ID_PREV4} ${RUN_ID_PREV5} ${RUN_ID_PREV6} ${RUN_ID_PREV7} ${RUN_ID_PREV8} ${RUN_ID_PREV9}"


# ------------------------------------------------ #
# python script to create the DRAGEN metrics plots #
# ------------------------------------------------ #
echo "; Running python script to create the DRAGEN metrics plots"
mkdir -p ${RESULTS_DIR_python}

# In this example the joint_sequencing_QC_file.tsv file has already been created, so the script will run only to create the plots
cp specs/plot_specs_Dragen.py specs/plot_specs.py
uv run plot_run_metrics.py 10 ${METRIC_PLOTS_DIR}/intermediate_metrics_files/joint_sequencing_QC_file.tsv ${RESULTS_DIR_python}/TSO500_run_metrics.pdf ${METRIC_PLOTS_DIR}/intermediate_metrics_files/master_metrics_table.tsv ${RUN_ID} False ${TSVs} ${RUNIDS}
uv run plot_run_metrics.py 10 ${METRIC_PLOTS_DIR}/intermediate_metrics_files/joint_sequencing_QC_file.tsv ${RESULTS_DIR_python}/TSO500_run_metrics.pdf ${METRIC_PLOTS_DIR}/intermediate_metrics_files/master_metrics_table.tsv ${RUN_ID} True ${TSVs} ${RUNIDS}


# ======================================= #
# Restore the original plot_specs.py file #
# ======================================= #
cp specs/plot_specs_minimal.py specs/plot_specs.py

# TSVs="${METRICS} ${METRICS_PREV1}"
# RUNIDS="${RUN_ID} ${RUN_ID_PREV1}"

# uv run plot_run_metrics.py D 2 ${METRIC_PLOTS_DIR}/intermediate_metrics_files/joint_sequencing_QC_file.tsv ${RESULTS_DIR_python}/TSO500_run_metrics.pdf ${METRIC_PLOTS_DIR}/intermediate_metrics_files/master_metrics_table.tsv ${RUN_ID} False ${TSVs} ${RUNIDS}
# uv run plot_run_metrics.py D 2 ${METRIC_PLOTS_DIR}/intermediate_metrics_files/joint_sequencing_QC_file.tsv ${RESULTS_DIR_python}/TSO500_run_metrics.pdf ${METRIC_PLOTS_DIR}/intermediate_metrics_files/master_metrics_table.tsv ${RUN_ID} True ${TSVs} ${RUNIDS}