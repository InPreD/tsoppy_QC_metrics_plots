# ========================================================= #
# python script to create LocalApp and Dragen metrics plots #
# ========================================================= #

# Configuration variables
MAX_RUNS="1 2 5 10"
WORKFLOWS="localapp dragen"

# ========================================= #
# Process different run counts: 1, 2, 5, 10 #
# ========================================= #
for RUN_COUNT in ${MAX_RUNS}; do
    RESULTS_DIR=results_python/TSO500_${RUN_COUNT}R
    mkdir -p ${RESULTS_DIR}
    
    for WORKFLOW in ${WORKFLOWS}; do
        RESULTS_FILE=${RESULTS_DIR}/TSO500_${RUN_COUNT}R_${WORKFLOW}_metric_plots.pdf
        
        echo "; Running plot_QC_metrics: ${RUN_COUNT} run(s) (${WORKFLOW})"
        uv run plot_QC_metrics.py --joint_QC      data/workflow_tables/${WORKFLOW}_plot_joint_sequencing_QC_file.tsv \
                                  --output_pdf    ${RESULTS_FILE}                                                    \
                                  --metrics_table data/workflow_tables/${WORKFLOW}_plot_metrics_table.tsv            \
                                  --max_runs      ${RUN_COUNT}                                                       \
                                  --workflow_type ${WORKFLOW}                                                        ;
    done
done


