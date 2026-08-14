# ========================================================= #
# python script to create LocalApp and Dragen metrics plots #
# ========================================================= #
RUN_ID=250822_A02134_0101_BHYFHLDRX5

# Configuration variables
RUN_COUNTS="1R 2R 5R 10R"
WORKFLOWS="LocalApp Dragen"

# ======================================== #
# Process all run counts: 1R, 2R, 5R, 10R #
# ======================================== #
for RUN_COUNT in ${RUN_COUNTS}; do
    RESULTS_DIR=results_python/TSO500_${RUN_COUNT}
    mkdir -p ${RESULTS_DIR}
    
    for WORKFLOW in ${WORKFLOWS}; do
        RESULTS_FILE=${RESULTS_DIR}/TSO500_${RUN_COUNT}_${WORKFLOW}_metric_plots.pdf
        
        echo "; Running plot_QC_metrics: ${RUN_COUNT} run(s) (${WORKFLOW})"
        uv run plot_QC_metrics.py --joint_QC      data/merged_workflow_tables/joint_sequencing_QC_file_${RUN_COUNT}.tsv \
                                  --output_pdf    ${RESULTS_FILE}                                                       \
                                  --metrics_table data/merged_workflow_tables/master_metrics_table_${RUN_COUNT}.tsv     \
                                  --run_id        ${RUN_ID}                                                             \
                                  --workflow_type ${WORKFLOW}                                                           ;
    done
done


