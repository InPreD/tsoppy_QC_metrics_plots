version_string <- "0.3:22-03-10"

# assign the command-line arguments to script variables
args <- commandArgs(trailingOnly = TRUE)

# command-line parameters:
# - 1: N (number of runs)
# - 2: joint run QC file
# - 3: output PDF file with plotting results
# - 4: output tsv file with merged metrics data
# - 5: label of the run that should be highlighted in the plots
# - 5 + 1, .., (5 + N): transposed metrics matrices for the individual runs
# - (5 + N + 1) - (5 + N + N): labels for the individual runs
run_count <- as.numeric(args[1])
joint_qc_file <- args[2]
output_pdf <- args[3]
merged_metrics_file <- args[4]
highlighted_run_id <- args[5]
create_plots <- args[6]

# exit if zero runs are submitted
if (run_count < 1) {
  stop("Zero runs provided as input? Exiting.")
}

# load the necessary libraries
library("ggplot2")
library("cowplot")
library("plyr")
theme_set(theme_cowplot())

joint_qc_table <- read.table(joint_qc_file, header = T, sep = "\t")
merged_tables <- read.table(args[7], header = T, sep = "\t",
                            stringsAsFactors = FALSE)
merged_tables$RUN <- args[7 + run_count]

print("  Merging input tables..")

if (run_count > 1) {
  for (run_number in 2:run_count) {
    run_table <- read.table(args[6 + run_number], header = T, sep = "\t",
                            stringsAsFactors = FALSE)
    run_table$RUN <- args[6 + run_count + run_number]
    print(paste("    iteration", run_number))
    merged_tables <- unique(rbind(merged_tables, run_table))
  }
}

print("  Writing merged input tables into a master metrics file..")
write.table(unique(merged_tables), merged_metrics_file, quote = FALSE,
            row.names = FALSE, sep = "\t")

if (create_plots == "False") {
  print("  Skipping metrics plotting. Exiting.")
  quit("no", 0)
}

# keep the DNA guideline values, RNA guideline values,
# the RNA data and the DNA data separate from each other
# (they are to be used differently in the plotting)
guideline_table <- merged_tables[merged_tables$SAMPLE_ID == "LSL_Guideline"
    | merged_tables$SAMPLE_ID == "USL_Guideline", ]
# both LSL and USL guidelines are specified
# for the CONTAMINATION_SCORE DNA metric
dna_guideline_table <- guideline_table[
    !is.na(guideline_table$DNA_CONTAMINATION_SCORE), ]
# both LSL and USL guidelines are specified
# for	the MEDIAN_CV_GENE_500X RNA metric
rna_guideline_table <- guideline_table[
    !is.na(guideline_table$RNA_MEDIAN_CV_GENE_500X), ]

data_table <- merged_tables[merged_tables$SAMPLE_ID != "LSL_Guideline"
    & merged_tables$SAMPLE_ID != "USL_Guideline", ]
# each sample with specified DNA TOTAL_PF_READS count
# is included in the DNA data table
dna_data_table <- data_table[!is.na(data_table$DNA_TOTAL_PF_READS), ]
# each sample with specified RNA TOTAL_PF_READS count
# is included in the RNA data table
rna_data_table <- data_table[!is.na(data_table$RNA_TOTAL_PF_READS), ]

dna_sample_count <- nrow(dna_data_table)
rna_sample_count <- nrow(rna_data_table)

# samples from the highlighted run will be labeled
# on the contamination scatter plot
dna_data_table$highlighted_run <- "False"
dna_data_table$highlighted_run[
    dna_data_table$RUN == highlighted_run_id] <- "True"
dna_data_table$contamination_label <- dna_data_table$SAMPLE_ID
dna_data_table$contamination_label[
    dna_data_table$RUN != highlighted_run_id] <- ""

print("  Creating the metrics plots..")
pdf(output_pdf, 15, 6, onefile = TRUE)

###############################
# Sequencer-generated metrics #
###############################

ggplot(joint_qc_table[!is.na(joint_qc_table$CLUSTERS_PASSING_FILTER), ],
       aes(x = RUN_ID, y = CLUSTERS_PASSING_FILTER)) +
    geom_col(aes(fill = RUN_ID), alpha = 0.8) +
    guides(fill = guide_legend("Run")) +
    xlab("\nRun ID") +
    ylab("Clusters passing filter (percent)") +
    coord_cartesian(xlim = c(1, nrow(joint_qc_table) + 1), ylim = c(0, 104)) +
    ggtitle(paste("[Sequencer flowcell run metric] Clusters passing filter ",
                  "(based on RunCompletionStatus.xml files)\n",
                  "(only displaying runs with non-NA values)", sep = "")) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                     hjust = 1, size = 7)) +
    background_grid(major = "xy")

ggplot(joint_qc_table[!is.na(joint_qc_table$ESTIMATED_YIELD), ],
       aes(x = CLUSTER_DENSITY, y = ESTIMATED_YIELD)) +
    geom_point(aes(color = paste(RUN_NUMBER, " - ", RUN_ID, sep = "")),
               alpha = 0.8, size = 4) +
    geom_text(aes(label = RUN_NUMBER), alpha = 0.8, size = 4, nudge_y = -5) +
    guides(color = guide_legend("Run")) +
    coord_cartesian(xlim = c(0, 400), ylim = c(0, 150)) +
    xlab("\nCluster density") +
    ylab("Estimated yield") +
    ggtitle(paste("[Sequencer flowcell run metric] Cluster density against",
                  " estimated yield ",
                  "(based on RunCompletionStatus.xml files)\n",
                  "(only displaying runs with non-NA values)", sep = "")) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                     hjust = 1, size = 7)) +
    background_grid(major = "xy", minor = "xy") +
    scale_x_continuous(breaks = 50 * c(0:8)) +
    scale_y_continuous(breaks = 10 * c(0:15))


#############################
# LocalApp core run metrics #
#############################

ggplot(joint_qc_table, aes(x = RUN_ID, y = PCT_PF_READS)) +
    geom_col(aes(fill = RUN_ID), alpha = 0.8) +
    guides(fill = guide_legend("Run")) +
    xlab("\nRun ID") +
    ylab("Percentage of PF reads") +
    coord_cartesian(xlim = c(1, nrow(joint_qc_table) + 1), ylim = c(0, 104)) +
    ggtitle(paste("[LocalApp core run metric] Percentage of reads ",
                  "passing instrument filters", sep = "")) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                     hjust = 1, size = 7)) +
    geom_hline(yintercept = 80, alpha = 0.3, col = "red") +
    annotate("text", nrow(joint_qc_table) + 1, 80,
             label = paste("LSL_Guideline: 80", sep = ""), angle = 90) +
    background_grid(major = "xy")

ggplot(joint_qc_table, aes(x = RUN_ID, y = PCT_Q30_R1)) +
    geom_col(aes(fill = RUN_ID), alpha = 0.8) +
    guides(fill = guide_legend("Run")) +
    xlab("\nRun ID") +
    ylab("Percentage of Q30 R1 reads") +
    coord_cartesian(xlim = c(1, nrow(joint_qc_table) + 1), ylim = c(0, 104)) +
    ggtitle(paste("[LocalApp core run metric] Percentage of R1 reads",
                  " with quality ", ">=", " 30", sep = "")) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                     hjust = 1, size = 7)) +
    geom_hline(yintercept = 80, alpha = 0.3, col = "red") +
    annotate("text", nrow(joint_qc_table) + 1, 80,
             label = paste("LSL_Guideline: 80", sep = ""), angle = 90) +
    background_grid(major = "xy")

ggplot(joint_qc_table, aes(x = RUN_ID, y = PCT_Q30_R2)) +
    geom_col(aes(fill = RUN_ID), alpha = 0.8) +
    guides(fill = guide_legend("Run")) +
    xlab("\nRun ID") +
    ylab("Percentage of Q30 R2 reads") +
    coord_cartesian(xlim = c(1, nrow(joint_qc_table) + 1), ylim = c(0, 104)) +
    ggtitle(paste("[LocalApp core run metric] Percentage of R2 reads",
                  " with quality ", ">=", " 30", sep = "")) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                     hjust = 1, size = 7)) +
    geom_hline(yintercept = 80, alpha = 0.3, col = "red") +
    annotate("text", nrow(joint_qc_table) + 1, 80,
             label = paste("LSL_Guideline: 80", sep = ""), angle = 90) +
    background_grid(major = "xy")


#############################
# LocalApp core DNA metrics #
#############################

if (dna_sample_count > 0) {

  # [Illumina core DNA metric] Contamination assessment by TSO 500 LocalApp
  usl_contamination_score <-
      dna_guideline_table$DNA_CONTAMINATION_SCORE[
          dna_guideline_table$SAMPLE_ID == "USL_Guideline"][1]
  usl_contamination_pval <-
      dna_guideline_table$DNA_CONTAMINATION_P_VALUE[
          dna_guideline_table$SAMPLE_ID == "USL_Guideline"][1]
  max_contamination_score <- max(5000, dna_data_table$DNA_CONTAMINATION_SCORE)

  print(ggplot(dna_data_table, aes(x = DNA_CONTAMINATION_SCORE,
                             y = DNA_CONTAMINATION_P_VALUE)) +
    geom_point(aes(color = RUN)) +
    guides(color = guide_legend("Run")) +
    xlab("Contamination score") +
    ylab("Contamination P-value") +
    coord_cartesian(xlim = c(0, max_contamination_score + 500),
                    ylim = c(0.0, 1.05)) +
    scale_x_continuous(
        breaks = 1000 * c(0:(max_contamination_score / 1000) + 1)) +
    scale_y_continuous(breaks = 0.1 * c(0:10)) +
    ggtitle(paste("[LocalApp core DNA metric] Contamination assessment ",
                  "by TSO 500 LocalApp\n(samples located within the red box",
                  " are deemed contaminated)", sep = "")) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                     hjust = 1, size = 7)) +
    background_grid(major = "xy", minor = "xy") +
    annotate("rect", xmin = usl_contamination_score,
             xmax = max_contamination_score + 250,
             ymin = usl_contamination_pval,
             ymax = 1, fill = "red", alpha = 0.3) +
    geom_hline(yintercept = usl_contamination_pval, alpha = 0.3, col = "red") +
    annotate("text", max_contamination_score + 500,
             usl_contamination_pval + 0.1,
             label = paste("USL_Guideline: ", usl_contamination_pval, sep = ""),
             angle = 90) +
    geom_vline(xintercept = usl_contamination_score, alpha = 0.3, col = "red") +
    annotate("text", usl_contamination_score + 200, 1.05,
             label = paste("USL_Guideline: ", usl_contamination_score,
                           sep = "")))

  print(ggplot(dna_data_table, aes(x = DNA_CONTAMINATION_SCORE,
                             y = DNA_CONTAMINATION_P_VALUE)) +
    geom_point(aes(color = highlighted_run)) +
    guides(color = guide_legend("Latest run")) +
    xlab("Contamination score") +
    ylab("Contamination P-value") +
    coord_cartesian(xlim = c(0, max_contamination_score + 500),
                    ylim = c(0.0, 1.05)) +
    geom_text(aes(label = contamination_label), angle = -25, nudge_y = 0.02,
              size = 2) +
    scale_x_continuous(
        breaks = 1000 * c(0:(max_contamination_score / 1000) + 1)) +
    scale_y_continuous(breaks = 0.1 * c(0:10)) +
    ggtitle(paste("[LocalApp core DNA metric] Contamination assessment",
                  " by TSO 500 LocalApp\n(samples located within the red box",
                  " are deemed contaminated)", sep = "")) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                     hjust = 1, size = 7)) +
    background_grid(major = "xy", minor = "xy") +
    annotate("rect", xmin = usl_contamination_score,
             xmax = max_contamination_score + 250,
             ymin = usl_contamination_pval,
             ymax = 1, fill = "red", alpha = 0.3) +
    geom_hline(yintercept = usl_contamination_pval, alpha = 0.3, col = "red") +
    annotate("text", max_contamination_score + 500,
             usl_contamination_pval + 0.1,
             label = paste("USL_Guideline: ", usl_contamination_pval, sep = ""),
             angle = 90) +
    geom_vline(xintercept = usl_contamination_score, alpha = 0.3, col = "red") +
    annotate("text", usl_contamination_score + 200, 1.05,
             label = paste("USL_Guideline: ", usl_contamination_score,
                           sep = "")))

  # [Illumina core DNA metric] Median exon coverage
  lsl_dna_median_exon_coverage <-
      dna_guideline_table$DNA_MEDIAN_EXON_COVERAGE[
          dna_guideline_table$SAMPLE_ID == "LSL_Guideline"][1]

  print(ggplot(dna_data_table, aes(x = SAMPLE_ID,
                                   y = DNA_MEDIAN_EXON_COVERAGE)) +
    geom_col(aes(fill = RUN), alpha = 0.8) +
    coord_cartesian(xlim = c(1, dna_sample_count + 1)) +
    xlab("Sample ID") +
    ylab("Median exon coverage") +
    scale_y_continuous(breaks = 100 * c(0:20)) +
    ggtitle("[LocalApp core DNA metric] Median exon coverage") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                     hjust = 1, size = 7)) +
    geom_hline(yintercept = lsl_dna_median_exon_coverage,
               alpha = 0.3, col = "red") +
    annotate("text", dna_sample_count + 1,
             lsl_dna_median_exon_coverage + 100,
             label = paste("LSL_Guideline: ",
                           lsl_dna_median_exon_coverage, sep = ""),
             angle = 90) +
    guides(fill = guide_legend("Run")) +
    background_grid(major = "xy", minor = "y") +
    scale_x_discrete(limits = dna_data_table$SAMPLE_ID))

  # [Illumina core DNA metric] Percentage of exons with coverage >= 50X
  lsl_dna_pct_exon_50x <-
      dna_guideline_table$DNA_PCT_EXON_50X[
          dna_guideline_table$SAMPLE_ID == "LSL_Guideline"][1]

  print(ggplot(dna_data_table, aes(x = SAMPLE_ID, y = DNA_PCT_EXON_50X)) +
    geom_col(aes(fill = RUN), alpha = 0.8) +
    coord_cartesian(xlim = c(1, dna_sample_count + 1), , ylim = c(0, 104)) +
    xlab("Sample ID") +
    ylab(paste("Percentage of exons with coverage ", ">=", " 50X",
               sep = "")) +
    scale_y_continuous(breaks = 10 * c(0:10)) +
    ggtitle(paste("[LocalApp core DNA metric] Percentage of exons",
                  " with coverage ", ">=", " 50X", sep = "")) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                     hjust = 1, size = 7)) +
    geom_hline(yintercept = lsl_dna_pct_exon_50x, alpha = 0.3, col = "red") +
    annotate("text", dna_sample_count + 1, lsl_dna_pct_exon_50x - 10,
             label = paste("LSL_Guideline: ", lsl_dna_pct_exon_50x, sep = ""),
             angle = 90) +
    guides(fill = guide_legend("Run")) +
    background_grid(major = "xy", minor = "y") +
    scale_x_discrete(limits = dna_data_table$SAMPLE_ID))

  # [Illumina core DNA metric] Median insert size
  lsl_median_insert_size <-
      dna_guideline_table$DNA_MEDIAN_INSERT_SIZE[
          dna_guideline_table$SAMPLE_ID == "LSL_Guideline"][1]

  print(ggplot(dna_data_table, aes(x = SAMPLE_ID, y = DNA_MEDIAN_INSERT_SIZE)) +
    geom_col(aes(fill = RUN), alpha = 0.8) +
    coord_cartesian(xlim = c(1, dna_sample_count + 1)) +
    xlab("Sample ID") +
    ylab("Median insert size") +
    scale_y_continuous(breaks = 50 * c(0:20)) +
    ggtitle("[LocalApp core DNA metric] Median insert size") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                     hjust = 1, size = 7)) +
    geom_hline(yintercept = lsl_median_insert_size, alpha = 0.3, col = "red") +
    annotate("text", dna_sample_count + 1, lsl_median_insert_size,
             label = paste("LSL_Guideline: ", lsl_median_insert_size, sep = ""),
             angle = 90) +
    guides(fill = guide_legend("Run")) +
    background_grid(major = "xy", minor = "y") +
    scale_x_discrete(limits = dna_data_table$SAMPLE_ID))

  # [Illumina core DNA metric] Coverage MAD
  usl_median_coverage_mad <-
      dna_guideline_table$DNA_COVERAGE_MAD[
          dna_guideline_table$SAMPLE_ID == "USL_Guideline"][1]

  print(ggplot(dna_data_table, aes(x = SAMPLE_ID, y = DNA_COVERAGE_MAD)) +
    geom_col(aes(fill = RUN), alpha = 0.8) +
    coord_cartesian(xlim = c(1, dna_sample_count + 1)) +
    xlab("Sample ID") +
    ylab("Coverage MAD (Median Absolute Deviation)") +
    #scale_y_continuous(breaks = 50*c(0:20)) +
    ggtitle(paste("[LocalApp core DNA metric] Coverage MAD",
                  "\n(median normalized deviation across all regions",
                  " used for CNV calling)", sep = "")) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                     hjust = 1, size = 7)) +
    geom_hline(yintercept = usl_median_coverage_mad, alpha = 0.3, col = "red") +
    annotate("text", dna_sample_count + 1, usl_median_coverage_mad,
             label = paste("USL_Guideline: ", usl_median_coverage_mad,
                           sep = ""),
             angle = 90) +
    guides(fill = guide_legend("Run")) +
    background_grid(major = "xy", minor = "y") +
    scale_x_discrete(limits = dna_data_table$SAMPLE_ID))

  # [Illumina core DNA metric] Median raw bin count per CNV target
  lsl_median_bin_count <-
      dna_guideline_table$DNA_MEDIAN_BIN_COUNT_CNV_TARGET[
          dna_guideline_table$SAMPLE_ID == "LSL_Guideline"][1]

  print(ggplot(dna_data_table, aes(x = SAMPLE_ID,
                             y = DNA_MEDIAN_BIN_COUNT_CNV_TARGET)) +
    geom_col(aes(fill = RUN), alpha = 0.8) +
    coord_cartesian(xlim = c(1, dna_sample_count + 1)) +
    xlab("Sample ID") +
    ylab("Median raw bin count per CNV target") +
    #scale_y_continuous(breaks = 50*c(0:20)) +
    ggtitle(paste("[LocalApp core DNA metric] Median raw bin count",
                  " per CNV target", sep = "")) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                     hjust = 1, size = 7)) +
    geom_hline(yintercept = lsl_median_bin_count, alpha = 0.3, col = "red") +
    annotate("text", dna_sample_count + 1, lsl_median_bin_count + 1,
             label = paste("LSL_Guideline: ", lsl_median_bin_count, sep = ""),
             angle = 90) +
    guides(fill = guide_legend("Run")) +
    background_grid(major = "xy", minor = "y") +
    scale_x_discrete(limits = dna_data_table$SAMPLE_ID))

  # [Illumina core DNA metric] Number of usable MSI sites
  lsl_usable_msi_sites <-
      dna_guideline_table$DNA_USABLE_MSI_SITES[
          dna_guideline_table$SAMPLE_ID == "LSL_Guideline"][1]

  print(ggplot(dna_data_table, aes(x = SAMPLE_ID, y = DNA_USABLE_MSI_SITES)) +
    geom_col(aes(fill = RUN), alpha = 0.8) +
    coord_cartesian(xlim = c(1, dna_sample_count + 1)) +
    xlab("Sample ID") +
    ylab("Usable MSI sites") +
    scale_y_continuous(breaks = 50 * c(0:20)) +
    ggtitle("[LocalApp core DNA metric] Number of usable MSI sites") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                     hjust = 1, size = 7)) +
    geom_hline(yintercept = lsl_usable_msi_sites, alpha = 0.3, col = "red") +
    annotate("text", dna_sample_count + 1, lsl_usable_msi_sites,
             label = paste("LSL_Guideline: ", lsl_usable_msi_sites, sep = ""),
             angle = 90) +
    guides(fill = guide_legend("Run")) +
    background_grid(major = "xy", minor = "y") +
    scale_x_discrete(limits = dna_data_table$SAMPLE_ID))

} else {
  print("  Were no DNA samples included in the plotted runs?")
}

#############################
# LocalApp core RNA metrics #
#############################

if (rna_sample_count > 0) {

  # [Illumina core RNA metric] Uniformity of coverage across transcripts
  usl_median_cv <-
      rna_guideline_table$RNA_MEDIAN_CV_GENE_500X[
          rna_guideline_table$SAMPLE_ID == "USL_Guideline"][1]

  print(ggplot(rna_data_table, aes(x = SAMPLE_ID,
                                   y = RNA_MEDIAN_CV_GENE_500X)) +
    geom_col(aes(fill = RUN), alpha = 0.8) +
    coord_cartesian(xlim = c(1, rna_sample_count + 1), ylim = c(0, 125)) +
    xlab("Sample ID") +
    ylab("Median CV for genes with median coverage >500") +
    scale_y_continuous(breaks = 10 * c(0:12)) +
    ggtitle(paste("[LocalApp core RNA metric] Uniformity of coverage",
                  " across transcripts", sep = "")) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                     hjust = 1, size = 7)) +
    geom_hline(yintercept = usl_median_cv, alpha = 0.3, col = "red") +
    annotate("text", rna_sample_count + 1, usl_median_cv,
             label = paste("USL_Guideline: ", usl_median_cv, sep = ""),
             angle = 90) +
    guides(fill = guide_legend("Run")) +
    background_grid(major = "xy", minor = "y") +
    scale_x_discrete(limits = rna_data_table$SAMPLE_ID))

  # [Illumina core RNA metric] Total on-target reads in millions
  lsl_total_on_target_reads <-
      (rna_guideline_table$RNA_TOTAL_ON_TARGET_READS[
          rna_guideline_table$SAMPLE_ID == "LSL_Guideline"][1]) / 1000000
  max_total_on_target_reads <-
      max(rna_data_table$RNA_TOTAL_ON_TARGET_READS) + 10

  print(ggplot(rna_data_table, aes(x = SAMPLE_ID,
                                   y = RNA_TOTAL_ON_TARGET_READS / 1000000)) +
    geom_col(aes(fill = RUN), alpha = 0.8) +
    coord_cartesian(xlim = c(1, rna_sample_count + 1)) +
    xlab("Sample ID") +
    ylab("On target reads (millions)") +
    scale_y_continuous(breaks = 10 * c(0:(max_total_on_target_reads / 10))) +
    ggtitle("[LocalApp core RNA metric] Total on-target reads in millions") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                     hjust = 1, size = 7)) +
    geom_hline(yintercept = lsl_total_on_target_reads, alpha = 0.3,
               col = "red") +
    annotate("text", rna_sample_count + 1, lsl_total_on_target_reads,
             label = paste("LSL_Guideline: ", lsl_total_on_target_reads,
                           sep = ""),
             angle = 90) +
    guides(fill = guide_legend("Run")) +
    background_grid(major = "xy", minor = "y") +
    scale_x_discrete(limits = rna_data_table$SAMPLE_ID))

  # [Illumina core RNA metric] Median insert size
  lsl_median_insert_size <-
      (rna_guideline_table$RNA_MEDIAN_INSERT_SIZE[
          rna_guideline_table$SAMPLE_ID == "LSL_Guideline"][1])
  max_median_insert_size <- max(rna_data_table$RNA_MEDIAN_INSERT_SIZE) + 10

  print(ggplot(rna_data_table, aes(x = SAMPLE_ID, y = RNA_MEDIAN_INSERT_SIZE)) +
    geom_col(aes(fill = RUN), alpha = 0.8) +
    coord_cartesian(xlim = c(1, rna_sample_count + 1)) +
    xlab("Sample ID") +
    ylab("Median insert size") +
    scale_y_continuous(breaks = 10 * c(0:(max_median_insert_size / 10))) +
    ggtitle("[LocalApp core RNA metric] Median insert size") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                     hjust = 1, size = 7)) +
    geom_hline(yintercept = lsl_median_insert_size, alpha = 0.3, col = "red") +
    annotate("text", rna_sample_count + 1, lsl_median_insert_size,
             label = paste("LSL_Guideline: ", lsl_median_insert_size, sep = ""),
             angle = 90) +
    guides(fill = guide_legend("Run")) +
    background_grid(major = "xy", minor = "y") +
    scale_x_discrete(limits = rna_data_table$SAMPLE_ID))
} else {
  print("  No RNA samples included in the plotted runs?")
}

#########################
# Ih-house core metrics #
#########################

if (dna_sample_count > 0) {

  # [In-house core metric] Pipeline completion (DNA)
  print(ggplot(dna_data_table, aes(x = SAMPLE_ID, y = COMPLETED_ALL_STEPS)) +
    geom_col(aes(fill = RUN), alpha = 0.8) +
    guides(fill = guide_legend("Run")) +
    xlab("Sample ID") +
    ylab("Pipeline completion") +
    ggtitle("[InPreD core DNA metric] Pipeline completion") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                     hjust = 1, size = 7)) +
    background_grid(major = "xy") +
    scale_x_discrete(limits = dna_data_table$SAMPLE_ID))
}

if (rna_sample_count > 0) {

  # [In-house core metric] Pipeline completion (RNA)
  print(ggplot(rna_data_table, aes(x = SAMPLE_ID, y = COMPLETED_ALL_STEPS)) +
    geom_col(aes(fill = RUN), alpha = 0.8) +
    guides(fill = guide_legend("Run")) +
    xlab("Sample ID") +
    ylab("Pipeline completion") +
    ggtitle("[InPreD core RNA metric] Pipeline completion") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                     hjust = 1, size = 7)) +
    background_grid(major = "xy") +
    scale_x_discrete(limits = rna_data_table$SAMPLE_ID))
}

if (dna_sample_count > 0) {

  # [In-house core metric] Number of reads in millions (DNA)
  print(ggplot(dna_data_table, aes(x = SAMPLE_ID,
                                   y = DNA_TOTAL_PF_READS / 1000000)) +
    geom_col(aes(fill = RUN), alpha = 0.8) +
    coord_cartesian(xlim = c(1, dna_sample_count + 1)) +
    guides(fill = guide_legend("Run")) +
    xlab("Sample ID") +
    ylab("Number of reads (millions)") +
    scale_y_continuous(breaks = 20 * c(0:10)) +
    ggtitle("[InPreD core DNA metric] Number of reads in millions") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                     hjust = 1, size = 7)) +
    geom_hline(yintercept = 80, alpha = 0.3, col = "red") +
    annotate("text", dna_sample_count + 1, 80,
             label = paste("Internal Guideline: ", 80, sep = ""), angle = 90) +
    background_grid(major = "xy", minor = "y") +
    scale_x_discrete(limits = dna_data_table$SAMPLE_ID))

  # [In-house core metric] Mean target coverage (DNA)
  print(ggplot(dna_data_table, aes(x = SAMPLE_ID,
                                   y = DNA_MEAN_TARGET_COVERAGE)) +
    geom_col(aes(fill = RUN), alpha = 0.8) +
    coord_cartesian(xlim = c(1, dna_sample_count + 1)) +
    guides(fill = guide_legend("Run")) +
    xlab("Sample ID") +
    ylab("Mean target coverage") +
    scale_y_continuous(breaks = 100 * c(0:20)) +
    ggtitle("[InPreD core DNA metric] Mean target coverage") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                     hjust = 1, size = 7)) +
    geom_hline(yintercept = 350, alpha = 0.3, col = "red") +
    annotate("text", dna_sample_count + 1, 350,
             label = paste("Internal Guideline: ", 350, sep = ""), angle = 90) +
    background_grid(major = "xy", minor = "y") +
    scale_x_discrete(limits = dna_data_table$SAMPLE_ID))
}


#########################
# Remaining DNA metrics #
#########################

if (dna_sample_count > 0) {

  # Mean family size (DNA)
  print(ggplot(dna_data_table, aes(x = SAMPLE_ID, y = DNA_MEAN_FAMILY_SIZE)) +
    geom_col(aes(fill = RUN), alpha = 0.8) +
    geom_text(aes(label = DNA_MEAN_FAMILY_SIZE), angle = 90, vjust = 0.5,
              hjust = 1, nudge_y = 2.0) +
    coord_cartesian(ylim = c(1, max(dna_data_table$DNA_MEAN_FAMILY_SIZE)
                           + 10)) +
    guides(fill = guide_legend("Run")) +
    xlab("Sample ID") +
    ylab("Mean family size") +
    ggtitle("[LocalApp DNA metric] Mean family size") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                     hjust = 1, size = 7)) +
    background_grid(major = "xy", minor = "y") +
    scale_x_discrete(limits = dna_data_table$SAMPLE_ID))

  # Number of reads (millions)/mean family size
  print(ggplot(dna_data_table, aes(x = SAMPLE_ID,
                  y = (DNA_TOTAL_PF_READS / 1000000) / DNA_MEAN_FAMILY_SIZE)) +
    geom_col(aes(fill = RUN), alpha = 0.8) +
    guides(fill = guide_legend("Run")) +
    xlab("Sample ID") +
    ylab("Number of reads (millions)/mean family size") +
    scale_y_continuous(breaks = 20 * c(0:10)) +
    ggtitle(paste("[LocalApp DNA metric] Number of reads divided by mean",
                  " read family size", sep = "")) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                     hjust = 1, size = 7)) +
    background_grid(major = "xy", minor = "y") +
    scale_x_discrete(limits = dna_data_table$SAMPLE_ID))

  # Percentage of aligned reads (DNA)
  print(ggplot(dna_data_table, aes(x = SAMPLE_ID, y = DNA_PCT_ALIGNED_READS)) +
    geom_col(aes(fill = RUN), alpha = 0.8) +
    coord_cartesian(xlim = c(1, dna_sample_count + 1), ylim = c(0, 104)) +
    guides(fill = guide_legend("Run")) +
    xlab("Sample ID") +
    ylab("Percentage of aligned reads") +
    scale_y_continuous(breaks = 10 * c(0:10)) +
    ggtitle("[LocalApp DNA metric] Percentage of aligned reads") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                     hjust = 1, size = 7)) +
    background_grid(major = "xy", minor = "y") +
    scale_x_discrete(limits = dna_data_table$SAMPLE_ID))

  # Median target coverage (DNA)
  print(ggplot(dna_data_table, aes(x = SAMPLE_ID,
                                   y = DNA_MEDIAN_TARGET_COVERAGE)) +
    geom_col(aes(fill = RUN), alpha = 0.8) +
    coord_cartesian(xlim = c(1, dna_sample_count + 1)) +
    guides(fill = guide_legend("Run")) +
    xlab("Sample ID") +
    ylab("Median target coverage") +
    ggtitle("[LocalApp DNA metric] Median target coverage") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                     hjust = 1, size = 7)) +
    background_grid(major = "xy", minor = "y") +
    scale_x_discrete(limits = dna_data_table$SAMPLE_ID))

  # Percentage of exons at 100X coverage (DNA)
  print(ggplot(dna_data_table, aes(x = SAMPLE_ID, y = DNA_PCT_EXON_100X)) +
    geom_col(aes(fill = RUN), alpha = 0.8) +
    coord_cartesian(xlim = c(1, dna_sample_count + 1), ylim = c(0, 104)) +
    guides(fill = guide_legend("Run")) +
    xlab("Sample ID") +
    ylab("Percentage of exons at 100X coverage") +
    ggtitle("[LocalApp DNA metric] Percentage of exons at 100X coverage") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                     hjust = 1, size = 7)) +
    background_grid(major = "xy", minor = "y") +
    scale_x_discrete(limits = dna_data_table$SAMPLE_ID))

  # Percentage of target at 100X coverage (DNA)
  print(ggplot(dna_data_table, aes(x = SAMPLE_ID, y = DNA_PCT_TARGET_100X)) +
    geom_col(aes(fill = RUN), alpha = 0.8) +
    coord_cartesian(xlim = c(1, dna_sample_count + 1), ylim = c(0, 104)) +
    guides(fill = guide_legend("Run")) +
    xlab("Sample ID") +
    ylab("Percentage of target at 100X coverage") +
    ggtitle("[LocalApp DNA metric] Percentage of target at 100X coverage") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                     hjust = 1, size = 7)) +
    background_grid(major = "xy", minor = "y") +
    scale_x_discrete(limits = dna_data_table$SAMPLE_ID))

  # Percentage of target at 250X coverage (DNA)
  print(ggplot(dna_data_table, aes(x = SAMPLE_ID, y = DNA_PCT_TARGET_250X)) +
    geom_col(aes(fill = RUN), alpha = 0.8) +
    coord_cartesian(xlim = c(1, dna_sample_count + 1), ylim = c(0, 104)) +
    guides(fill = guide_legend("Run")) +
    xlab("Sample ID") +
    ylab("Percentage of target at 250X coverage") +
    ggtitle("[LocalApp DNA metric] Percentage of target at 250X coverage") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                     hjust = 1, size = 7)) +
    background_grid(major = "xy", minor = "y") +
    scale_x_discrete(limits = dna_data_table$SAMPLE_ID))

  # Percentage of target at 0.4X mean coverage (DNA)
  print(ggplot(dna_data_table, aes(x = SAMPLE_ID,
                                   y = DNA_PCT_TARGET_0.4X_MEAN)) +
    geom_col(aes(fill = RUN), alpha = 0.8) +
    coord_cartesian(xlim = c(1, dna_sample_count + 1), ylim = c(0, 104)) +
    guides(fill = guide_legend("Run")) +
    xlab("Sample ID") +
    ylab("Percentage of target at 0.4 * mean coverage") +
    ggtitle(paste("[LocalApp DNA metric] Percentage of target",
                  " at 0.4 * mean coverage", sep = "")) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                     hjust = 1, size = 7)) +
    background_grid(major = "xy", minor = "y") +
    scale_x_discrete(limits = dna_data_table$SAMPLE_ID))

  # Percentage of UQ reads passing filters (DNA)
  print(ggplot(dna_data_table, aes(x = SAMPLE_ID, y = DNA_PCT_PF_UQ_READS)) +
    geom_col(aes(fill = RUN), alpha = 0.8) +
    coord_cartesian(xlim = c(1, dna_sample_count + 1), ylim = c(0, 104)) +
    guides(fill = guide_legend("Run")) +
    xlab("Sample ID") +
    ylab("Percentage of UQ reads passing filters") +
    ggtitle("[LocalApp DNA metric] Percentage of UQ reads passing filters") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                     hjust = 1, size = 7)) +
    background_grid(major = "xy", minor = "y") +
    scale_x_discrete(limits = dna_data_table$SAMPLE_ID))

  # Read enrichment (DNA)
  print(ggplot(dna_data_table, aes(x = SAMPLE_ID,
                                   y = DNA_PCT_READ_ENRICHMENT)) +
    geom_col(aes(fill = RUN), alpha = 0.8) +
    coord_cartesian(xlim = c(1, dna_sample_count + 1), ylim = c(0, 104)) +
    guides(fill = guide_legend("Run")) +
    xlab("Sample ID") +
    ylab("Read enrichment (percent)") +
    ggtitle("[LocalApp DNA metric] Read enrichment") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                     hjust = 1, size = 7)) +
    background_grid(major = "xy", minor = "y") +
    scale_x_discrete(limits = dna_data_table$SAMPLE_ID))

  # Contamination estimate (DNA)
  print(ggplot(dna_data_table, aes(x = SAMPLE_ID,
                                   y = DNA_PCT_CONTAMINATION_EST)) +
    geom_col(aes(fill = RUN), alpha = 0.8) +
    coord_cartesian(xlim = c(1, dna_sample_count + 1), ylim = c(0, 104)) +
    geom_text(aes(label = DNA_PCT_CONTAMINATION_EST), angle = 90, vjust = 0.5,
              hjust = 1, nudge_y = 10) +
    guides(fill = guide_legend("Run")) +
    xlab("Sample ID") +
    ylab("Contamination estimate (percent)") +
    ggtitle("[LocalApp DNA metric] Contamination estimate") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                     hjust = 1, size = 7)) +
    background_grid(major = "xy", minor = "y") +
    scale_x_discrete(limits = dna_data_table$SAMPLE_ID))

  # Usable UMI reads (DNA)
  print(ggplot(dna_data_table, aes(x = SAMPLE_ID,
                                   y = DNA_PCT_USABLE_UMI_READS)) +
    geom_col(aes(fill = RUN), alpha = 0.8) +
    coord_cartesian(xlim = c(1, dna_sample_count + 1), ylim = c(0, 104)) +
    guides(fill = guide_legend("Run")) +
    xlab("Sample ID") +
    ylab("Usable UMI reads (percent)") +
    ggtitle("[LocalApp DNA metric] Usable UMI reads") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                     hjust = 1, size = 7)) +
    background_grid(major = "xy", minor = "y") +
    scale_x_discrete(limits = dna_data_table$SAMPLE_ID))

  # Percentage of chimeric reads (DNA)
  print(ggplot(dna_data_table, aes(x = SAMPLE_ID, y = DNA_PCT_CHIMERIC_READS)) +
    geom_col(aes(fill = RUN), alpha = 0.8) +
    coord_cartesian(xlim = c(1, dna_sample_count + 1)) +
    guides(fill = guide_legend("Run")) +
    xlab("Sample ID") +
    ylab("Percentage of chimeric reads") +
    ggtitle("[LocalApp DNA metric] Percentage of chimeric reads") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                     hjust = 1, size = 7)) +
    background_grid(major = "xy", minor = "y") +
    scale_x_discrete(limits = dna_data_table$SAMPLE_ID))
}

#########################
# Remaining RNA metrics #
#########################

if (rna_sample_count > 0) {

  # Percentage of chimeric reads (RNA)
  print(ggplot(rna_data_table, aes(x = SAMPLE_ID, y = RNA_PCT_CHIMERIC_READS)) +
    geom_col(aes(fill = RUN), alpha = 0.8) +
    coord_cartesian(xlim = c(1, rna_sample_count + 1)) +
    guides(fill = guide_legend("Run")) +
    xlab("Sample ID") +
    ylab("Percentage of chimeric reads") +
    ggtitle("[LocalApp RNA metric] Percentage of chimeric reads") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                     hjust = 1, size = 7)) +
    background_grid(major = "xy", minor = "y") +
    scale_x_discrete(limits = rna_data_table$SAMPLE_ID))

  # Total filter-passing reads in millions (RNA)
  print(ggplot(rna_data_table, aes(x = SAMPLE_ID,
                                   y = RNA_TOTAL_PF_READS / 1000000)) +
    geom_col(aes(fill = RUN), alpha = 0.8) +
    coord_cartesian(xlim = c(1, rna_sample_count + 1)) +
    guides(fill = guide_legend("Run")) +
    xlab("Sample ID") +
    ylab("Total reads that pass quality filters (millions)") +
    ggtitle(paste("[LocalApp RNA metric]",
                  " Total filter-passing reads in millions",
                  "\n(the input is down-sampled to 30 million reads",
                  " [not read pairs])", sep = "")) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                     hjust = 1, size = 7)) +
    background_grid(major = "xy", minor = "y") +
    scale_x_discrete(limits = rna_data_table$SAMPLE_ID))

  # Percentage of on-target reads (RNA)
  print(ggplot(rna_data_table, aes(x = SAMPLE_ID,
                                   y = RNA_PCT_ON_TARGET_READS)) +
    geom_col(aes(fill = RUN), alpha = 0.8) +
    coord_cartesian(xlim = c(1, rna_sample_count + 1), ylim = c(0, 100)) +
    guides(fill = guide_legend("Run")) +
    xlab("Sample ID") +
    ylab("Percentage of on-target reads") +
    ggtitle("[LocalApp RNA metric] Percentage of on-target reads") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                     hjust = 1, size = 7)) +
    background_grid(major = "xy", minor = "y") +
    scale_x_discrete(limits = rna_data_table$SAMPLE_ID))

  # Scaled median gene coverage (RNA)
  print(ggplot(rna_data_table, aes(x = SAMPLE_ID,
                                   y = RNA_SCALED_MEDIAN_GENE_COVERAGE)) +
    geom_col(aes(fill = RUN), alpha = 0.8) +
    coord_cartesian(xlim = c(1, rna_sample_count + 1)) +
    guides(fill = guide_legend("Run")) +
    xlab("Sample ID") +
    ylab("Scaled median gene coverage") +
    ggtitle("[LocalApp RNA metric] Scaled median gene coverage") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                     hjust = 1, size = 7)) +
    background_grid(major = "xy", minor = "y") +
    scale_x_discrete(limits = rna_data_table$SAMPLE_ID))
}
dev.off()
