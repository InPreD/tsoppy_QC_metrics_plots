###################################################################
## Author: Tonje (2020-11-05); edited by Daniel (2021-04-14)     ##
##                                                               ##
## This script uses a transposed Illumina metric file as input   ##
## and prints out a QC metrics figure for a selected DNA sample. ##
###################################################################

version_string <- "0.3:22-02-23"

args <- commandArgs(trailingOnly = TRUE)

# test whether the correct number of command-line arguments was supplied
# - if not, exit with a help/error message
if (length(args) != 6) {
  stop(paste("Exactly 6 arguments must be supplied (",
             "1: path to the input TSV file; ",
             "2: path to the output PDF file; ",
             "3: path to the output PNG file; ",
             "4: ID of a sample which should be selected from the input file; ",
             "5: sample ID string displayed in the output",
             "6: InPreD nomenclature",
             ")\n", sep = ""), call. = FALSE)
}

# assign the command-line arguments to variables
input_tsv <- args[1]            # path to the input TSV file
output_pdf <- args[2]           # path to the output PDF file
output_png <- args[3]           # path to the output PNG file
input_id <- args[4]             # ID of a selected input file sample
output_id <- args[5]            # sample ID displayed in the output
inpred_nomenclature <- args[6]  # "True"|"False": InPreD ID nomenclature usage

##################################################
##################################################
# load libraries
library("reshape2")
library("graphics")

##################################################
##################################################
# select colors
color_fail <- "#E69F00"
color_pass <- "#009E73"

#################################################
#################################################
## construction of data from the input file:
qc_variables <- c("COMPLETED_ALL_STEPS",
                  "DNA_CONTAMINATION_SCORE",
                  "DNA_CONTAMINATION_P_VALUE",
                  "DNA_MEDIAN_EXON_COVERAGE",
                  "DNA_PCT_EXON_50X",
                  "DNA_MEDIAN_INSERT_SIZE",
                  "DNA_PCT_ALIGNED_READS",
                  "DNA_USABLE_MSI_SITES",
                  "DNA_MEDIAN_BIN_COUNT_CNV_TARGET",
                  "DNA_COVERAGE_MAD")

qc_values <- data.frame("qc_variables" = qc_variables,
        "lowest_possible_value" = c(NA, 0, 0, 0, 0, 0, 0, 0, 0, 0),
        "threshold" = c(NA, 3106, 0.049, 150, 90, 70, 70, 40, 1, 0.21),
        "highest_possible_value" = c(NA, NA, 1, 1000, 100, 200, 100,
                                     130, 20, 0.3),
        "obs_values" = NA)

# construction of data as it will be plotted:
# - "CONTAMINATION" is a YES/NO variable combining two Illumina QC metrics
#   (DNA_CONTAMINATION_SCORE and DNA_CONTAMINATION_P_VALUE)
y_legend <- c("PIPELINE COMPLETED",
              "CONTAMINATION",
              "MEDIAN EXON COVERAGE",
              "EXON BASES WITH COVERAGE >=50 (%)",
              "MEDIAN INSERT SIZE",
              "ALIGNED READS (%)",
              "USABLE MSI SITES",
              "MEDIAN RAW BIN COUNT PER CNV TARGET",
              "CNV COVERAGE MAD")

plotting_data <- data.frame("yText" = y_legend,
                            "yText_T" =  c("",
                                           "",
                                           ">= 150   [small variants, TMB]",
                                           ">= 90   [small variants, TMB]",
                                           ">= 70   [small variants, TMB]",
                                           ">= 70",
                                           ">= 40   [MSI]",
                                           ">= 1   [CNV]",
                                           "<= 0.21   [CNV]"),
                            "QCtype" = c("TrueFalse",
                                        "TrueFalse",
                                        "FailPass",
                                        "FailPass",
                                        "FailPass",
                                        "FailPass",
                                        "FailPass",
                                        "FailPass",
                                        "PassFail"),
                            "obs_values" = NA)

#################################################
#################################################
##### read the input data:
run_table <- read.table(input_tsv, header = T,
                        sep = "\t", stringsAsFactors = FALSE)

############################################
############################################
# sample type information:
sample_type <- "Sample type unspecified"

if (inpred_nomenclature == "True") {
  if (nchar(output_id) == 19) {
    sample_type_code <- unlist(strsplit(output_id, split = ""))[13]
    if (sample_type_code == "M") {
      sample_type <- "(metastasis)"
    } else if (sample_type_code == "T") {
      sample_type <- "(tumor)"
    } else if (sample_type_code == "N") {
      sample_type <- "(normal)"
    } else if (sample_type_code == "C") {
      sample_type <- "(cell-line)"
    } else if (sample_type_code == "L") {
      sample_type <- "(liquid sample)"
    } else if (sample_type_code == "P") {
      sample_type <- "(primary tumor, naive)"
    } else if (sample_type_code == "p") {
      sample_type <- "(primary tumor, post-treatment)"
    } else if (sample_type_code == "R") {
      sample_type <- "(regional metastasis, naive)"
    } else if (sample_type_code == "r") {
      sample_type <- "(regional metastasis, post-treatment)"
    } else if (sample_type_code == "D") {
      sample_type <- "(distal metastasis, naive)"
    } else if (sample_type_code == "d") {
      sample_type <- "(distal metastasis, post-treatment)"
    } else {
      sample_type <- "Sample type unknown"
    }
  }
}

# select the observed data from the input file:
qc_values$obs_values <- as.vector(as.matrix(run_table[
                                         run_table$SAMPLE_ID == input_id,
                                         as.character(qc_values$qc_variables)]))

# construct the plotting data:
n <- nrow(plotting_data)
plotting_data$obs_values[1] <- qc_values$obs_values[1]
# - contamination
if (qc_values$obs_values[2] > qc_values$threshold[2]
    && qc_values$obs_values[3] > qc_values$threshold[3]) {
  plotting_data$obs_values[2] <- 1
} else {
  plotting_data$obs_values[2] <- 0
}
plotting_data$obs_values[3:n] <- qc_values$obs_values[4:(n + 1)]

## plot
for (device in c("pdf", "png")) {

  if (device == "pdf") {
    pdf(output_pdf, width = 9)
  } else {
    png(output_png, width = 9, height = 7, units = "in", res = 150)
  }

  par(mar = c(5, 22, 5, 2))
  plot(x = 1, y = 1, xlim = c(-0.2, 3.5), ylim = c(0, n), xlab = "", ylab = "",
       yaxt = "n", xaxt = "n", col = "white",
       main = paste("QC for DNA sample \"", output_id, "\"\n",
                    sample_type, sep = ""))
  abline(v = 2.5)
  axis(side = 1, at = c(1, 3), labels = c("Observed values", "Conclusions"),
       las = 1, tick = FALSE)
  axis(side = 2, at = ((n - 1):0) + 0.4 - 0.1, labels = plotting_data$yText,
       las = 2, tick = FALSE, font = 2)
  axis(side = 2, at = ((n - 1):0) + 0.10 - 0.1, labels = plotting_data$yText_T,
       las = 2, cex.axis = 0.8, tick = FALSE)
  #text(x = -0.5, y = -2,
  #     label = "Overall conclusion and comments (with date and initials):",
  #     xpd = NA, font = 2)
  ## empty conclusion boxes to the far right:
  for (i in 1:n) {
    rect(xleft = 2.8, xright = 3.2,
         ytop = n - i + 0.5, ybottom = n - i)
  }

  ## start plotting the "TrueFalse" QC:
  if (plotting_data$obs_values[1] == 1) {
    item_color <- color_pass
    item_label <- "YES"
  }
  if (plotting_data$obs_values[1] == 0) {
    item_color <- color_fail
    item_label <- "NO"
  }
  rect(xleft = 0.8, xright = 1.2, ytop = n - 0.5, ybottom = n - 1,
       col = item_color)
  text(x = 1, y = n - 0.75, labels = item_label)

  if (plotting_data$obs_values[2] == 1) {
    item_color <- color_fail
    item_label <- "YES"
  }
  if (plotting_data$obs_values[2] == 0) {
    item_color <- color_pass
    item_label <- "NO"
  }
  rect(xleft = 0.8, xright = 1.2, ytop = n - 1.5, ybottom = n - 2,
       col = item_color)
  text(x = 1, y = n - 1.75, labels = item_label)

  ## plotting the continuous variables:
  #we need the low/T/high settings now
  plotting_data$lowest_possible_value <- c(NA,
                                           NA,
                                      qc_values$lowest_possible_value[-c(1:3)])
  plotting_data$threshold <- c(NA, NA, qc_values$threshold[-c(1:3)])
  plotting_data$highest_possible_value <- c(NA,
                                            NA,
                                      qc_values$highest_possible_value[-c(1:3)])

  for (i in 3:n) {
    if (plotting_data[i, "QCtype"] == "FailPass") {
      color_left <- color_fail
      color_right <- color_pass
    }
    if (plotting_data[i, "QCtype"] == "NoThold") {
      color_left <- "white"
      color_right <- "white"
    }
    if (plotting_data[i, "QCtype"] == "PassFail") {
      color_left <- color_pass
      color_right <- color_fail
    }

    # scale all plotted values to range [0, 2]
    # draw the boxes (all values from x=0 to x=2).
    a <- plotting_data[i, -c(1:3)] /
                                plotting_data[i, "highest_possible_value"] * 2
    a$observed_value_color <- "white"
    a$observed_value_color[a$obs_values > a$highest_possible_value] <- "red"
    rect(xleft = 0, xright = a$threshold, ytop = n - i + 0.5, ybottom = n - i,
         col = color_left)
    rect(xleft = a$threshold, xright = 2, ytop = n - i + 0.5, ybottom = n - i,
         col = color_right)
    rect(xleft = min(a$obs_values, 2) - 0.01,
         xright = min(a$obs_values, 2) + 0.01,
         ytop = n - i + 0.6, ybottom = n - i - 0.1,
         col = a$observed_value_color)

    #write the min/max values:
    text(x = - 0.1, y = n - i + 0.25,
         labels = plotting_data[i, "lowest_possible_value"])
    text(x = 2.2, y = n - i + 0.25,
         labels = plotting_data[i, "highest_possible_value"])
    text(x = min(a$obs_values, 2), y = n - i + 0.75,
         labels = plotting_data[i, "obs_values"])
  }

dev.off()
}
