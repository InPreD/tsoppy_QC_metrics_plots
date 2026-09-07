###################################################################
## Author: Tonje (2020-12-04); edited by Daniel (2021-04-07)     ##
##                                                               ##
## This script uses a transposed Illumina metric file as input   ##
## and prints out a QC metrics figure for a selected RNA sample. ##
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
# construction of data from the input file:
qc_variables <- c("COMPLETED_ALL_STEPS",
                  "RNA_MEDIAN_CV_GENE_500X",
                  "RNA_TOTAL_ON_TARGET_READS",
                  "RNA_MEDIAN_INSERT_SIZE",
                  "RNA_PCT_CHIMERIC_READS",
                  "RNA_PCT_ON_TARGET_READS",
                  "RNA_SCALED_MEDIAN_GENE_COVERAGE",
                  "RNA_TOTAL_PF_READS")

qc_values <- data.frame("qc_variables" = qc_variables,
             "lowest_possible_value"  = c(NA,   0,  0,   0,   0,   0,    0,  0),
             "threshold"              = c(NA,  93,  9,  80,   0,   0,    0,  0),
             "highest_possible_value" = c(NA, 130, 30, 200, 0.3, 100, 7000, 30),
             "obs_values" = NA)


# construction of data as it will be plotted:
y_legend <- c("PIPELINE COMPLETED",
             "CV FOR GENES AT MEDIAN COVERAGE >500",
             "TOTAL ON-TARGET READS (MILLIONS)",
             "MEDIAN INSERT SIZE",
             "CHIMERIC READS (%)",
             "ON-TARGET READS (%)",
             "SCALED MEDIAN GENE COVERAGE",
             "TOTAL READS PASSING FILTERS (MILLIONS)")

plotting_data <- data.frame("yText" = y_legend,
                            "yText_T" = c("",
                                          "<= 93   [fusions, splice variants]",
                                          ">= 9   [fusions, splice variants]",
                                          ">= 80   [fusions, splice variants]",
                                          "No reference value",
                                          "No reference value",
                                          "No reference value",
                                          "No reference value"),
                            "QCtype" = c("TrueFalse",
                                         "PassFail",
                                         "FailPass",
                                         "FailPass",
                                         "NoThold",
                                         "NoThold",
                                         "NoThold",
                                         "NoThold"),
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
      sample_type <- "Sample type: unknown"
    }
  }
}


# select the observed data from the input file:
qc_values$obs_values <- as.vector(as.matrix(run_table[
                                        run_table$SAMPLE_ID == input_id,
                                        as.character(qc_values$qc_variables)]))

n <- nrow(plotting_data)
plotting_data$obs_values <- qc_values$obs_values
plotting_data$obs_values[3] <- round(plotting_data$obs_values[3] / 1000000, 1)
plotting_data$obs_values[8] <- round(plotting_data$obs_values[8] / 1000000, 1)

# plot
for (device in c("pdf", "png")) {

  if (device == "pdf") {
    pdf(output_pdf, width = 9)
  } else {
    png(output_png, width = 9, height = 7, units = "in", res = 150)
  }

  par(mar = c(5, 22, 5, 2))
  plot(x = 1, y = 1, xlim = c(-0.2, 3.5), ylim = c(0, n), xlab = "", ylab = "",
       yaxt = "n", xaxt = "n", col = "white",
       main = paste("QC for RNA sample \"", output_id, "\"\n",
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

  # empty conclusion boxes to the far right:
  for (i in 1:n) {
    rect(xleft = 2.8, xright = 3.2, ytop = n - i + 0.5, ybottom = n - i)
  }

  # start plotting the "TrueFalse" QC:
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

  # plotting the continuous variables:
  # we need the low/T/high settings now
  plotting_data$lowest_possible_value  <- qc_values$lowest_possible_value
  plotting_data$threshold              <- qc_values$threshold
  plotting_data$highest_possible_value <- qc_values$highest_possible_value

  for (i in 2:n) {
  # select colors:
    if (plotting_data[i, "QCtype"] == "FailPass") {
      color_left <- color_fail; color_right <- color_pass
    }
    if (plotting_data[i, "QCtype"] == "NoThold") {
      color_left <- "white"; color_right <- "white"
    }
    if (plotting_data[i, "QCtype"] == "PassFail") {
      color_left <- color_pass; color_right <- color_fail
    }

    # scale all plotted values to range [0, 2]
    # draw the boxes:
    a <- plotting_data[i, -c(1:3)] /
                                 plotting_data[i, "highest_possible_value"] * 2
    a$obs_value_color <- "white"
    a$obs_value_color[a$obs_values > a$highest_possible_value] <- "red"
    rect(xleft = 0, xright = a$threshold, ytop = n - i + 0.5,
         ybottom = n - i, col = color_left)
    rect(xleft = a$threshold, xright = 2, ytop = n - i + 0.5, ybottom = n - i,
         col = color_right)
    rect(xleft = min(a$obs_values, 2) - 0.01,
         xright = min(a$obs_values, 2) + 0.01,
         ytop = n - i + 0.6, ybottom = n - i - 0.1,
         col = a$obs_value_color)

    # write the min/max values:
    text(x = -0.1, y = n - i + 0.25,
         labels = plotting_data[i, "lowest_possible_value"])
    text(x = 2.2, y = n - i + 0.25,
         labels = plotting_data[i, "highest_possible_value"])
    text(x = min(a$obs_values, 2), y = n - i + 0.75,
         labels = plotting_data[i, "obs_values"])
  }

dev.off()
}
