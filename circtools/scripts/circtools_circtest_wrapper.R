#!/usr/bin/env Rscript

# Copyright (C) 2017 Tobias Jakobi
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

#install.packages("/home/tjakobi/repos/dieterichlab/CircTest/", repos = NULL, type="source")


library(CircTest)

# pre load libraries so we don't get messages later:
suppressMessages(library(aod))
suppressMessages(library(ggplot2))
suppressMessages(library(plyr))

# may not be optimal, but we do not want warning for now
options(warn=-1)

# fix 'could not find zip' error from openxlsx
Sys.setenv(R_ZIPCMD = "/usr/bin/zip")
library(openxlsx)

# fixed params:
MAX_LINES <- 50000

# let's now get the command line arguments

# quiet mode
options(echo=FALSE)
args <- commandArgs(trailingOnly = TRUE)

# read arguments
message("Raw args[1]: '", args[1], "'")
arg_dcc_data <- trimws(gsub("['\"\r\n]", "", args[1]))
message("Sanitized arg_dcc_data: '", arg_dcc_data, "'")
message("basename(arg_dcc_data): '", basename(arg_dcc_data), "'")

arg_replictes <- as.integer(args[2]) # integer
arg_condition_list <- strsplit(args[3],",")[[1]] # list of strings
arg_condition_columns <- lapply(strsplit(args[4],","), as.numeric) # list of integers
arg_condition_columns <- unlist(arg_condition_columns)
arg_output_name <- args[5] # string
arg_max_fdr <-as.numeric(args[6]) # float
arg_max_plots <- as.integer(args[7]) # string
arg_filter_sample <- as.integer(args[8]) # integer
arg_filter_count <- as.integer(args[9]) # integer
arg_groups <-   unlist(lapply(strsplit(args[10],","), as.numeric)) # list of strings
arg_output_label <- args[11] # string
arg_percent_filter <-as.numeric(args[12]) # float
arg_only_negative <-as.logical(args[13]) # float
arg_add_header <-as.logical(args[14]) # float
arg_range <-as.numeric(args[15]) # float
arg_colour_mode <- args[16]

# check colour mode
if (arg_colour_mode != "colour" & arg_colour_mode != "bw" ) {
    print(arg_colour_mode)
    message("Please specify the colour model: colour or bw")
    quit()
}


run_CircTest = function(CircRNACount, LinearCount, CircCoordinates, groups, indicators, label, filename, filer.sample,
                        filter.count, max_fdr, max.plots, replicates, percent_filter, only_negative_direction, header, range, colour_mode) {

    message("Filtering circRNA counts")
    CircRNACount_filtered <- Circ.filter(circ = CircRNACount, linear = LinearCount,
    Nreplicates = replicates, filter.sample = filer.sample, filter.count =  filter.count, percentage = percent_filter)

    message("Filtering circRNA coordinates")
    CircCoordinates_filtered <- CircCoordinates[rownames(CircRNACount_filtered),]

    message("Filtering linear RNA counts")
    LinearCount_filtered <- LinearCount[rownames(CircRNACount_filtered),]

    message("running circTest")

    data = Circ.test(CircRNACount_filtered, LinearCount_filtered, CircCoordinates_filtered, group = groups, alpha=max_fdr)

    message("Generating plots")

    max <- nrow(data$summary_table)

    # if we havbe too many results, cut them down to the user threshold
    if (nrow(data$summary_table) > max.plots) {
        max <- max.plots
    }

    # only if we have any results
    if (nrow(data$summary_table) > 0){
        pdf(file = paste(filename,".pdf",sep=""), width= 8.2, height=11.69, title="circtools circtest analysis")

        if (only_negative_direction){
            for (i in rownames(data$summary_table)){
                if(data$summary_table[i,]$direction < 0){
                    invisible(capture.output(Circ.ratioplot(CircRNACount_filtered, LinearCount_filtered, CircCoordinates_filtered,
                    plotrow = i, size = 24, gene_column = 4, groupindicator1 = indicators,
                    x = "", y = "", lab_legend = label, y_axis_range = range, colour_mode = colour_mode)))
                }
            }
        } else {
            for (i in rownames(data$summary_table[1 : max,])){
                invisible(capture.output(Circ.ratioplot(CircRNACount_filtered, LinearCount_filtered, CircCoordinates_filtered,
                plotrow = i, size = 24, gene_column = 4, groupindicator1 = indicators,
                x = "", y = "", lab_legend = label, y_axis_range = range, colour_mode = colour_mode)))
            }
        }
        dev.off()

        message("creating Excel sheets")

        ## Significant result show in a summary table

        write.table(na.omit(data$summary_table[1 : MAX_LINES,]),
                    file = paste(filename, ".csv", sep = "") ,
                    quote = FALSE,
                    sep = "\t",
                    eol = "\n",
                    row.names = FALSE,
                    col.names = header)

        wb <- createWorkbook()

        addWorksheet(wb, sheetName = "Significant circles")
        writeDataTable(wb, sheet = 1, x = data$summary_table[1 : MAX_LINES,])

        idx <- rownames(data$summary_table[1 : MAX_LINES,])

        addWorksheet(wb, sheetName = "Circle Counts")
        writeDataTable(wb, sheet = 2, x = CircRNACount[idx,])

        addWorksheet(wb, sheetName = "Linear Counts")
        writeDataTable(wb, sheet = 3, x = LinearCount[idx,])

        saveWorkbook(wb, paste(filename, ".xlsx", sep = "") , overwrite = TRUE)

    } else {
        message("No canidates to plot, exiting.")
    }
}

###########################################################

## load complete data set
## load complete data set
# --- HDF5 OR directory detection ---
if (grepl("\\.h5$|\\.hdf5$", basename(arg_dcc_data), ignore.case = TRUE)) {
  message("Detected HDF5 input file: ", arg_dcc_data)
  suppressMessages(library(rhdf5))
  
  keys <- h5ls(arg_dcc_data)$name
  message("Available datasets: ", paste(keys, collapse = ", "))

  safe_h5_read <- function(file, key) {
    if (key %in% h5ls(file)$name) {
      df <- as.data.frame(h5read(file, key))
      return(df)
    } else {
      stop(paste("Dataset", key, "not found in", file))
    }
  }

  message("Loading CircRNACount from HDF5")
  CircRNACount <- safe_h5_read(arg_dcc_data, "CircRNACount")

  message("Loading LinearCount from HDF5")
  LinearCount <- safe_h5_read(arg_dcc_data, "LinearCount")

  message("Loading CircCoordinates from HDF5")
  CircCoordinates <- safe_h5_read(arg_dcc_data, "CircCoordinates")

  message("Loaded all datasets from HDF5")

} else {
  # --- Classic directory input ---
  message("Loading CircRNACount File")
  CircRNACount <- read.delim(file.path(arg_dcc_data, "CircRNACount"), header = TRUE)
  
  message("Loading LinearRNACount FIle")
  LinearCount <- read.delim(file.path(arg_dcc_data, "LinearCount"), header = TRUE)
  
  message("Loading CircCoordinates File")
  CircCoordinates <- read.delim(file.path(arg_dcc_data, "CircCoordinates"), header = TRUE)
}

# --- Ensure Strand column exists in LinearCount and CircRNACount ---
# --- Standardize and order columns for all datasets ---

# canonical coordinate order
# --- Standardize and order columns for all datasets ---
coord_cols <- c("Chr", "Start", "End", "Strand")

order_columns <- function(df, name) {
  # add missing coordinate columns *before* reordering
  if (!"Chr" %in% colnames(df)) {
    message(paste("➕ Adding missing Chr column to", name))
    df$Chr <- NA
  }
  if (!"Start" %in% colnames(df)) {
    message(paste("➕ Adding missing Start column to", name))
    df$Start <- NA
  }
  if (!"End" %in% colnames(df)) {
    message(paste("➕ Adding missing End column to", name))
    df$End <- NA
  }
  if (!"Strand" %in% colnames(df)) {
    message(paste("➕ Adding missing Strand column to", name))
    df$Strand <- "+"
  }

  # now reorder cleanly
  present_coords <- intersect(coord_cols, colnames(df))
  remaining_cols <- setdiff(colnames(df), present_coords)
  df <- df[, c(present_coords, remaining_cols), drop = FALSE]

  message(paste("✅ Ordered columns for", name, "→", paste(colnames(df)[1:6], collapse = ", "), "..."))
  return(df)
}

CircCoordinates <- order_columns(CircCoordinates, "CircCoordinates")
CircRNACount    <- order_columns(CircRNACount, "CircRNACount")
LinearCount     <- order_columns(LinearCount, "LinearCount")


group_length <- length(arg_condition_list)

message(paste(group_length, " different groups provided", sep=""))

final_grouping <- unlist(lapply(seq(1, length(arg_condition_columns)), function(x) {return(arg_condition_list[arg_groups[x]])}))

# call the main function to run circTest
run_CircTest(
    CircRNACount[, c(1 : 3, arg_condition_columns)], # we always need the first 3 columns
    LinearCount[, c(1 : 3, arg_condition_columns)], # we always need the first 3 columns
    CircCoordinates,
    arg_groups,
    final_grouping,
    arg_output_label,
    arg_output_name,
    arg_filter_sample,
    arg_filter_count,
    arg_max_fdr,
    arg_max_plots,
    arg_replictes,
    arg_percent_filter,
    arg_only_negative,
    arg_add_header,
    arg_range,
    arg_colour_mode
)
