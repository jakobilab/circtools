#!/usr/bin/env Rscript

# Copyright (C) 2025 Tobias Jakobi
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

args <- commandArgs(trailingOnly = TRUE)
base_path <- args[1]

# we need these packages
pkgs <- c(
    "aod",
    "amap",
    "ballgown",
    "devtools",
    "biomaRt",
    "data.table",
    "edgeR",
    "GenomicFeatures",
    "GenomicRanges",
    "ggbio",
    "ggfortify",
    "ggplot2",
    "gplots",
    "ggrepel",
    "gridExtra",
    "openxlsx",
    "plyr",
    "Hmisc",
    "reshape2",
    "devtools",
    "kableExtra",
    "formattable",
    "dplyr",
    "RColorBrewer",
    "devtools",
    'GenomicRanges',
    'ggplot2',
    'BSgenome',
    'IRanges',
    'S4Vectors',
    'Biostrings',
    "openxlsx",
    "ggrepel",
    "aod",
    "plyr"
)

countdown <- function(from)
{
  cat(from)
  while(from!=0)
  {
    Sys.sleep(1)
    from <- from - 1
    cat("\r",from)
  }
}
message("")

# check if devtools is already installed
pkgs <- pkgs[!pkgs %in% installed.packages()[,1]]

minorVersion <- as.numeric(strsplit(version[['minor']], '')[[1]][[1]])
majorVersion <- as.numeric(strsplit(version[['major']], '')[[1]][[1]])

message("")
message("This script will automatically install R packages required by circtools.")
message("")

message(paste("Detected R version ", majorVersion, ".", version[['minor']], "\n", sep=""))

message("Detected library paths:")
for (path in .libPaths()){
    message(paste0("-> ",path))
}
message("")

for (package in pkgs){
    message(paste("Need to install package", package))
}

if (
    majorVersion >= 4
    || ( majorVersion == 3 && minorVersion >= 6 )
){
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("ggplot2")
        install.packages("BiocManager")


        if (length(pkgs) > 0)
            BiocManager::install(pkgs)

} else {
    source("https://bioconductor.org/biocLite.R")
    biocLite()


    if (length(pkgs) > 0)
        biocLite(pkgs)
}

# load devtools library
library(devtools)

message("")
message("Now installing R CircTest and primex R packages.")
message("")

install.packages(paste0(base_path,"/contrib/primex"),
                         repos = NULL,
                          type = "source")


install.packages(paste0(base_path,"/contrib/circtest"),
                         repos = NULL,
                          type = "source")


