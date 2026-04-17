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

pkgs <- c(
  "aod", "amap", "ballgown", "devtools", "biomaRt", "data.table", "edgeR",
  "GenomicFeatures", "GenomicRanges", "ggbio", "ggfortify", "ggplot2",
  "gplots", "ggrepel", "gridExtra", "openxlsx", "plyr",
  "reshape2", "kableExtra", "formattable", "dplyr", "RColorBrewer",
  "BSgenome", "IRanges", "S4Vectors", "Biostrings", "readr"
)

# Remove already installed packages
pkgs <- pkgs[!pkgs %in% installed.packages()[, 1]]

minorVersion <- as.numeric(strsplit(version[['minor']], '')[[1]][[1]])
majorVersion <- as.numeric(strsplit(version[['major']], '')[[1]][[1]])

message("")
message("This script will automatically install R packages required by circtools.\n")
message(paste("Detected R version ", majorVersion, ".", version[['minor']], "\n", sep=""))
message("Detected library paths:")
for (path in .libPaths()) message(paste0("-> ", path))
message("")

for (package in pkgs) {
  message(paste("Need to install package", package))
}

if (majorVersion >= 4 || (majorVersion == 3 && minorVersion >= 6)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos="https://cloud.r-project.org")
  }

  #bypass both
  bioc_version <- BiocManager::version()
  bioc_minor <- as.numeric(strsplit(as.character(bioc_version), "\\.")[[1]][2])
  if (bioc_minor <= 19) {
    message("Bioconductor <= 3.19 detected — installing BiocGenerics and S4Vectors from Bioc 3.20 tarballs...")
    lib_path <- .libPaths()[1]

    install.packages(
      "https://bioconductor.org/packages/3.20/bioc/src/contrib/BiocGenerics_0.54.0.tar.gz",
      repos = NULL, type = "source", lib = lib_path
    )
    install.packages(
      "https://bioconductor.org/packages/3.20/bioc/src/contrib/S4Vectors_0.44.0.tar.gz",
      repos = NULL, type = "source", lib = lib_path
    )

    if ("package:BiocGenerics" %in% search()) detach("package:BiocGenerics", unload = TRUE, force = TRUE)
    if ("package:S4Vectors" %in% search()) detach("package:S4Vectors", unload = TRUE, force = TRUE)

    library(BiocGenerics, lib.loc = lib_path)
    library(S4Vectors, lib.loc = lib_path)

    message(paste("BiocGenerics version now:", packageVersion("BiocGenerics")))
    message(paste("S4Vectors version now:", packageVersion("S4Vectors")))

    # Remove both from pkgs so BiocManager doesn't overwrite them
    pkgs <- pkgs[!pkgs %in% c("BiocGenerics", "S4Vectors")]
  }

  if (length(pkgs) > 0) {
    BiocManager::install(pkgs, ask = FALSE, update = FALSE)
  }

} else {
  source("https://bioconductor.org/biocLite.R")
  biocLite()
  if (length(pkgs) > 0) biocLite(pkgs)
}

# --- Archive packages ---
message("\nInstalling archived R packages...")

install.packages("https://cran.r-project.org/src/contrib/Archive/Hmisc/Hmisc_4.6-0.tar.gz",
                 repos = NULL, type = "source")
install.packages("https://cran.r-project.org/src/contrib/Archive/GGally/GGally_2.1.2.tar.gz",
                 repos = NULL, type = "source")
install.packages("https://cran.r-project.org/src/contrib/Archive/ggstats/ggstats_0.3.0.tar.gz",
                 repos = NULL, type = "source")


# --- Local source installs ---
message("\nInstalling local R packages (primex, circtest)...")

install.packages(paste0(base_path, "/contrib/primex"), repos = NULL, type = "source")
install.packages(paste0(base_path, "/contrib/circtest"), repos = NULL, type = "source")

message("\n✅ All R dependencies for circtools are installed.")