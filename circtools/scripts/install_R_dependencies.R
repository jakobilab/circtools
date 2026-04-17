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

minorVersion <- as.numeric(strsplit(version[['minor']], '')[[1]][[1]])
majorVersion <- as.numeric(strsplit(version[['major']], '')[[1]][[1]])

message("")
message("This script will automatically install R packages required by circtools.\n")
message(paste("Detected R version ", majorVersion, ".", version[['minor']], "\n", sep=""))
message("Detected library paths:")
for (path in .libPaths()) message(paste0("-> ", path))
message("")

# ── Step 1: Bootstrap BiocManager and upgrade Bioconductor FIRST ──────────────
if (majorVersion >= 4 || (majorVersion == 3 && minorVersion >= 6)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }

  bioc_version <- BiocManager::version()
  bioc_minor <- as.numeric(strsplit(as.character(bioc_version), "\\.")[[1]][2])

  if (bioc_minor <= 19) {
    message("Bioconductor <= 3.19 detected — upgrading to 3.20 before installing packages...")
    BiocManager::install(version = "3.20", ask = FALSE, update = TRUE)
    message(paste("Bioconductor upgraded. BiocGenerics version now:", packageVersion("BiocGenerics")))
  }
}

# ── Step 2: Now determine what needs installing ───────────────────────────────
pkgs <- c(
  "aod", "amap", "ballgown", "devtools", "biomaRt", "data.table", "edgeR",
  "GenomicFeatures", "GenomicRanges", "ggbio", "ggfortify", "ggplot2",
  "gplots", "ggrepel", "gridExtra", "openxlsx", "plyr",
  "reshape2", "kableExtra", "formattable", "dplyr", "RColorBrewer",
  "BSgenome", "IRanges", "S4Vectors", "Biostrings", "readr"
)

# Remove already installed packages
pkgs <- pkgs[!pkgs %in% installed.packages()[, 1]]

for (package in pkgs) {
  message(paste("Need to install package", package))
}

# ── Step 3: Install packages against the now-correct Bioc version ─────────────
if (majorVersion >= 4 || (majorVersion == 3 && minorVersion >= 6)) {
  if (length(pkgs) > 0) {
    BiocManager::install(pkgs, ask = FALSE, update = FALSE)
  }
} else {
  source("https://bioconductor.org/biocLite.R")
  biocLite()
  if (length(pkgs) > 0) biocLite(pkgs)
}

# ── Step 4: Archive packages ──────────────────────────────────────────────────
message("\nInstalling archived R packages...")

install.packages("https://cran.r-project.org/src/contrib/Archive/Hmisc/Hmisc_4.6-0.tar.gz",
                 repos = NULL, type = "source")
install.packages("https://cran.r-project.org/src/contrib/Archive/GGally/GGally_2.1.2.tar.gz",
                 repos = NULL, type = "source")
install.packages("https://cran.r-project.org/src/contrib/Archive/ggstats/ggstats_0.3.0.tar.gz",
                 repos = NULL, type = "source")

# ── Step 5: Local source installs ────────────────────────────────────────────
message("\nInstalling local R packages (primex, circtest)...")

install.packages(paste0(base_path, "/contrib/primex"), repos = NULL, type = "source")
install.packages(paste0(base_path, "/contrib/circtest"), repos = NULL, type = "source")

message("\n✅ All R dependencies for circtools are installed.")