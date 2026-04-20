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

# ── Step 1: Bootstrap BiocManager, then force BiocGenerics upgrade FIRST ──────
if (majorVersion >= 4 || (majorVersion == 3 && minorVersion >= 6)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }

  biocgenerics_ok <- tryCatch({
    current <- packageVersion("BiocGenerics")
    current >= package_version("0.53.2")
  }, error = function(e) FALSE)

  if (!biocgenerics_ok) {
    message("BiocGenerics < 0.53.2 or missing — installing 'generics' then BiocGenerics from Bioc 3.21...")
    lib_path <- .libPaths()[1]

    # Install 'generics' from CRAN first (new dependency of BiocGenerics 0.54.x)
    if (!requireNamespace("generics", quietly = TRUE)) {
      install.packages("generics", repos = "https://cloud.r-project.org", lib = lib_path)
    }

    install.packages(
      "BiocGenerics",
      repos = "https://bioconductor.org/packages/3.21/bioc",
      lib   = lib_path,
      type  = "source"
    )

    if ("package:BiocGenerics" %in% search()) detach("package:BiocGenerics", unload = TRUE, force = TRUE)
    library(BiocGenerics, lib.loc = lib_path)
    message(paste("BiocGenerics version now:", packageVersion("BiocGenerics")))
  } else {
    message(paste("BiocGenerics", packageVersion("BiocGenerics"), "already satisfies >= 0.53.2, skipping."))
  }
}

# ── Step 2: Determine what needs installing ────────────────────────────────────
bioc_pkgs <- c(
  "ballgown", "biomaRt", "edgeR",
  "GenomicFeatures", "GenomicRanges", "ggbio",
  "BSgenome", "IRanges", "S4Vectors", "Biostrings"
)

cran_pkgs <- c(
  "aod", "amap", "devtools", "data.table", "ggfortify", "ggplot2",
  "gplots", "ggrepel", "gridExtra", "openxlsx", "plyr",
  "reshape2", "kableExtra", "formattable", "dplyr", "RColorBrewer",
  "readr"
)

# Filter out already installed
bioc_pkgs <- bioc_pkgs[!bioc_pkgs %in% installed.packages()[, 1]]
cran_pkgs <- cran_pkgs[!cran_pkgs %in% installed.packages()[, 1]]

for (package in c(bioc_pkgs, cran_pkgs)) {
  message(paste("Need to install package", package))
}

# ── Step 3: Install packages ───────────────────────────────────────────────────
if (majorVersion >= 4 || (majorVersion == 3 && minorVersion >= 6)) {

  if (length(bioc_pkgs) > 0) {
    message("Installing Bioconductor packages...")
    BiocManager::install(bioc_pkgs, ask = FALSE, update = FALSE)
  }

  if (length(cran_pkgs) > 0) {
    message("Installing CRAN packages...")
    install.packages(cran_pkgs, repos = "https://cloud.r-project.org")
  }

} else {
  source("https://bioconductor.org/biocLite.R")
  biocLite()
  if (length(c(bioc_pkgs, cran_pkgs)) > 0) biocLite(c(bioc_pkgs, cran_pkgs))
}

# ── Step 4: Verify all packages installed correctly ───────────────────────────
all_pkgs <- c(bioc_pkgs, cran_pkgs)
missing <- all_pkgs[!all_pkgs %in% installed.packages()[, 1]]
if (length(missing) > 0) {
  stop(paste("ERROR: The following packages failed to install:",
             paste(missing, collapse = ", ")))
}

# ── Step 5: Archive packages ───────────────────────────────────────────────────
message("\nInstalling archived R packages...")

install.packages("https://cran.r-project.org/src/contrib/Archive/Hmisc/Hmisc_4.6-0.tar.gz",
                 repos = NULL, type = "source")
install.packages("https://cran.r-project.org/src/contrib/Archive/GGally/GGally_2.1.2.tar.gz",
                 repos = NULL, type = "source")
install.packages("https://cran.r-project.org/src/contrib/Archive/ggstats/ggstats_0.3.0.tar.gz",
                 repos = NULL, type = "source")

# ── Step 6: Local source installs ─────────────────────────────────────────────
message("\nInstalling local R packages (primex, circtest)...")

install.packages(paste0(base_path, "/contrib/primex"), repos = NULL, type = "source")
install.packages(paste0(base_path, "/contrib/circtest"), repos = NULL, type = "source")

message("\n✅ All R dependencies for circtools are installed.")