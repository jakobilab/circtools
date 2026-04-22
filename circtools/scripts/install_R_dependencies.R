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
  "GenomicFeatures", "GenomicRanges", "ggfortify", "ggplot2",
  "gplots", "ggrepel", "gridExtra", "openxlsx", "plyr",
  "reshape2", "kableExtra", "formattable", "dplyr", "RColorBrewer",
  "BSgenome", "IRanges", "S4Vectors", "Biostrings", "readr",
  "Hmisc", "biovizBase", "ggbio"
)

# Remove already installed packages
pkgs <- pkgs[!pkgs %in% installed.packages()[, 1]]

minorVersion <- as.numeric(strsplit(version[['minor']], '')[[1]][[1]])
majorVersion <- as.numeric(strsplit(version[['major']], '')[[1]][[1]])

# Determine explicit library path (first writable path)
lib_path <- .libPaths()[1]

message("")
message("This script will automatically install R packages required by circtools.\n")
message(paste("Detected R version ", majorVersion, ".", version[['minor']], "\n", sep=""))
message(paste("Installing all packages to:", lib_path))
message("Detected library paths:")
for (path in .libPaths()) message(paste0("-> ", path))
message("")

for (package in pkgs) {
  message(paste("Need to install package", package))
}

if (majorVersion >= 4 || (majorVersion == 3 && minorVersion >= 6)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos="https://cloud.r-project.org", lib = lib_path)
  }

  # --- Step 1: Pre-install deps for scales 1.3.0 before pinning it ---
  # scales must be installed from archive BEFORE BiocManager runs, otherwise
  # BiocManager pulls in current scales (>= 1.4.0) which removed trans_breaks,
  # trans_format, and math_format used in circtools_quickcheck_wrapper.R
  message("\nPre-installing dependencies for pinned scales 1.3.0...")
  scales_deps <- c(
    "munsell", "R6", "viridisLite", "labeling", "farver", "glue",
    "cli", "lifecycle", "rlang", "RColorBrewer"
  )
  scales_deps <- scales_deps[!scales_deps %in% installed.packages()[, 1]]
  if (length(scales_deps) > 0) {
    BiocManager::install(scales_deps, ask = FALSE, update = FALSE, lib = lib_path)
  }

  still_missing <- scales_deps[!scales_deps %in% installed.packages()[, 1]]
  if (length(still_missing) > 0) {
    stop(paste("Could not install scales dependencies:", paste(still_missing, collapse = ", ")))
  }

  # --- Step 2: Pin scales 1.3.0 ---
  message("\nInstalling pinned scales 1.3.0...")
  tryCatch({
    install.packages("https://cran.r-project.org/src/contrib/Archive/scales/scales_1.3.0.tar.gz",
                     repos = NULL, type = "source", lib = lib_path)
  }, error = function(e) {
    stop(paste("scales archive install failed:", e$message))
  })
  if (!"scales" %in% installed.packages()[, 1]) {
    stop("scales archive install failed - non-zero exit status")
  }
  message(paste("scales installed, version:", packageVersion("scales")))

  # --- Step 3: Main BiocManager install ---
  # scales is already pinned so BiocManager won't upgrade it
  message("\nInstalling core packages via BiocManager...")
  if (length(pkgs) > 0) {
    BiocManager::install(pkgs, ask = FALSE, update = FALSE, lib = lib_path)
  }

} else {
  source("https://bioconductor.org/biocLite.R")
  biocLite()
  if (length(pkgs) > 0) biocLite(pkgs)
}

# --- Verify all core package installs succeeded ---
message("\nVerifying core package installations...")
core_pkgs <- c(
  "aod", "amap", "ballgown", "devtools", "biomaRt", "data.table", "edgeR",
  "GenomicFeatures", "GenomicRanges", "ggbio", "ggfortify", "ggplot2",
  "gplots", "ggrepel", "gridExtra", "openxlsx", "plyr",
  "reshape2", "kableExtra", "formattable", "dplyr", "RColorBrewer",
  "BSgenome", "IRanges", "S4Vectors", "Biostrings", "readr",
  "Hmisc", "biovizBase", "ggbio", "scales"
)
failed_pkgs <- core_pkgs[!core_pkgs %in% installed.packages()[, 1]]
if (length(failed_pkgs) > 0) {
  message("Warning: The following packages were not installed, retrying:")
  for (p in failed_pkgs) message(paste0("  -> ", p))
  BiocManager::install(failed_pkgs, ask = FALSE, update = FALSE, lib = lib_path)
  still_missing <- failed_pkgs[!failed_pkgs %in% installed.packages()[, 1]]
  if (length(still_missing) > 0) {
    stop(paste("Could not install required packages:",
               paste(still_missing, collapse = ", ")))
  }
} else {
  message("All core packages installed successfully")
}


# --- Local source installs ---
message("\nInstalling local R packages (primex, CircTest)...")

primex_path <- paste0(base_path, "/contrib/primex")
circtest_path <- paste0(base_path, "/contrib/circtest")

message(paste("primex path:", primex_path))
message(paste("circtest path:", circtest_path))
message(paste("primex exists:", dir.exists(primex_path)))
message(paste("circtest exists:", dir.exists(circtest_path)))

# Install primex with error propagation
tryCatch({
  install.packages(primex_path, repos = NULL, type = "source", lib = lib_path)
}, error = function(e) {
  stop(paste("primex install.packages() threw an error:", e$message))
})

# Verify primex before proceeding to CircTest
if (!"primex" %in% installed.packages(lib.loc = lib_path)[, 1]) {
  stop("primex installation failed - aborting before CircTest install")
}
message("primex installed successfully")

# Make primex available on the search path for CircTest's build
library(primex, lib.loc = lib_path)

# Install CircTest with error propagation
tryCatch({
  install.packages(
    circtest_path,
    repos = NULL,
    type = "source",
    lib = lib_path,
    INSTALL_opts = c("--no-multiarch", "--with-keep.source")
  )
}, error = function(e) {
  stop(paste("CircTest install.packages() threw an error:", e$message))
})

# Verify CircTest - if present but broken, try loading it for a descriptive error
if (!"CircTest" %in% installed.packages(lib.loc = lib_path)[, 1]) {
  tryCatch(
    library(CircTest, lib.loc = lib_path),
    error = function(e) stop(paste("CircTest failed to load:", e$message))
  )
  stop("CircTest installation failed (package not found after install)")
}
message("CircTest installed successfully")

message("\nAll R dependencies for circtools are installed.")