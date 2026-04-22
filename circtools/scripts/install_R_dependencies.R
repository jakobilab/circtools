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

# NOTE: Hmisc, GGally, ggstats, ggbio, biovizBase are excluded from the main
# BiocManager pass. They are pinned archive versions installed from source.
# biovizBase and ggbio depend on Hmisc and are installed after it.
# Install order: Hmisc -> ggstats -> GGally (GGally imports ggstats at lazy-load time)

pkgs <- c(
  "aod", "amap", "ballgown", "devtools", "biomaRt", "data.table", "edgeR",
  "GenomicFeatures", "GenomicRanges", "ggfortify", "ggplot2",
  "gplots", "ggrepel", "gridExtra", "openxlsx", "plyr",
  "reshape2", "kableExtra", "formattable", "dplyr", "RColorBrewer",
  "BSgenome", "IRanges", "S4Vectors", "Biostrings", "readr"
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


  message("\nPre-installing all dependencies for archived R packages...")
  archive_deps <- c(
    # Hmisc 4.6-0 Imports
    "latticeExtra", "Formula", "base64enc", "htmltools", "htmlTable",
    "viridis", "cluster", "foreign", "gtable", "nnet", "rpart",
    # Hmisc 4.6-0 Depends (non-base)
    "lattice", "survival",

    # ggstats 0.5.0 Imports
    "broom.helpers", "cli", "magrittr", "patchwork", "purrr",
    "stringr", "forcats", "lifecycle", "rlang", "scales", "tidyr",

    # GGally 2.2.1 Imports (beyond what ggstats already covers)
    "progress",

    # Also in pkgs but required before Step 2 (archive installs) —
    # safe to list here too since pkgs filters already-installed packages
    "ggplot2", "dplyr", "gridExtra", "plyr", "RColorBrewer", "data.table"
  )
  archive_deps <- archive_deps[!archive_deps %in% installed.packages()[, 1]]
  if (length(archive_deps) > 0) {
    BiocManager::install(archive_deps, ask = FALSE, update = FALSE, lib = lib_path)
  }

  # Verify all archive deps installed before proceeding
  still_missing_deps <- archive_deps[!archive_deps %in% installed.packages()[, 1]]
  if (length(still_missing_deps) > 0) {
    stop(paste("Could not install archive dependencies:",
               paste(still_missing_deps, collapse = ", ")))
  }

  # --- Step 2: Install pinned archive packages ---
  # URL tarballs require repos=NULL. Dependencies are pre-installed in Step 1.
  message("\nInstalling archived R packages...")

  tryCatch({
    install.packages("https://cran.r-project.org/src/contrib/Archive/Hmisc/Hmisc_4.6-0.tar.gz",
                     repos = NULL, type = "source", lib = lib_path)
  }, error = function(e) {
    stop(paste("Hmisc archive install failed:", e$message))
  })
  if (!"Hmisc" %in% installed.packages()[, 1]) {
    stop("Hmisc archive install failed - non-zero exit status")
  }
  message(paste("Hmisc installed, version:", packageVersion("Hmisc")))

  # ggstats BEFORE GGally — GGally 2.2.1 imports ggstats at lazy-load time
  tryCatch({
    install.packages("https://cran.r-project.org/src/contrib/Archive/ggstats/ggstats_0.5.0.tar.gz",
                     repos = NULL, type = "source", lib = lib_path)
  }, error = function(e) {
    stop(paste("ggstats archive install failed:", e$message))
  })
  if (!"ggstats" %in% installed.packages()[, 1]) {
    stop("ggstats archive install failed - non-zero exit status")
  }
  message(paste("ggstats installed, version:", packageVersion("ggstats")))

  tryCatch({
    install.packages("https://cran.r-project.org/src/contrib/Archive/GGally/GGally_2.2.1.tar.gz",
                     repos = NULL, type = "source", lib = lib_path)
  }, error = function(e) {
    stop(paste("GGally archive install failed:", e$message))
  })
  if (!"GGally" %in% installed.packages()[, 1]) {
    stop("GGally archive install failed - non-zero exit status")
  }
  message(paste("GGally installed, version:", packageVersion("GGally")))

  # --- Step 3: Install biovizBase and ggbio now that Hmisc is pinned ---
  message("\nInstalling biovizBase and ggbio (depend on pinned Hmisc)...")
  if (!"biovizBase" %in% installed.packages()[, 1]) {
    BiocManager::install("biovizBase", ask = FALSE, update = FALSE, lib = lib_path)
  }
  if (!"ggbio" %in% installed.packages()[, 1]) {
    BiocManager::install("ggbio", ask = FALSE, update = FALSE, lib = lib_path)
  }

  # --- Step 4: Main BiocManager install ---
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
  "Hmisc", "GGally", "ggstats"
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
