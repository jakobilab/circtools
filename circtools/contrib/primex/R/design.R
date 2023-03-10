
#' A helper to define template parameters
#'
#' @param settings a list with predefined parameters
#' @param granges A `GRanges` object with `seq` column. It can include one or
#'   two elements. In the case of two, the sequences are concatenated and
#'   one of the primers will be forced two overlap the splice junction.
#' @param seq (if gr is NULL) a character vector of length 1 or 2.
#'   In the case of two sequences, they are treated as consequent exons and
#'   either of two primers will be forced two overlap the splice junction.
#' @param seqId a character name for the sequence
#' @return a named list with all template parameters (with defaults for the ones
#'  not specified in  the `seqOpts`  or in the function arguments. 
#'  The function arguments, i.e. `seq1` etc. overwrite the ones provided in
#'  the `seqOpts`.
#' @export
#' @examples 
#' seqOpts <- seqSettings(seqId = "seq1", seq = "AATCTGAATCGCGCTTAAAGCTA")
#'
seqSettings <- function(settings = NULL,
                        granges = NULL,
                        seqId = NULL,
                        seq = NULL) {
  if (!is.null(granges)) {
    stopifnot(inherits(granges, "GRanges"))
    stopifnot(length(granges) %in% 1:2)
    stopifnot(!is.null(mcols(granges)$seq))
    if (is.null(seqId))
      seqId <- paste(granges$exon_id, collapse = "|")
    seq <- S4Vectors::mcols(granges)$seq
  }
  # defined on the basis of the Primer3 html doc
  defaults <- list(
    SEQUENCE_EXCLUDED_REGION = NULL,
    SEQUENCE_INCLUDED_REGION = NULL,
    SEQUENCE_PRIMER_REVCOMP = NULL,
    SEQUENCE_FORCE_LEFT_END = NULL,
    SEQUENCE_INTERNAL_EXCLUDED_REGION = NULL,
    SEQUENCE_QUALITY = NULL,
    SEQUENCE_FORCE_LEFT_START = NULL,
    SEQUENCE_INTERNAL_OLIGO = NULL,
    SEQUENCE_START_CODON_POSITION = NULL,
    SEQUENCE_FORCE_RIGHT_END = NULL,
    SEQUENCE_OVERLAP_JUNCTION_LIST = NULL,
    SEQUENCE_TARGET = NULL,
    SEQUENCE_FORCE_RIGHT_START = NULL,
    SEQUENCE_PRIMER = NULL,
    SEQUENCE_TEMPLATE = NULL,
    SEQUENCE_ID = NULL,
    SEQUENCE_PRIMER_PAIR_OK_REGION_LIST = NULL
  )
  defaults[names(settings)] <- settings
  defaults$SEQUENCE_ID <- ifelse(is.null(seqId), "sequence1", seqId)
  defaults$SEQUENCE_TEMPLATE <- paste(seq, collapse = "")
  if (length(seq) == 2)
    defaults$SEQUENCE_OVERLAP_JUNCTION_LIST <- nchar(seq[1])
  defaults
}

#' Run primer design
#'
#' @param seqOpts parameters for the given sequence, 
#'   i.e. template sequence, its name, junction position
#' @param primerOpts (optional) thermodynamical and other primer-specific
#'  parameters. Can be created using \link{p3Settings}.
#' @param returnStats (default: TRUE) if the "EXPLAIN" fields with 
#'   the run statistics should be returned. This option overwrites 
#'   the one provided in the `primerOpts$PRIMER_EXPLAIN_FLAG`.
#' @param ... path parameters for the `\link{runPrimer3}` call.
#'
#' @return a data.frame with the designed primers
#' 
#' @details By default, if `primer3Config` is not secified/NULL,
#'   it is inferred on the basis of the `primer3Path` as a subdirectory 
#'   `primer3_config`
#' @export
#' @examples 
#' # make up a sequence
#' seq1 <- paste(c("A", "T", "G", "C")[sample(4, 300, replace = TRUE)], 
#'               collapse = "")
#' seqOpts <- seqSettings(seqId = "seq1", seq = seq1)
#' res <- design(seqOpts)
#' res$primers
#'
design <- function(seqOpts, primerOpts  = NULL, returnStats = TRUE, ...) {
  # we put all in a one input file
  allOpts <- seqOpts
  allOpts[names(primerOpts)] <- primerOpts
  allOpts$PRIMER_EXPLAIN_FLAG <- ifelse(returnStats, "1", "0")
  result <- runPrimer3(settings = allOpts, ...)
  if (!inherits(result, "try-error")) 
    result <- extractPrimers(result) 
  result
}

#' Invokes Primer3 for a given list of options
#'
#' @param settings a named list with all the run options.
#' @param path a list(primer3=, config=) with the path to the Primer3 executable
#'   and the configuration file. By default, it uses the files provided with the
#'   package installation.
#'
#' @return a named list
#' @export
#'
runPrimer3 <- function(settings,path = list(primer3 = NULL, config  = NULL)) {
  path$primer3 <- pickPrimer3Exec(path$primer3)
  path$config  <- pickPrimer3Config(path$config, path$primer3)
  settings <- p3Settings(settings)
  settings$PRIMER_THERMODYNAMIC_PARAMETERS_PATH <- path$config
  inputFile <- tempfile()
  writeLines(listToP3(settings), inputFile)
  result <- try(system(paste(path$primer3, inputFile), intern = TRUE))
  unlink(inputFile)
  result
}

.errorIfNotExists <- function(x, type = "") {
  if (!file.exists(x)) 
    stop("File ", type, ": ", x)
}

# transform Primer3 output to a data.frame describing pairs
# (`primers` list item) . Other  output values are returned 
# in the `option` list item.
extractPrimers <- function(result) {
  result <- p3ToList(result)
  primersInexes <- grep("(PRIMER_LEFT|PRIMER_RIGHT|PRIMER_PAIR)_\\d+",
                        names(result))
  primers <- NULL
  if (length(primersInexes) > 0) {
    primers <- result[primersInexes]
    result <- result[-primersInexes]
    pairNumbers <- vapply(strsplit(names(primers), "_"),
                          `[[`, 3, FUN.VALUE = character(1))
    primers <- do.call(rbind, split(unlist(primers),pairNumbers))
    colnames(primers) <- sub("_0_","_", x = colnames(primers))
    colnames(primers) <- sub("_0$","", x = colnames(primers))
    primers <- data.frame(primers, stringsAsFactors = FALSE)
  } 
  list(primers = primers, options = result)
}

listToP3 <- function(x) {
  x <- x[!vapply(x, is.null, logical(1))]
  x <- paste(names(x), x, sep = "=")
  # must have "=" sign in the very end of the file
  c(x, "=")
}

p3ToList <- function(lines) {
  lines <- grep(".+=", lines, value = TRUE) 
  lines <- strsplit(lines, "=")
  result <- lapply(lines, `[`, 2)
  names(result) <- lapply(lines, `[`, 1)
  result
}

#' Get/define Primer3 settings
#'
#' @param settings user provided list of non-default settings
#' @param defaultsFile (optional) a path to the file with default settings
#'   (e.g. "primer3/primer3_v1_1_4_default_settings.txt")
#'   
#'
#' @return a named list of the settings
#' @export
#'
p3Settings <- function(settings = NULL, defaultsFile = NULL) {
  if (is.null(defaultsFile)) {
    pkgDir <- find.package('primex', .libPaths())
    defaultsFile <- normalizePath(
      file.path(pkgDir,"primer3", "primer3_v1_1_4_default_settings.txt"))
    if (!file.exists(defaultsFile))
      stop("File with the default settings is not found! \n",
           defaultsFile)
  }
  settingsTemplate <- readLines(defaultsFile)
  defaults <- p3ToList(settingsTemplate)
  if ("" %in% names(settings)) {
    stop("Some arguments have empty names.")
  }
  defaults[names(settings)] <- settings
  defaults
}

pickPrimer3Exec <- function(file) {
  if (is.null(file)) {
    executables <- list(
      Linux =  list(x86_64 = "primer3_core_Linux_64",
                    i686   = "primer3_core_Linux_32"),
      Darwin = list(x86_64 = "primer3_core_Darwin_64",
                    i686   = NULL) #"primer3_core_Linux_i686"),
    )
    info <- as.list(Sys.info())
    machine <- ifelse(grepl("64", info$machine), "x86_64", "i686")
    file <- executables[[info$sysname]][[machine]]
    pkgDir <- find.package('primex', .libPaths())
    file <- normalizePath(file.path(pkgDir, "primer3", file))
  }
  .errorIfNotExists(file, "for the Primer3 executable  is not found.")
  file
}

pickPrimer3Config <- function(primer3Config, primer3Path) {
  if (is.null(primer3Config)) {
    primer3Config <- normalizePath(file.path(dirname(primer3Path),
                                            "primer3_config"))
    primer3Config <- paste0(primer3Config, .Platform$file.sep)
  }
  .errorIfNotExists(primer3Config, " for the Primer3 config is not found!")
  primer3Config
}

#' Returns diagnostics information from the run
#' 
#' For this to work, `PRIMER_EXPLAIN_FLAG` in the must be set to "1".
#'
#' @param result a value returned by \link{runPrimer3}
#'
#' @return a list with statistics of the run for the left, right primers
#'  and their pairs.
#' @export
#'
#' @examples
#' # Too short sequence
#' seq1 <- paste(c("A", "T", "G", "C")[sample(4, 100, replace = TRUE)], 
#'               collapse = "")
#' seqOpts <- seqSettings(seqId = "seq1", seq = seq1)
#' res <- design(seqOpts)
#' diagnose(res)
#' 
diagnose <- function(result) {
  expl <- grep("EXPLAIN", names(result$options), value = TRUE)
  expl <- grep("FLAG", expl, value = TRUE,invert = TRUE)
  result$options[expl]
}


