
#' Get exon pairs
#'
#' @param exonsByTx a list of exons returned by ensembl exonsBy tx
#'
#' @return a list of character vectors
#' @export
#'
extractExonPairs <- function(exonsByTx) {
  stopifnot(allInherit(exonsByTx, "GRanges"))
  lapply(exonsByTx, function(x)  splitToPairs(x))
}

# works for mcols and data.frames
assertColumns <- function(x, ...) {
  columns <- as.character(...)
  missedColumns <- columns[!columns %in% names(S4Vectors::mcols(x))]
  if (length(missedColumns) > 0)
    stop(paste("Missing columns: ",
               paste(missedColumns, collapse = ", ")))
}

# if zero-length return TRUE
allInherit <- function(x, className) {
  all(vapply(x, inherits, className, FUN.VALUE = logical(1)))
}

# we believe that there is at list one exon in a transcript
splitToPairs <- function(x) {
  lapply(seq_along(x[-1]), function(i) {
    c(x[i], x[i + 1])
  })
}

#' Given a list of exons returns candidates for primer design
#' 
#' The transcripts of single exon are excluded!
#'
#' @param exonsByTx exon GRanges by transcripts
#' @param tolerance the minimal tolerable difference in splice junctions 
#'   coordinates (0: return everything, 1: only unique splice junctions).
#'
#' @return a GRangesList
#' @export
#'
selectPairs <- function(exonsByTx, tolerance = 1) {
  expairsbyTx <- extractExonPairs(exonsByTx)
  filterByDistance(expairsbyTx, tolerance)
}



# returns a filtered list
filterByDistance <- function(expairs, minDiff) {
  txIds    <- rep(names(expairs), vapply(expairs, length, integer(1)))
  sjCoords <- do.call(rbind, lapply(unlist(expairs), pair2sj))
  # check if second minimal distance (after self == 0) is above tolerance
  d <- as.matrix(stats::dist(sjCoords, method = "minkowski", p = 1))
  diag(d) <- Inf
  minDist <- apply(d, 1, min)
  tolerated <- split(minDist >= minDiff, txIds)
  Map(`[`, expairs, tolerated)
}

# take end of the first and the start of the second exon
# assume that they are ordered by the exon_rank
pair2sj <- function(p) {
  x <- as.data.frame(p)
  if (x$strand[1] == "-") {
    c(x$end[2], x$start[1])
  } else {
    c(x$end[1], x$start[2])
  }
}


#' Return sequences for the exon pairs
#'
#' @param exonPairs a list of GRanges with the splice junction exon pairs.
#'   For a single pair a single GRanges object may be provided.
#' @param src whether a BSgenome instance, e.g. BSgenome.Hsapiens.NCBI.GRCh38, 
#'   XStringSet object or a FaFile object
#'
#' @return a list of GRanges with additional `seq` data column with
#'   the retrieved sequence.
#' @export
#' @import GenomicRanges
#'
addSeq <- function(exonPairs, src) {
  stopifnot(allInherit(exonPairs, "GRanges"))
  lapply(exonPairs, assertColumns, "exon_id")
  if (!is.list(exonPairs))
    exonPairs <- list(exonPairs)
  # extract all exons by exon_id
  allExons <- do.call(c, exonPairs)
  allExons <- allExons[!duplicated(allExons$exon_id)]
  if (inherits(src, "FaFile")) {
    exseqs <- Rsamtools::getSeq(x = src, allExons)
  } else  if (inherits(src, "BSgenome") || inherits(src, "DNAStringSet")) {
    exseqs <- BSgenome::getSeq(src, allExons)
  } else {
    stop("src must be FaFile, DNAStringSet or BSgenome object")
  }
  exseqs <- as.character(exseqs)
  names(exseqs) <- allExons$exon_id
  for (i in seq_along(exonPairs)) {
      ids <- exonPairs[[i]]$exon_id
      S4Vectors::mcols(exonPairs[[i]])$seq <- exseqs[ids]
  }
  setSJNames(exonPairs)
}

setSJNames <- function(x) {
  pairIds <- lapply(seq_along(x), function(i) {
     paste(x[[i]]$exon_id, collapse = "|")
  })
  names(x) <- pairIds
  x
}