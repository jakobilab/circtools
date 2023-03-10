.toGenome <- function(x, granges) {
  granges <- as.data.frame(granges)
  up   <- granges[1,, drop = FALSE]
  down <- granges[2,, drop = FALSE]
  strand <- unique(granges$strand)
  inUpstream <- x <= up$width
  if (strand == '+') {
    x[ inUpstream] <- x[ inUpstream] + up$start - 1
    x[!inUpstream] <- x[!inUpstream] - up$width + down$start - 1
  } else {
    x[ inUpstream] <- -x[ inUpstream] + up$end + 1
    x[!inUpstream] <- -x[!inUpstream] + down$end + up$width + 1
  }
  x
}

extractCoords <- function(coords) {
  coords <- strsplit(coords, ",") %>% lapply(as.integer)
  stopifnot(all(lapply(coords, length) == 2))
  do.call(rbind, coords) %>% 
    as.data.frame() %>% 
    stats::setNames(c("start", "width"))
}

#' Return primer coordinates.
#' 
#' In Primer3 LEFT and RIGHT are 0-based.
#' RIGHT start is 5'-end of the primer, we coerce it to IRanges style.
#'
#' @param primers 
#' @noRd
primersToPairs <- function(primers) {
  lefts <- extractCoords(primers$PRIMER_LEFT)
  lefts$start <- lefts$start + 1
  rights <- extractCoords(primers$PRIMER_RIGHT)
  rights$start <- rights$start - rights$width + 2
  lapply(seq_along(lefts), function(i) {
    cbind(rbind(lefts[i,], rights[i,]), 
          direction = c("forward", "reverse"))
  })
}


splitPrimer <- function(x, upExonWidth) {
  res <- lapply(1:2, function(i) {
    y <- x[i,]
    if (y$start <= upExonWidth && y$start + y$width - 1 > upExonWidth) {
      y <- rbind(y, y, make.row.names = FALSE)
      y$width[1] <- upExonWidth - y$start[1] + 1
      y$start[2] <- upExonWidth + 1
      y$width[2] <- y$width[2] - y$width[1]
    }
    y
  })
  do.call(rbind, res)
}

insertJunctions <- function(pairs, exons) {
  force(pairs)
  lapply(pairs, splitPrimer, IRanges::width(exons)[1])
}

pairToGenome <- function(pair, exons) {
  ends   <- .toGenome(pair$start + pair$width - 1, granges = exons)
  starts <- .toGenome(pair$start, granges = exons)
  res <- do.call(rbind, Map(range, starts, ends))
  seqnames <- S4Vectors::Rle(GenomicRanges::seqnames(exons[1]))
  strands  <- S4Vectors::Rle(GenomicRanges::strand(exons[1]))
  res <- GenomicRanges::GRanges(seqnames,
                         IRanges::IRanges(start = res[, 1], end = res[, 2]),
                         strand = strands)
  S4Vectors::mcols(res) <- pair[!names(pair) %in% c("start", "width")]
  res
}

#' Transforms output of design to GRanges
#'
#' @param primers the "primers" item of the `\link{design}` output
#' @param exons the GRanges of the exons
#'
#' @return GRangesList
#' @export
#'
toGRanges <- function(primers, exons) {
  pairs <- primers %>%
    primersToPairs() %>%
    insertJunctions(exons) 
  GenomicRanges::GRangesList(lapply(pairs, pairToGenome, exons))
}
