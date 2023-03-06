#' Plot exon structure 
#'
#' @param ... gene model, primers or other intervals in GRanges
#'
#' @export
#' @import  ggplot2
plotSegments <- function(...) {
  grLists <- list(...)
  # remove metadata
  grLists <- lapply(grLists, unlist)
  for (i in seq_along(grLists)) {
      S4Vectors::mcols(grLists[[i]]) <- NULL
      grLists[[i]] <- split(grLists[[i]], names(grLists[[i]]))
  }
  grLists <- do.call(c, grLists)
  segs <- gr2segments(unname(grLists))
  # connecting lines
  lines <- do.call(rbind, 
                   lapply(split(segs, segs$track_id),
                          function(x) {
                            data.frame(start = min(x$start), end = max(x$end))
                          }))
  lines$track_id <- rownames(lines)
  q <- ggplot2::ggplot() + 
    ggplot2::geom_segment(data = segs,
                 mapping = ggplot2::aes(
                   x = start,
                   xend = end,
                   y = track_id,
                   yend = track_id),
                 size = 8) + 
    ggplot2::geom_segment(
      data = lines,
      mapping = ggplot2::aes(
        x = start,
        xend = end,
        y = track_id,
        yend = track_id
      ), size = 1, colour = "black"
    ) +  
    ggplot2::geom_vline(xintercept = unique(c(segs$start, segs$end)),
               colour = "grey80") + 
    ggplot2::theme_classic() + 
    ggplot2::theme(axis.text.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank()) +
    ggplot2::ylab("")  +
    ggplot2::xlab("")
  print(q)
}

#' Split primers GRanges into a list by their direction
#'
#' @param grPrimers a GRanges with the forward and reverse primers 
#'
#' @return a list with separate GRanges for every direction
#' @export
#'
primersToList <- function(grPrimers) {
  result <- lapply(names(grPrimers), function(primersName) {
    if (!is.null(grPrimers[[primersName]]$direction)) {
      primers <- split(grPrimers[[primersName]], grPrimers[[primersName]]$direction)
      names(primers) <- paste(primersName, names(primers), sep = "-")
      primers
    } else {
      grPrimers[[primersName]]
    }
  })
  do.call(c, result)
}



gr2segments <- function(gr, normalise = TRUE) {
  gr <- do.call(c, as.list(gr))
  gr$track_id <- names(gr)
  gr <- unname(gr)
  gr <- as.data.frame(gr)
  if (normalise) {
    norm <- normaliseData(gr = gr)
    gr[c("start", "end")] <- norm$gr
  }
  gr 
}

#' Substitutes coordinates by their ranks (among all points).
#'
#' @param x an integer vector (e.g. starts)
#' @param y an integer vector (e.g. ends)
#'
#' @return a list with the ranked coordinates, list(x_ranks, y_ranks)
#'
#' @noRd
unifyDiff <- function(x,y) {
  points <- cbind(x,y)
  o <- order(points)
  res <- rle(points[o]) 
  res$values <- seq_along(unique(res$values))
  points[o] <- inverse.rle(res)
  list(points[seq_along(x)],
       points[seq_along(y) + length(x)])
}

#' Apply unifyDiff on a list of intervals
#'
#' @param ... data.frames with `start` and `end` columns
#'
#' @return
#'
#' @noRd
normaliseData <- function(...){
  dat <- list(...)
  toProcess <- vapply(dat, Negate(is.null), logical(1))
  columns <- c("start", "end")
  positions <- do.call(rbind,
                       lapply(names(dat)[toProcess], function(x) {
                         cbind(id = x, dat[[x]][, columns])
                       }))
  result <- as.data.frame(
    do.call(cbind,
            unifyDiff(positions$start, positions$end))
  )
  names(result) <- columns
  dat[toProcess] <- split(result, positions$id)
  dat
}


  
