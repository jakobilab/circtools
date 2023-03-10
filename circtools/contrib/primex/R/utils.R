
mockDNAStringSet <- function(seqnames, lens) {
  res <- vapply(lens,
                function(len)
                  paste(sample(
                    c("A", "C", "G", "T"),
                    size = len,
                    replace = TRUE
                  ),
                  collapse = ""), character(1))
  names(res) <- seqnames
  Biostrings::DNAStringSet(res)
}

#' Get Primer3 version 
#'
#' @export
#'
p3Version <- function() {
  
}