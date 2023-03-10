#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`


#' Set primer size
#'
#' @param settings a settings list created by `p3Settings`
#' @param min minimal length (PRIMER_MIN_SIZE)
#' @param optimal optimal length (PRIMER_OPT_SIZE)
#' @param max maximal length (PRIMER_MAX_SIZE)
#'
#' @return a settings list with updated values
#' @export
#'
#' @examples
#' p3Settings() %>%
#'   primerSize(min = 18, optimal = 20, max = 22) %>%
#'   str()
primerSize <- function(settings, min = NULL, optimal = NULL, max = NULL) {
    stopifnot(min <= optimal, optimal <= max)
    settings["PRIMER_MIN_SIZE"] <- list(min)
    settings["PRIMER_OPT_SIZE"] <- list(optimal)
    settings["PRIMER_MAX_SIZE"] <- list(max)
    settings
}

#' Set ptimer melting temperatures
#'
#' @param settings a settings list created by `p3Settings`
#' @param min minimal temperature (PRIMER_MAX_TM)
#' @param optimal optimal temperature (PRIMER_OPT_TM)
#' @param max maximal temperature (PRIMER_MIN_TM)
#'
#' @return a settings list with updated values
#' @export
#'
#' @examples
#' p3Settings() %>%
#'   primerTm(min = 58, optimal = 60, max = 62) %>%
#'   str()
primerTm <- function(settings, min = NULL, optimal = NULL, max = NULL) {
  stopifnot(min <= optimal, optimal <= max)
  settings["PRIMER_MAX_TM"] <- list(max)
  settings["PRIMER_OPT_TM"] <- list(optimal)
  settings["PRIMER_MIN_TM"] <- list(min)
  settings
}

#' Set product size 
#'
#' @inheritParams primerSize
#' @param range a vector with two numbers, c(minimal, maximal) product size
#' (inclusive).
#'
#' @return a settings list with updated values
#' @export
#'
#' @examples
#' str(p3Settings() %>% productSize(c(100,300)))
#' 
productSize <- function(settings, range) {
  stopifnot(length(range) == 2) 
  stopifnot(range[1] <= range[2])
  range <- paste(as.integer(range), collapse = "-")
  settings["PRIMER_PRODUCT_SIZE_RANGE"] <- range
  settings
}