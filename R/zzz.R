## usethis namespace: start
#' @useDynLib discent, .registration = TRUE
## usethis namespace: end
NULL

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL

#------------------------------------------------
# unload dll when package is unloaded
#' @noRd
.onUnload <- function(libpath) {
  library.dynam.unload("discent", libpath)  # nocov
}
