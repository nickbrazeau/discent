#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL

#------------------------------------------------
# update progress bar.
# pb_list = list of progress bar objects
# name = name of this progress bar
# i = new value of bar
# max_i = max value of bar (close when reach this value)
#' @noRd
update_progress <- function(pb_list, name, i, max_i) {
  utils::setTxtProgressBar(pb_list[[name]], i)
  if (i == max_i) {
    close(pb_list[[name]])
  }
}

#------------------------------------------------
#' @title Logit Transformation
#' @param p numeric vector of values in (0,1) (probability scale)
#' @description Transforms probability values to the logit scale using the standard
#'   logit formula: \eqn{\text{logit}(p) = \log\left(\frac{p}{1-p}\right)}
#' @return Numeric vector of values on the real line (logit scale)
#' @details This function transforms probabilities (0,1) to the logit scale (real line).
#'   Values of 0 and 1 will produce -Inf and +Inf respectively. The transformation
#'   is used internally in DISC to ensure inbreeding coefficients remain non-negative.
#' @seealso \code{\link{expit}}
#' @export
logit <- function(p){
  # out
   log(p/(1-p))
}

#------------------------------------------------
#' @title Expit (Inverse Logit) Transformation
#' @param p numeric vector of values in logit space
#' @description Converts logit-scale values back to probability scale using the standard
#'   expit formula: \eqn{\text{expit}(p) = \frac{1}{1 + e^{-p}}}
#' @return Numeric vector of values in (0,1) (probability scale)
#' @details This function is the inverse of \code{\link{logit}}. It transforms values
#'   from the logit scale (real line) back to probabilities (0,1).
#' @seealso \code{\link{logit}}
#' @export
expit <- function(p){
  # out
  1/(1+exp(-p))
}

#------------------------------------------------
#' @title Expand Pairwise Distance Matrix
#' @param discdf dataframe with columns: "smpl1", "smpl2", "deme1", "deme2", "gendist", "geodist"
#' @description Internal function to create symmetric pairwise distance matrix by adding
#'   reverse pairs (i,j) -> (j,i). This ensures all pairwise comparisons are represented
#'   in both directions for gradient calculations.
#' @return Expanded dataframe with both (i,j) and (j,i) pairs, duplicates removed
#' @details This function is used internally to ensure symmetric contribution of all
#'   deme pairs to gradient calculations, which is essential for the improved gradient
#'   computation that includes fgrad[j] contributions.
#' @keywords internal
#' @noRd

expand_pairwise <- function(discdf){
  #......................
  # assertions
  #......................
  assert_dataframe(discdf)
  assert_eq(colnames(discdf), c("smpl1", "smpl2", "deme1", "deme2", "gendist", "geodist"))

  #......................
  # expansion
  #......................
  discdfexpand <- discdf
  colnames(discdfexpand) <- c("smpl2", "smpl1", "deme2", "deme1", "gendist", "geodist")
  discdfexpand <- rbind.data.frame(discdf, discdfexpand) # now have all pairwise possibilities
  discdfexpand <- discdfexpand[!duplicated(discdfexpand), ]
  # out
  discdfexpand
}


#............................................................
# DISCresult S3 Class Overloading
#...........................................................
#' @title Check if DISCresult S3 Class
#' @description Tests whether an object is of class DISCresult
#' @param x Object to test
#' @return Logical value: TRUE if object is of class DISCresult, FALSE otherwise
#' @export
is.DISCresult <- function(x) {
  inherits(x, "is.DISCresult")
}

#' @title Print DISCresult S3 Class
#' @description S3 method for printing DISCresult objects with summary information
#' @inheritParams is.DISCresult
#' @param ... Further arguments passed to or from other methods
#' @return Invisibly returns the input object. Called for side effect of printing
#' @export
print.DISCresult <- function(x, ...) {

  # print summary only
  cat(crayon::red("Final DISC Range:"),  paste(round(min(x$Final_Fis),2), round(max(x$Final_Fis),2), sep = " - "), "\n")
  cat(crayon::blue("Final Migration Rate:"), round(x$Final_m, 3), "\n")

  invisible(x)
}

#' @title Summary of DISCresult S3 Class
#' @description S3 method for summarizing DISCresult objects
#' @param object DISCresult object from \code{\link{disc}} function
#' @param ... Further arguments passed to or from other methods
#' @return A list containing tidied final inbreeding coefficients and migration rate
#' @export
summary.DISCresult <- function(object, ...) {
  # send summary only
  tidyout.DISCresult(object)
}


#............................................................
# organizing disc output
#...........................................................
#' @title Tidy Output Generic Method
#' @description Generic method for tidying output from DISC analysis
#' @param x Object to tidy (typically a DISCresult)
#' @return Method-specific tidied output
#' @export
tidyout <- function(x) {
  UseMethod("tidyout")
}

#' @title Tidy DISCresult Output
#' @description S3 method for tidying DISCresult objects into user-friendly format
#' @inheritParams is.DISCresult
#' @return A list with two elements:
#'   \itemize{
#'     \item \code{Final_Fis}: Dataframe with deme names and final inbreeding coefficients
#'     \item \code{Final_M}: Final migration rate estimate
#'   }
#' @export
tidyout.DISCresult <- function(x) {
  #......................
  # clean up
  #......................
  fis <- dplyr::bind_cols(x$deme_key, x$Final_Fis) %>%
    dplyr::select(-key) %>%
    magrittr::set_colnames(c("Deme", "DISC"))

  # out
  list(Final_Fis = fis,
       Final_M = x$Final_m)
}
