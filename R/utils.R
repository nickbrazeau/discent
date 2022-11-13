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
  setTxtProgressBar(pb_list[[name]], i)
  if (i == max_i) {
    close(pb_list[[name]])
  }
}

#------------------------------------------------
#' @title Simple function for expanding a pairwise matrix
#' @param discdf dataframe; that is specific to the DISCent framework. Must have colnames of: "smpl1", "smpl2", "deme1", "deme2", "gendist", "geodist" in correct order
#' @details This is an internal function that largely lacks generalizability. Do not recommend external use.
#' @export

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
  return(discdfexpand)
}


