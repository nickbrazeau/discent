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
#' @title logit transformation
#' @noRd
# no export because simple
logit <- function(p){
  return(log(p/(1-p)))
}

#------------------------------------------------
#' @title expit transformation
#' @noRd
# no export because simple
expit <- function(p){
  return(1/(1+exp(-p)))
}

#------------------------------------------------
#' @title Simple function for expanding a pairwise matrix
#' @noRd
# no export because lacks generalizability
expand_pairwise <- function(y){
  yexpand <- y
  colnames(yexpand) <- c("smpl2", "smpl1", "locat2", "locat1", "gendist", "geodist")
  yexpand <- rbind.data.frame(y, yexpand) # now have all pairwise possibilities
  yexpand <- yexpand[!duplicated(yexpand), ] # remove duplicate selfs
  return(yexpand)
}

