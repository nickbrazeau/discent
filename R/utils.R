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
  return( log(p/(1-p)) )
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
  colnames(yexpand) <- c("smpl2", "smpl1", "deme2", "deme1", "gendist", "geodist")
  yexpand <- rbind.data.frame(y, yexpand) # now have all pairwise possibilities
  yexpand <- yexpand[!duplicated(yexpand), ]
  return(yexpand)
}


#------------------------------------------------
#' @title Truncated Normal Distrubtion
#' @noRd
rnorm_interval <- function(mean, sd, a=0, b=1) {

  # draw raw value relative to a
  ret <- rnorm(1, mean, sd) - a
  # reflect off boundries at 0 and (b-a)
  if (ret < 0 || ret > (b-a)) {
    # use multiple reflections to bring into range [-(b-a), 2(b-a)]
    while (ret < -(b-a)) {
      ret <- ret + 2*(b-a)
    }
    while (ret > 2*(b-a)) {
      ret <- ret - 2*(b-a)
    }
    # use one more reflection to bring into range [0,(b-a)]
    if (ret < 0) {
      ret <- -ret
    }
    if (ret > (b-a)) {
      ret <- 2*(b-a) - ret
    }
  }
  # no longer relative to a
  ret <- ret + a
  return(ret)
}
