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

#------------------------------------------------
#' @title Simple function for simulating a distance matrix
#' @param nDemes integer; number of demes to consider
#' @param rateDist numeric; rate in exponential distribution (for drawing distances)
#' @details Draws distances from a squared-normal distribution
#' @returns Distance matrix of class \code{dist}
#' @export

drawDistMat <- function(nDemes = 5, rateDist = 1e-3) {
  ddists <- rexp(n = ncol(combn(nDemes, 2)), rateDist)
  # catch any 0s
  ddists[ddists == 0] <- 1
  distmat <- matrix(NA, nrow = nDemes, ncol = nDemes)
  distmat[ lower.tri(distmat) ] <- ddists
  distmat[ upper.tri(distmat) ] <- t(distmat)[ lower.tri(distmat) ]
  # selves to 0
  diag(distmat) <- 0
  return(as.dist(distmat))
}
