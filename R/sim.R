#' @title Simulate pairwise Identity by Descent accounting for Isolation by Distance
#' @param demesize numeric vector; Contains the number of individuals residing within
#' the respective deme.
#' @param distmat numeric matrix; A distance matrix corresponding to the distances
#' between demes. The \eqn{m x m} matrix should be of the same length as the number of
#' demes.
#' @param rate numeric; Expected rate of decay for isolation by distance.
#' In effect, a scalar for the distance matrix.
#' @param Ft numeric; Inbreeding coefficient parameter in the population
#' (\emph{i.e.} the total population, \eqn{N_e}).
#' @description
#' \deqn{ \frac{1}{N_i + N_j}(1 - P(D \leq d)) + (1 - \frac{1}{N_i + N_j})(1 - P(D \leq d)) }
#'
#'@export

sim_IBDIBD <- function(demesize, distmat, rate, Ft) {
  #......................
  # assertions
  #......................
  assert_numeric(demesize)
  assert_numeric(distmat)
  assert_numeric(rate)
  assert_numeric(Ft)
  assert_single_bounded(Ft, left = 0, right = 1)
  if (inherits(distmat, "dist")) {
    assert_eq(length(demesize), nrow(as.matrix(distmat)),
              message = "Distance matrix must contain a row and column
                       for each deme")
  } else {
    assert_eq(length(demesize), nrow(distmat),
              message = "Distance matrix must contain a row and column
                       for each deme")
    assert_eq(nrow(distmat), ncol(distmat),
              message = "If not a distance matrix (class: 'dist'), then
              distance matrix must be square")
    assert_eq(distmat[ lower.tri(distmat) ], distmat[ upper.tri(distmat) ],
              message = "If not a distance matrix (class: 'dist'), then
                          must be a symmetric matrix")
  }


  #......................
  # simulation
  #......................
  # distance to long
  demedistlong <- broom::tidy(as.dist(distmat))
  demedistlong <- dplyr::bind_rows(demedistlong,
                                   tibble::tibble(item1 = sort(unique(c(demedistlong$item1, demedistlong$item2))),
                                                  item2 = sort(unique(c(demedistlong$item1, demedistlong$item2))),
                                                  distance = 0)
  ) %>%
    magrittr::set_colnames(c("deme1", "deme2", "geodist"))
  # individuals to long
  inds <- 1:cumsum(demesize)[length(demesize)]
  # trackers for demes
  d1 <- tibble::tibble(smpl1 = inds, deme1 = as.factor(rep(1:length(demesize), demesize)),
                       demesize1 = rep(demesize, demesize))
  d2 <- tibble::tibble(smpl2 = inds, deme2 = as.factor(rep(1:length(demesize), demesize)),
                       demesize2 = rep(demesize, demesize))
  # combinations and joins
  combinds <- tibble::as_tibble(t(combn(inds, 2)), .name_repair = "minimal") %>%
    magrittr::set_colnames(c("smpl1", "smpl2")) %>%
    dplyr::left_join(., d1, by = "smpl1") %>%
    dplyr::left_join(., d2, by = "smpl2") %>%
    dplyr::left_join(., demedistlong, by = c("deme1", "deme2"))

  # function for drawing mean IBD based iso by dist
  draw_mean_Ftibd <- function(demesize1, demesize2, geodist, rate, Ft) {
    ret <- 1/(demesize1 + demesize2) * (1 - pexp(geodist, rate)) +
      ( (1 - 1/(demesize1 + demesize2)) *  (1 - pexp(geodist, rate)) * Ft )
    return(ret)
  }
  combinds$mft <- purrr::pmap_dbl(combinds[,c("demesize1", "demesize2", "geodist")],
                                  draw_mean_Ftibd,
                                  rate = rate, Ft = Ft)

  # draw value of realized IBD from mean distribution
  draw_realized_ibdibd <- function(mft) {
    ribdibd <- discent:::rnorm_interval(mean = mft, sd = mft * (1-mft),
                                        a = 0, b = 1)
    return(ribdibd)
  }

  # tidy up and out
  combinds <- combinds %>%
    dplyr::mutate(gendist = purrr::map_dbl(mft, draw_realized_ibdibd)) %>%
    dplyr::select(c("smpl1", "smpl2", "deme1", "deme2", "gendist", "geodist"))
  return(combinds)

}
