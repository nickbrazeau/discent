## .................................................................................
## Purpose: IBD and IBDistance simulator that is not generalizable but provides
## tody data
## .................................................................................
#...........................................................
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
  inds <- seq_len(cumsum(demesize)[length(demesize)])
  # trackers for demes
  d1 <- tibble::tibble(smpl1 = inds, deme1 = as.factor(rep(seq_len(length(demesize)), demesize)),
                       demesize1 = rep(demesize, demesize))
  d2 <- tibble::tibble(smpl2 = inds, deme2 = as.factor(rep(seq_len(length(demesize)), demesize)),
                       demesize2 = rep(demesize, demesize))
  # combinations and joins
  combinds <- tibble::as_tibble(t(combn(inds, 2)), .name_repair = "minimal") %>%
    magrittr::set_colnames(c("smpl1", "smpl2")) %>%
    dplyr::left_join(., d1, by = "smpl1") %>%
    dplyr::left_join(., d2, by = "smpl2") %>%
    dplyr::left_join(., demedistlong, by = c("deme1", "deme2"))

  # function for drawing mean IBD based iso by dist
  draw_mean_Ftibd <- function(demesize1, demesize2, geodist, rate, Ft) {
    1/(demesize1 + demesize2) * (1 - pexp(geodist, rate)) +
      ( (1 - 1/(demesize1 + demesize2)) *  (1 - pexp(geodist, rate)) * Ft )
  }
  combinds$mft <- purrr::pmap_dbl(combinds[,c("demesize1", "demesize2", "geodist")],
                                  draw_mean_Ftibd,
                                  rate = rate, Ft = Ft)

  # draw value of realized IBD from mean distribution
  draw_realized_ibdibd <- function(mft) {
    rnorm_interval(mean = mft, sd = mft * (1-mft),
                   a = 0, b = 1)
  }

  # tidy up and out
  combinds %>%
    dplyr::mutate(gendist = purrr::map_dbl(mft, draw_realized_ibdibd)) %>%
    dplyr::select(c("smpl1", "smpl2", "deme1", "deme2", "gendist", "geodist"))


}

# Deme A and B are close but A is very far from C, while B is closer (eg a right triangle)
set.seed(48)
IBD_simulation_data <- sim_IBDIBD(demesize = c(3,3,4), distmat = matrix(c(0,500,1000,
                                                                          500,0,750,
                                                                          100,750,0), nrow = 3),
                                  rate = 1e-3, Ft = 0.3)

# drop within deme
IBD_simulation_data <- IBD_simulation_data[IBD_simulation_data$geodist > 0, ]

#............................................................
# out
#...........................................................

usethis::use_data(IBD_simulation_data, overwrite = TRUE)
