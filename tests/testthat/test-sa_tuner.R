test_that("Simulated Annealer Tuner works as expected", {
  #............................................................
  # simulator that is not generalizable
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
    # #......................
    # # assertions
    # #......................
    # assert_numeric(demesize)
    # assert_numeric(distmat)
    # assert_numeric(rate)
    # assert_numeric(Ft)
    # assert_single_bounded(Ft, left = 0, right = 1)
    # if (inherits(distmat, "dist")) {
    #   assert_eq(length(demesize), nrow(as.matrix(distmat)),
    #             message = "Distance matrix must contain a row and column
    #                    for each deme")
    # } else {
    #   assert_eq(length(demesize), nrow(distmat),
    #             message = "Distance matrix must contain a row and column
    #                    for each deme")
    #   assert_eq(nrow(distmat), ncol(distmat),
    #             message = "If not a distance matrix (class: 'dist'), then
    #           distance matrix must be square")
    #   assertthat::are_equal(distmat[ lower.tri(distmat) ], distmat[ upper.tri(distmat) ],
    #                         message = "If not a distance matrix (class: 'dist'), then
    #                       must be a symmetric matrix")
    # }


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
      ribdibd <- rnorm_interval(mean = mft, sd = mft * (1-mft),
                                a = 0, b = 1)
      return(ribdibd)
    }

    # tidy up and out
    combinds <- combinds %>%
      dplyr::mutate(gendist = purrr::map_dbl(mft, draw_realized_ibdibd)) %>%
      dplyr::select(c("smpl1", "smpl2", "deme1", "deme2", "gendist", "geodist"))
    return(combinds)

  }


  #............................................................
  # actual test
  #...........................................................
  # sim dat
  distmatsim <- matrix(runif(9, min = 0, max = 1e2), nrow = 3)
  diag(distmatsim) <- 0
  distmatsim[lower.tri(distmatsim)] <- t(distmatsim)[lower.tri(distmatsim)]
  dat <- sim_IBDIBD(demesize = c(2,2,2), distmat = distmatsim,
                    rate = 1e-2, Ft = 0.3)
  dat <- dat %>%
    dplyr::filter(deme1 != deme2)
  # start params
  our_start_params <- rep(0.2, 3)
  names(our_start_params) <- 1:3
  our_start_params <- c(our_start_params, "m" = 1e-5)

  #......................
  # simulated annealer
  #......................
  ret <- discent::find_grad_params(discdat = dat,
                                   initstart_params = our_start_params,
                                   initf_learningrate = 1e-5,
                                   initm_learningrate = 1e-10,
                                   momentum = 0.9,
                                   discsteps = 1e1,

                                   initTemp = 1,
                                   annealstep = 1e1,
                                   SLratio = 0.8,
                                   FMSratio = 0.8,
                                   FMLratio = 0.5,
                                   demeSwitchSize = 3,

                                   fstartmin = 0.1,
                                   fstartmax = 0.9,
                                   fstartsteps = 100,

                                   mstartmin = 1e-6,
                                   mstartmax = 1,
                                   mstartsteps = 100,

                                   flearnmin = 1e-15,
                                   flearnmax = 1e-2,
                                   flearnsteps = 100,

                                   mlearnmin = 1e-20,
                                   mlearnmax = 1e-8,
                                   mlearnsteps = 100)

  testthat::expect_length(ret, 2)
  testthat::expect_vector(ret$optimalParams$start_params)
  testthat::expect_gte(min(ret$optimalParams$start_params[!grepl("^m$", names(ret$optimalParams$start_params))]),
                       0)
  testthat::expect_type(ret$optimalParams$start_params[grepl("^m$", names(ret$optimalParams$start_params))],
                        "double")
  testthat::expect_gte(ret$optimalParams$f_learn, 0)
  testthat::expect_gte(ret$optimalParams$m_learn, 0)

})
