test_that("M is properly bounded", {

  #............................................................
  # simulator that is not generalizable
  #...........................................................
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
  # set seed
  set.seed(48)
  # sim data
  dat <- sim_IBDIBD(demesize = c(12,12), distmat = matrix(c(0,1e3,1e3,0), nrow = 2),
                    rate = 1e-2, Ft = 0.3)
  # start params
  our_start_params <- rep(0.2, 2)
  names(our_start_params) <- 1:2
  our_start_params <- c(our_start_params, "m" = 1e3)
  # run model w/ reasonable bounds
  inputdisc <- dat %>%
    dplyr::filter(deme1 != deme2)
  mod1 <- deme_inbreeding_spcoef(discdat = inputdisc,
                                 start_params = our_start_params,
                                 f_learningrate = 1e-3,
                                 m_learningrate = 1e-5,
                                 m_lowerbound = 999,
                                 m_upperbound = 1001,
                                 b1 = 0.9,
                                 b2 = 0.999,
                                 e = 1e-8,
                                 steps = 1e2,
                                 report_progress = TRUE,
                                 return_verbose = TRUE)
  testthat::expect_gte(min(mod1$m_run), 999)
  testthat::expect_lte(max(mod1$m_run), 1001)

  #......................
  # now show that M won't move if bounded completely
  #......................
  mod2 <- deme_inbreeding_spcoef(discdat = inputdisc,
                                start_params = our_start_params,
                                f_learningrate = 1e-3,
                                m_learningrate = 1e-5,
                                m_lowerbound = 99,
                                m_upperbound = 1e5,
                                b1 = 0.9,
                                b2 = 0.999,
                                e = 1e-8,
                                steps = 1e2,
                                report_progress = TRUE,
                                return_verbose = TRUE)
  testthat::expect_gte(min(mod2$m_run[2:length(mod2$m_run)]), 99) # offset for init param
  testthat::expect_lte(max(mod2$m_run), 1e5)
})
