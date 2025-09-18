test_that("Fi gradient by hand", {
  # sim data
  data("IBD_simulation_data", package = "discent")
  dat <- IBD_simulation_data

  # start params
  our_start_params <- rep(0.2, 3)
  names(our_start_params) <- 1:3
  our_start_params <- c(our_start_params, "m" = 1e3)

  # manual gradient
  fgrad <- function(gendist, geodist, fi, fj, m) {
    -gendist * exp(-geodist/m) + ((fi+fj)/2) * exp(-2*geodist/m)
  }

  # tidy date and calculate gradient for F1
  input <- dat %>%
    dplyr::filter(deme1 != deme2)


  input <- input %>%
    dplyr::mutate(gendist = ifelse(gendist > 0.999, 0.999,
                                   ifelse(gendist < 0.001, 0.001,
                                          gendist))) %>% # reasonable bounds on logit
    dplyr::mutate(gendist = logit(gendist))


  # With fgrad[j] contribution, F1 gets gradient from all pairs involving deme 1
  # Calculate F1 gradient contributions from all relevant pairs
  f1retgrad <- 0

  # F1 as first deme (pairs 1-2, 1-3)
  f1_pairs_as_i <- input %>% dplyr::filter(deme1 == 1, deme2 %in% c(2,3))
  if(nrow(f1_pairs_as_i) > 0) {
    f1retgrad <- f1retgrad + sum(purrr::pmap_dbl(f1_pairs_as_i[,c("gendist", "geodist")], fgrad,
                                                  fi = logit(0.2),
                                                  fj = logit(0.2),
                                                  m = 1e3))
  }

  # F1 as second deme (pairs 2-1, 3-1) - same contribution due to symmetry
  f1_pairs_as_j <- input %>% dplyr::filter(deme2 == 1, deme1 %in% c(2,3))
  if(nrow(f1_pairs_as_j) > 0) {
    f1retgrad <- f1retgrad + sum(purrr::pmap_dbl(f1_pairs_as_j[,c("gendist", "geodist")], fgrad,
                                                  fi = logit(0.2),
                                                  fj = logit(0.2),
                                                  m = 1e3))
  }

  # now run model
  inputdisc <- dat %>%
    dplyr::filter(deme1 != deme2)
  ret <- disc(discdat = inputdisc,
              start_params = our_start_params,
              learningrate = 1e-3,
              lambda = 0,
              b1 = 0.9,
              b2 = 0.999,
              e = 1e-8,
              steps = 1e3,
              normalize_geodist = F,
              report_progress = T,
              return_verbose = T)
  # back out gradient for F1
  # NBset to 0 at first iter, so the additional adam term cancels out
  discF1 <- ret$fi_gradtraj[2,1]

  # test out
  testthat::expect_equal(f1retgrad, discF1)

})


test_that("M gradient by hand", {
  # sim data
  data("IBD_simulation_data", package = "discent")
  dat <- IBD_simulation_data
  # start params
  our_start_params <- rep(0.2, 3)
  names(our_start_params) <- 1:3
  our_start_params <- c(our_start_params, "m" = 1e3)

  # manual gradient
  mgrad <- function(gendist, geodist, fi, fj, m) {
    (-2*gendist*geodist)/(m^2) * ((fi+fj)/2) * exp(-geodist/m) +
      (2*geodist)/(m^2) * ((fi^2 + 2*fi*fj + fj^2)/4) * exp(-2*geodist/m)
  }

  # tidy date and calculate gradient for M
  input <- dat %>%
    dplyr::mutate(gendist = ifelse(gendist > 0.999, 0.999,
                                   ifelse(gendist < 0.001, 0.001,
                                          gendist))) %>% # reasonable bounds on logit
    dplyr::mutate(gendist = logit(gendist))

  # run grad by hand
  Mretgrad <- sum(purrr::pmap_dbl(input[,c("gendist", "geodist")], mgrad,
                                  fi = logit(0.2),
                                  fj = logit(0.2),
                                  m = 1e3)) # from start params

  # now run model
  inputdisc <- dat %>%
    dplyr::filter(deme1 != deme2)
  ret <- disc(discdat = inputdisc,
              start_params = our_start_params,
              learningrate = 1e-3,
              lambda = 0,
              b1 = 0.9,
              b2 = 0.999,
              e = 1e-8,
              steps = 1e3,
              normalize_geodist = F,
              report_progress = T,
              return_verbose = T)
  # back out gradient for M
  # NB  set to 0 at first iter, so the additional adam term cancels out
  discM <- ret$m_gradtraj[2]

  # test out
  testthat::expect_equal(Mretgrad, discM)

})

