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
    dplyr::filter(deme1 != deme2) %>%
    dplyr::filter(deme1 == 1 | deme2 == 1)


  input <- input %>%
    dplyr::mutate(gendist = ifelse(gendist > 0.999, 0.999,
                                   ifelse(gendist < 0.001, 0.001,
                                          gendist))) %>% # reasonable bounds on logit
    dplyr::mutate(gendist = logit(gendist))


  f1retgrad <- sum(purrr::pmap_dbl(input[,c("gendist", "geodist")], fgrad,
                                   fi = logit(0.2),
                                   fj = logit(0.2),
                                   m = 1e3)) # from start params

  # now run model
  inputdisc <- dat %>%
    dplyr::filter(deme1 != deme2)
  ret <- deme_inbreeding_spcoef(discdat = inputdisc,
                                start_params = our_start_params,
                                f_learningrate = 1e-3,
                                m_learningrate = 1e-5,
                                b1 = 0.9,
                                b2 = 0.999,
                                e = 1e-8,
                                steps = 1e3,
                                standardize_geodist = F,
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
  # note, need to expand matrix because comparisons are made for each deme - e.g. full distance matrix (not just combinations: lower or upper tri)
  datexpand <- dat
  colnames(datexpand) <- c("smpl2", "smpl1", "deme2", "deme1", "gendist", "geodist")
  input <- dplyr::bind_rows(dat, datexpand) %>%
    dplyr::filter(deme1 != deme2)
  # logit transform
  input <- input %>%
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
  ret <- deme_inbreeding_spcoef(discdat = inputdisc,
                                start_params = our_start_params,
                                f_learningrate = 1e-3,
                                m_learningrate = 1e-5,
                                b1 = 0.9,
                                b2 = 0.999,
                                e = 1e-8,
                                steps = 1e3,
                                standardize_geodist = F,
                                report_progress = T,
                                return_verbose = T)
  # back out gradient for M
  # NB  set to 0 at first iter, so the additional adam term cancels out
  discM <- ret$m_gradtraj[2]

  # test out
  testthat::expect_equal(Mretgrad, discM)

})

