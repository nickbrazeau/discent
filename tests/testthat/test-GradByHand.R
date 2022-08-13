test_that("Fi gradient by hand", {
  # set seed
  set.seed(48)
  # sim data
  dat <- discent::sim_IBDIBD(demesize = c(12,12), distmat = matrix(c(0,1e3,1e3,0), nrow = 2),
                             rate = 1e-2, Ft = 0.3)
  # start params
  our_start_params <- rep(0.2, 2)
  names(our_start_params) <- 1:2
  our_start_params <- c(our_start_params, "m" = 1e-5)

  # manual gradient
  fgrad <- function(gendist, geodist, fi, fj, m) {
    -gendist * exp(-geodist * m) + ((fi+fj)/2)*exp(-2*geodist*m)
  }

  # tidy date and calculate gradient for F1
  input <- dat %>%
    dplyr::filter(deme1 != deme2) %>%
    dplyr::filter(deme1 == 1 | deme2 == 1)
  f1retgrad <- sum(purrr::pmap_dbl(input[,c("gendist", "geodist")], fgrad,
                                   fi = 0.2, fj = 0.2, m = 1e-5)) # from start params

  # now run model
  inputdisc <- dat %>%
    dplyr::filter(deme1 != deme2)
  ret <- discent::deme_inbreeding_spcoef(discdat = inputdisc,
                                         start_params = our_start_params,
                                         f_learningrate = 1e-5,
                                         m_learningrate = 1e-10,
                                         steps = 1e3,
                                         report_progress = T,
                                         return_verbose = T)
  # back out gradient for F1
  discF1 <- ret$fi_update[2,1]/1e-5

  # test out
  testthat::expect_equal(f1retgrad, discF1)

})


test_that("M gradient by hand", {
  # set seed
  set.seed(48)
  # sim data
  dat <- discent::sim_IBDIBD(demesize = c(2,2), distmat = matrix(c(0,1e3,1e3,0), nrow = 2),
                             rate = 1e-2, Ft = 0.3)
  # start params
  our_start_params <- rep(0.2, 2)
  names(our_start_params) <- 1:2
  our_start_params <- c(our_start_params, "m" = 1e-5)

  # manual gradient
  mgrad <- function(gendist, geodist, fi, fj, m) {
    2*gendist*geodist * ((fi+fj)/2)*exp(-geodist*m) -
      2*geodist*((fi^2 + 2*fi*fj + fj^2)/4)*exp(-2*geodist*m)
  }

  # tidy date and calculate gradient for M
  # note, need to expand matrix because comparisons are made for each deme - e.g. full distance matrix (not just combinations: lower or upper tri)
  datexpand <- dat
  colnames(datexpand) <- c("smpl2", "smpl1", "deme2", "deme1", "gendist", "geodist")
  input <- dplyr::bind_rows(dat, datexpand) %>%
    dplyr::filter(deme1 != deme2)
  Mretgrad <- sum(purrr::pmap_dbl(input[,c("gendist", "geodist")], mgrad,
                                  fi = 0.2, fj = 0.2, m = 1e-5)) # from start params

  # now run model
  inputdisc <- dat %>%
    dplyr::filter(deme1 != deme2)
  ret <- discent::deme_inbreeding_spcoef(discdat = inputdisc,
                                         start_params = our_start_params,
                                         f_learningrate = 1e-5,
                                         m_learningrate = 1e-10,
                                         steps = 1e3,
                                         report_progress = T,
                                         return_verbose = T)
  # back out gradient for M
  discM <- ret$m_update[2]/1e-10

  # test out
  testthat::expect_equal(Mretgrad, discM)

})
