test_that("Fi gradient by hand", {
  # set seed
  set.seed(48)
  # sim data
  dat <- discent::sim_IBDIBD(demesize = c(2,2), distmat = matrix(c(0,100000,100000,0), nrow = 2),
                      rate = 1e-2, Ft = 0.3)
  # start params
  our_start_params <- rep(0.2, 2)
  names(our_start_params) <- 1:2
  our_start_params <- c(our_start_params, "m" = 1e-5)

  # manual gradient
  fgrad <- function(gendist, geodist, fi, fj, m) {
    -gendist * exp(-geodist * m) + ((fi+fj)/2)*exp(2*geodist*m)
  }

  # tidy date and calculate gradient for F1
  input <- dat %>%
    dplyr::filter(locat1 != locat2) %>%
    dplyr::filter(locat1 == 1 | locat2 == 1)
  f1retgrad <- sum(purrr::pmap_dbl(input[,c("gendist", "geodist")], fgrad,
              fi = 0.2, fj = 0.2, m = 1e-5)) # from start params

  # now run model
  inputdisc <- dat %>%
    dplyr::filter(locat1 != locat2)
  ret <- discent::deme_inbreeding_spcoef(K_gendist_geodist = inputdisc,
                                         start_params = our_start_params,
                                         f_learningrate = 1e-5,
                                         m_learningrate = 1e-10,
                                         steps = 1e4,
                                         report_progress = T,
                                         return_verbose = T)
  # back out gradient for F1
  discF1 <- ret$fi_update[2,1]/1e-5

  # test out
  testthat::expect_equal(f1retgrad, discF1)

})


test_that("M gradient by hand", {
  #TODO
})
