test_that("model runs works", {
  # set seed
  set.seed(48)
  # sim data
  dat <- discent::sim_IBDIBD(demesize = c(5,5), distmat = matrix(c(0,1e4,1e4,0), nrow = 2),
                             rate = 1e-3, Ft = 0.5)
  # start params
  our_start_params <- rep(0.2, 2)
  names(our_start_params) <- 1:2
  our_start_params <- c(our_start_params, "m" = 1e-3)
  # run model
  inputdisc <- dat %>%
    dplyr::filter(deme1 != deme2)
  mod <- discent::deme_inbreeding_spcoef(discdat = inputdisc,
                                         start_params = our_start_params,
                                         f_learningrate = 1e-5,
                                         m_learningrate = 1e-10,
                                         momentum = 0.9,
                                         steps = 1e2,
                                         report_progress = TRUE)
  testthat::expect_length(mod, 4)
})
