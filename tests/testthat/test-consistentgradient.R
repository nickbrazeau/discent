test_that("model has deterministic results from same start", {

  # sim data
  data("IBD_simulation_data", package = "discent")
  inputdisc <- IBD_simulation_data
  # start params
  our_start_params <- rep(0.2, 3)
  names(our_start_params) <- 1:3
  our_start_params <- c(our_start_params, "m" = 1e3)

  # given same start parameters and data and iters
  # model should always follow the same gradient/trajectory
  mod1 <- deme_inbreeding_spcoef_vanilla(discdat = inputdisc,
                                 start_params = our_start_params,
                                 f_learningrate = 1e-5,
                                 m_learningrate = 1e-1,
                                 b1 = 0.9,
                                 b2 = 0.999,
                                 e = 1e-8,
                                 steps = 1e4,
                                 normalize_geodist = F,
                                 report_progress = F,
                                 return_verbose = F)

  mod2 <- deme_inbreeding_spcoef_vanilla(discdat = inputdisc,
                                 start_params = our_start_params,
                                 f_learningrate = 1e-5,
                                 m_learningrate = 1e-1,
                                 b1 = 0.9,
                                 b2 = 0.999,
                                 e = 1e-8,
                                 steps = 1e4,
                                 normalize_geodist = F,
                                 report_progress = F,
                                 return_verbose = F)

  mod3 <- deme_inbreeding_spcoef_vanilla(discdat = inputdisc,
                                 start_params = our_start_params,
                                 f_learningrate = 1e-5,
                                 m_learningrate = 1e-1,
                                 b1 = 0.9,
                                 b2 = 0.999,
                                 e = 1e-8,
                                 steps = 1e4,
                                 normalize_geodist = F,
                                 report_progress = F,
                                 return_verbose = F)

  mod4 <- deme_inbreeding_spcoef_vanilla(discdat = inputdisc,
                                 start_params = our_start_params,
                                 f_learningrate = 1e-5,
                                 m_learningrate = 1e-1,
                                 b1 = 0.9,
                                 b2 = 0.999,
                                 e = 1e-8,
                                 steps = 5e4, # different number of steps!!
                                 normalize_geodist = F,
                                 report_progress = F,
                                 return_verbose = F)

  testthat::expect_equal(mod1,mod2)
  testthat::expect_equal(mod1,mod3)
  testthat::expect_equal(mod2,mod2) # trans, so know true, but for complete
  testthat::expect_false(identical(mod2,mod4)) # more steps, different results
})


