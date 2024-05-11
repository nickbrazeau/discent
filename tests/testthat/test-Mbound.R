test_that("M is properly bounded", {
  # sim data
  data("IBD_simulation_data", package = "discent")
  dat <- IBD_simulation_data
  # start params
  our_start_params <- rep(0.2, 3)
  names(our_start_params) <- 1:3
  our_start_params <- c(our_start_params, "m" = 1e3)
  # run model w/ EXTREMELY TIGHT bounds
  inputdisc <- dat %>%
    dplyr::filter(deme1 != deme2)
  mod1 <- deme_inbreeding_spcoef_vanilla(discdat = inputdisc,
                                         start_params = our_start_params,
                                         f_learningrate = 1e-3,
                                         m_learningrate = 1e-1,
                                         m_lowerbound = 999,
                                         m_upperbound = 1001,
                                         b1 = 0.9,
                                         b2 = 0.999,
                                         e = 1e-8,
                                         steps = 1e2,
                                         report_progress = TRUE,
                                         normalize_geodist = FALSE,
                                         return_verbose = TRUE)
  testthat::expect_gte(min(mod1$m_run), 999)
  testthat::expect_lte(max(mod1$m_run), 1001)

})
