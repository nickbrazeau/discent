test_that("model runs", {

  # sim data
  data("IBD_simulation_data", package = "discent")
  dat <- IBD_simulation_data
  # start params
  our_start_params <- rep(0.2, 3)
  names(our_start_params) <- 1:3
  our_start_params <- c(our_start_params, "m" = 3)
  # run model
  inputdisc <- dat %>%
    dplyr::filter(deme1 != deme2)
  mod <- deme_inbreeding_spcoef(discdat = inputdisc,
                                start_params = our_start_params,
                                f_learningrate = 1e-5,
                                m_learningrate = 1e-1,
                                b1 = 0.9,
                                b2 = 0.999,
                                e = 1e-8,
                                steps = 1e3,
                                normalize_geodist = T,
                                report_progress = T,
                                return_verbose = F)
  testthat::expect_length(mod, 4)
})

