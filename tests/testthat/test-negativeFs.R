test_that("No More negative Fs with logit", {
  ## Note: Produced negative value with seed on previous versions

  # sim data
  data("IBD_simulation_data", package = "discent")
  input <- IBD_simulation_data
  input <- input %>%
    dplyr::filter(deme1 != deme2)

  #......................
  # run discent
  #......................
  our_start_params <- rep(0.2, 3)
  names(our_start_params) <- 1:3
  our_start_params <- c(our_start_params, "m" = 1e3)
  ret <- deme_inbreeding_spcoef_vanilla(discdat = input,
                                        start_params = our_start_params,
                                        f_learningrate = 1e-3,
                                        m_learningrate = 1e-5,
                                        b1 = 0.9,
                                        b2 = 0.999,
                                        e = 1e-8,
                                        steps = 1e4,
                                        normalize_geodist = F,
                                        report_progress = T,
                                        return_verbose = T)

  # Fis all positive
  testthat::expect_gte(min(ret$Final_Fis), 0)

})
