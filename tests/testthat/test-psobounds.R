test_that("PSO respects bounds", {

  # sim data
  data("IBD_simulation_data", package = "discent")
  dat <- IBD_simulation_data
  # run model
  inputdisc <- dat %>%
    dplyr::filter(deme1 != deme2)
  mod <- deme_inbreeding_spcoef_pso(discdat = inputdisc,
                                    m_lowerbound = 100,
                                    m_upperbound = 101,
                                    fi_lowerbound= 0.1,
                                    fi_upperbound = 0.11,
                                    flearn_lowerbound = 1e-3,
                                    flearn_upperbound = 1e-1,
                                    mlearn_lowerbound = 1e-3,
                                    mlearn_upperbound = 100,
                                    c1 = 0.1,
                                    c2 = 0.1,
                                    w = 0.25,
                                    b1 = 0.9,
                                    b2 = 0.999,
                                    e = 1e-8,
                                    steps = 1e3,
                                    searchsteps = 1e2,
                                    swarmsize = 5,
                                    swarmsteps = 10,
                                    normalize_geodist = F,
                                    report_sd_progress = F,
                                    return_verbose = T)
  #......................
  # check that bounds are respect
  #......................
  testthat::expect_gte(min(mod$m_run), 100)
  testthat::expect_lte(max(mod$m_run), 101)
  testthat::expect_gte(min(mod$fi_run[,1]), 0.1)
  testthat::expect_lte(max(mod$fi_run[,1]), 0.11)
  testthat::expect_gte(min(mod$fi_run[,2]), 0.1)
  testthat::expect_lte(max(mod$fi_run[,2]), 0.11)
  testthat::expect_gte(min(mod$fi_run[,3]), 0.1)
  testthat::expect_lte(max(mod$fi_run[,3]), 0.11)


})

