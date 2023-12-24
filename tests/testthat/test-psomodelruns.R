test_that("PSO modelruns", {

  # sim data
  data("IBD_simulation_data", package = "discent")
  dat <- IBD_simulation_data
  # run model
  inputdisc <- dat %>%
    dplyr::filter(deme1 != deme2)
  mod <- deme_inbreeding_spcoef_pso(discdat = inputdisc,
                                    m_lowerbound = 1e-10,
                                    m_upperbound = Inf,
                                    fi_lowerbound= 1e-3,
                                    fi_upperbound = 0.3,
                                    flearn_lowerbound = 1e-10,
                                    flearn_upperbound = 1e-2,
                                    mlearn_lowerbound = 1e-15,
                                    mlearn_upperbound = 1e-8,
                                    c1 = 0.1,
                                    c2 = 0.1,
                                    w = 0.25,
                                    b1 = 0.9,
                                    b2 = 0.999,
                                    e = 1e-8,
                                    steps = 1e3,
                                    searchsteps = 1e1,
                                    swarmsize = 5,
                                    swarmsteps = 10,
                                    normalize_geodist = F,
                                    report_progress = F,
                                    return_verbose = F)
  testthat::expect_length(mod, 6)
})

