test_that("PSO model output is consistent w/ Fstart velocity calculations", {

  # sim data
  data("IBD_simulation_data", package = "discent")
  dat <- IBD_simulation_data
  # run model
  inputdisc <- dat %>%
    dplyr::filter(deme1 != deme2)
  mod0 <- deme_inbreeding_spcoef_pso(discdat = inputdisc,
                                    m_lowerbound = 0.1,
                                    m_upperbound = 1e3,
                                    fi_lowerinit = 1e-3,
                                    fi_upperinit = 0.3,
                                    learn_lowerinit = 1e-10,
                                    learn_upperinit = 1e-2,
                                    c1 = 2,
                                    c2 = 2,
                                    w = 0.8,
                                    b1 = 0.9,
                                    b2 = 0.999,
                                    e = 1e-8,
                                    swarmsize = 5,
                                    swarmmoves = 10,
                                    particlesteps = 1e1,
                                    finalsteps = 1e3,
                                    normalize_geodist = F,
                                    return_verbose = T)

  Fstart <- mod0$swarm %>%
    dplyr::filter(particle_int ==1) %>%
    dplyr::select(dplyr::contains("Fstart"))
  chck <- Fstart$particle_Poscurr_Fstart + dplyr::lead(Fstart$particle_Veloccurr_Fstart)
  testthat::expect_equal(chck[1:9], Fstart$particle_Poscurr_Fstart[2:10])
})

test_that("PSO model output is consistent w/ Mstart velocity calculations", {

  # sim data
  data("IBD_simulation_data", package = "discent")
  dat <- IBD_simulation_data
  # run model
  inputdisc <- dat %>%
    dplyr::filter(deme1 != deme2)
  mod0 <- deme_inbreeding_spcoef_pso(discdat = inputdisc,
                                     m_lowerbound = 0.1,
                                     m_upperbound = 1e3,
                                     fi_lowerinit = 1e-3,
                                     fi_upperinit = 0.3,
                                     flearn_lowerinit = 1e-10,
                                     flearn_upperinit = 1e-2,
                                     mlearn_lowerinit = 1e-5,
                                     mlearn_upperinit = 1e-1,
                                     c1 = 2,
                                     c2 = 2,
                                     w = 0.8,
                                     b1 = 0.9,
                                     b2 = 0.999,
                                     e = 1e-8,
                                     swarmsize = 5,
                                     swarmmoves = 10,
                                     particlesteps = 1e1,
                                     finalsteps = 1e3,
                                     normalize_geodist = F,
                                     return_verbose = T)

  Mstart <- mod0$swarm %>%
    dplyr::filter(particle_int ==1) %>%
    dplyr::select(dplyr::contains("Mstart"))
  chck1 <- Mstart$particle_Poscurr_Mstart + dplyr::lead(Mstart$particle_Veloccurr_Mstart)
  chck1 <- chck1[1:9]
  chck2 <- Mstart$particle_Poscurr_Mstart[2:10]
  # need to remove ones with bounds bc those should be right
  mknas <- which(chck2 %in% c(0.1, 1e3))
  chck1[mknas] <- NA
  chck2[mknas] <- NA
  testthat::expect_equal(chck1, chck2)
})


test_that("PSO model output is consistent w/ Flearn velocity calculations", {

  # sim data
  data("IBD_simulation_data", package = "discent")
  dat <- IBD_simulation_data
  # run model
  inputdisc <- dat %>%
    dplyr::filter(deme1 != deme2)
  mod0 <- deme_inbreeding_spcoef_pso(discdat = inputdisc,
                                     m_lowerbound = 0.1,
                                     m_upperbound = 1e3,
                                     fi_lowerinit = 1e-3,
                                     fi_upperinit = 0.3,
                                     flearn_lowerinit = 1e-10,
                                     flearn_upperinit = 1e-2,
                                     mlearn_lowerinit = 1e-5,
                                     mlearn_upperinit = 1e-1,
                                     c1 = 2,
                                     c2 = 2,
                                     w = 0.8,
                                     b1 = 0.9,
                                     b2 = 0.999,
                                     e = 1e-8,
                                     swarmsize = 5,
                                     swarmmoves = 10,
                                     particlesteps = 1e1,
                                     finalsteps = 1e3,
                                     normalize_geodist = F,
                                     return_verbose = T)

  Flearn <- mod0$swarm %>%
    dplyr::filter(particle_int ==1) %>%
    dplyr::select(dplyr::contains("Flearn"))
  chck <- Flearn$particle_Poscurr_Flearn + dplyr::lead(Flearn$particle_Veloccurr_Flearn)
  testthat::expect_equal(chck[1:9], Flearn$particle_Poscurr_Flearn[2:10])
})


test_that("PSO model output is consistent w/ Mlearn velocity calculations", {

  # sim data
  data("IBD_simulation_data", package = "discent")
  dat <- IBD_simulation_data
  # run model
  inputdisc <- dat %>%
    dplyr::filter(deme1 != deme2)
  mod0 <- deme_inbreeding_spcoef_pso(discdat = inputdisc,
                                     m_lowerbound = 0.1,
                                     m_upperbound = 1e3,
                                     fi_lowerinit = 1e-3,
                                     fi_upperinit = 0.3,
                                     flearn_lowerinit = 1e-10,
                                     flearn_upperinit = 1e-2,
                                     mlearn_lowerinit = 1e-5,
                                     mlearn_upperinit = 1e-1,
                                     c1 = 2,
                                     c2 = 2,
                                     w = 0.8,
                                     b1 = 0.9,
                                     b2 = 0.999,
                                     e = 1e-8,
                                     swarmsize = 5,
                                     swarmmoves = 10,
                                     particlesteps = 1e1,
                                     finalsteps = 1e3,
                                     normalize_geodist = F,
                                     return_verbose = T)

  Mlearn <- mod0$swarm %>%
    dplyr::filter(particle_int ==1) %>%
    dplyr::select(dplyr::contains("Mlearn"))
  chck <- Mlearn$particle_Poscurr_Mlearn + dplyr::lead(Mlearn$particle_Veloccurr_Mlearn)
  testthat::expect_equal(chck[1:9], Mlearn$particle_Poscurr_Mlearn[2:10])
})




test_that("PSO model position doesn't move if velocity zero", {

  # sim data
  data("IBD_simulation_data", package = "discent")
  dat <- IBD_simulation_data
  # run model
  inputdisc <- dat %>%
    dplyr::filter(deme1 != deme2)
  mod0 <- deme_inbreeding_spcoef_pso(discdat = inputdisc,
                                     m_lowerbound = 0.1,
                                     m_upperbound = 1e3,
                                     fi_lowerinit = 1e-3,
                                     fi_upperinit = 0.3,
                                     flearn_lowerinit = 1e-10,
                                     flearn_upperinit = 1e-2,
                                     mlearn_lowerinit = 1e-5,
                                     mlearn_upperinit = 1e-1,
                                     c1 = 0,
                                     c2 = 0,
                                     w = 0,
                                     b1 = 0.9,
                                     b2 = 0.999,
                                     e = 1e-8,
                                     swarmsize = 5,
                                     swarmmoves = 10,
                                     particlesteps = 1e1,
                                     finalsteps = 1e3,
                                     normalize_geodist = F,
                                     return_verbose = T)

  testthat::expect_equal(length(unique(mod0$swarm$particle_Posbest_Fstart)), 5) # no velocity, so it is just initial start params for the particles
})


test_that("PSO model final run respects migration bounds in swarm run", {

  # sim data
  data("IBD_simulation_data", package = "discent")
  dat <- IBD_simulation_data
  # run model
  inputdisc <- dat %>%
    dplyr::filter(deme1 != deme2)
  mod1 <- deme_inbreeding_spcoef_pso(discdat = inputdisc,
                                     m_lowerbound = 999,
                                     m_upperbound = 1001,
                                     fi_lowerinit = 1e-3,
                                     fi_upperinit = 0.3,
                                     flearn_lowerinit = 1e-10,
                                     flearn_upperinit = 1e-2,
                                     mlearn_lowerinit = 1e-5,
                                     mlearn_upperinit = 1e-1,
                                     c1 = 0.1,
                                     c2 = 0.1,
                                     w = 0.25,
                                     b1 = 0.9,
                                     b2 = 0.999,
                                     e = 1e-8,
                                     swarmsize = 5,
                                     swarmmoves = 10,
                                     particlesteps = 1e1,
                                     finalsteps = 1e3,
                                     normalize_geodist = F,
                                     return_verbose = T)
  testthat::expect_gte(min(mod1$swarm$particle_Posbest_Mstart), 999)
  testthat::expect_lte(max(mod1$swarm$particle_Posbest_Mstart), 1001)
})



test_that("PSO model final run respects migration bounds in final run", {

  # sim data
  data("IBD_simulation_data", package = "discent")
  dat <- IBD_simulation_data
  # run model
  inputdisc <- dat %>%
    dplyr::filter(deme1 != deme2)
  mod1 <- deme_inbreeding_spcoef_pso(discdat = inputdisc,
                                    m_lowerbound = 999,
                                    m_upperbound = 1001,
                                    fi_lowerinit = 1e-3,
                                    fi_upperinit = 0.3,
                                    flearn_lowerinit = 1e-10,
                                    flearn_upperinit = 1e-2,
                                    mlearn_lowerinit = 1e-5,
                                    mlearn_upperinit = 1e-1,
                                    c1 = 0.1,
                                    c2 = 0.1,
                                    w = 0.25,
                                    b1 = 0.9,
                                    b2 = 0.999,
                                    e = 1e-8,
                                    swarmsize = 5,
                                    swarmmoves = 10,
                                    particlesteps = 1e1,
                                    finalsteps = 1e3,
                                    normalize_geodist = F,
                                    return_verbose = T)
  testthat::expect_gte(min(mod1$m_run), 999)
  testthat::expect_lte(max(mod1$m_run), 1001)
})

