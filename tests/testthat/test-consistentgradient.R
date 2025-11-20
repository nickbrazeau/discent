test_that("cost decreases monotonically", {

  # sim data
  data("IBD_simulation_data", package = "discent")
  inputdisc <- IBD_simulation_data
  # start params
  our_start_params <- rep(0.2, 3)
  names(our_start_params) <- 1:3
  our_start_params <- c(our_start_params, "m" = 1e3)

  # given same start parameters and data and iters
  # model should always follow the same gradient/trajectory
  mod1 <- disc(discdat = inputdisc,
               start_params = our_start_params,
               learningrate = 1e-5,
               b1 = 0.9,
               b2 = 0.999,
               e = 1e-8,
               steps = 1e4,
               report_progress = F,
               return_verbose = F)
  # monotonic decrease general
  testthat::expect_lt(mod1$cost[1e2], (mod1$cost[1]))
  testthat::expect_lt(mod1$cost[1e3], (mod1$cost[1e2]))
  testthat::expect_lt(mod1$cost[1e4], (mod1$cost[1e3]))
})

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
  mod1 <- disc(discdat = inputdisc,
               start_params = our_start_params,
               learningrate = 1e-5,
               b1 = 0.9,
               b2 = 0.999,
               e = 1e-8,
               steps = 1e4,
               report_progress = F,
               return_verbose = F)

  mod2 <- disc(discdat = inputdisc,
               start_params = our_start_params,
               learningrate = 1e-5,
               b1 = 0.9,
               b2 = 0.999,
               e = 1e-8,
               steps = 1e4,
               report_progress = F,
               return_verbose = F)

  mod3 <- disc(discdat = inputdisc,
               start_params = our_start_params,
               learningrate = 1e-5,
               b1 = 0.9,
               b2 = 0.999,
               e = 1e-8,
               steps = 1e4,
               report_progress = F,
               return_verbose = F)

  mod4 <- disc(discdat = inputdisc,
               start_params = our_start_params,
               learningrate = 1e-5,
               b1 = 0.9,
               b2 = 0.999,
               e = 1e-8,
               steps = 5e4, # different number of steps!!
               report_progress = F,
               return_verbose = F)

  testthat::expect_equal(mod1,mod2)
  testthat::expect_equal(mod1,mod3)
  testthat::expect_equal(mod2,mod2) # trans, so know true, but for complete
  testthat::expect_false(identical(mod2,mod4)) # more steps, different results
})


