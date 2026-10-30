# Extracted from test-diagnostics.R:25

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "discent", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
data("IBD_simulation_data", package = "discent")
inputdisc <- IBD_simulation_data
our_start_params <- rep(0.2, 3)
names(our_start_params) <- 1:3
our_start_params <- c(our_start_params, "m" = 1e3)
mod1 <- disc(discdat = inputdisc,
               start_params = our_start_params,
               learningrate = 1e-4,
               b1 = 0.9,
               b2 = 0.999,
               e = 1e-8,
               steps = 1e5,
               report_progress = F,
               return_verbose = F)
testthat::expect_equal(mod1$Hessian, t(mod1$Hessian), tolerance=1e-10)
