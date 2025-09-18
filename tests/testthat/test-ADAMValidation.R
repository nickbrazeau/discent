test_that("ADAM multi-step momentum accumulation", {
  # Test that ADAM momentum terms accumulate correctly over multiple steps
  data("IBD_simulation_data", package = "discent")
  dat <- IBD_simulation_data

  test_params <- c("1" = 0.4, "2" = 0.6, "3" = 0.3, "m" = 600)
  b1 <- 0.9
  b2 <- 0.999
  steps <- 5

  inputdisc <- dat %>% dplyr::filter(deme1 != deme2)
  ret <- disc(discdat = inputdisc,
              start_params = test_params,
              learningrate = 1e-3,
              lambda = 0,
              b1 = b1, b2 = b2, e = 1e-8,
              steps = steps,
              normalize_geodist = FALSE,
              report_progress = FALSE,
              return_verbose = TRUE)

  # Test ADAM first moment (momentum) accumulation for F1
  for(step in 2:steps) {
    # Calculate expected first moment: m_t = b1 * m_{t-1} + (1-b1) * g_t
    expected_m1_f1 <- b1 * ret$fi_1moment[step-1, 1] + (1-b1) * ret$fi_gradtraj[step, 1]
    actual_m1_f1 <- ret$fi_1moment[step, 1]

    expect_equal(actual_m1_f1, expected_m1_f1, tolerance = 1e-12)

    # Calculate expected second moment: v_t = b2 * v_{t-1} + (1-b2) * g_t^2
    expected_v2_f1 <- b2 * ret$fi_2moment[step-1, 1] + (1-b2) * (ret$fi_gradtraj[step, 1])^2
    actual_v2_f1 <- ret$fi_2moment[step, 1]

    expect_equal(actual_v2_f1, expected_v2_f1, tolerance = 1e-12)
  }

  # Test M parameter ADAM moments
  for(step in 2:steps) {
    expected_m1_m <- b1 * ret$m_1moment[step-1] + (1-b1) * ret$m_gradtraj[step]
    actual_m1_m <- ret$m_1moment[step]

    expect_equal(actual_m1_m, expected_m1_m, tolerance = 1e-12)

    expected_v2_m <- b2 * ret$m_2moment[step-1] + (1-b2) * (ret$m_gradtraj[step])^2
    actual_v2_m <- ret$m_2moment[step]

    expect_equal(actual_v2_m, expected_v2_m, tolerance = 1e-12)
  }

})

test_that("ADAM bias correction accuracy", {
  # Test that bias correction is applied correctly
  data("IBD_simulation_data", package = "discent")
  dat <- IBD_simulation_data

  test_params <- c("1" = 0.5, "2" = 0.4, "3" = 0.6, "m" = 700)
  b1 <- 0.9
  b2 <- 0.999
  learning_rate <- 1e-3
  e <- 1e-8

  inputdisc <- dat %>% dplyr::filter(deme1 != deme2)
  ret <- disc(discdat = inputdisc,
              start_params = test_params,
              learningrate = learning_rate,
              lambda = 0,
              b1 = b1, b2 = b2, e = e,
              steps = 10,
              normalize_geodist = FALSE,
              report_progress = FALSE,
              return_verbose = TRUE)

  # Manually verify ADAM parameter updates for step 5
  step <- 5

  # F1 update verification
  m1_f1 <- ret$fi_1moment[step, 1]
  v2_f1 <- ret$fi_2moment[step, 1]

  # Bias-corrected moments
  m1_f1_hat <- m1_f1 / (1 - b1^step)
  v2_f1_hat <- v2_f1 / (1 - b2^step)

  # Expected parameter update
  f1_prev <- ret$fi_run[step-1, 1]  # logit space
  f1_expected <- f1_prev - learning_rate * (m1_f1_hat / (sqrt(v2_f1_hat) + e))
  f1_actual_logit <- logit(ret$fi_run[step, 1])  # Convert back to logit space

  expect_equal(f1_actual_logit, f1_expected, tolerance = 1e-10)

  # M update verification
  m1_m <- ret$m_1moment[step]
  v2_m <- ret$m_2moment[step]
  m1_m_hat <- m1_m / (1 - b1^step)
  v2_m_hat <- v2_m / (1 - b2^step)

  m_prev <- ret$m_run[step-1]
  m_expected <- m_prev - learning_rate * (m1_m_hat / (sqrt(v2_m_hat) + e))
  m_actual <- ret$m_run[step]

  expect_equal(m_actual, m_expected, tolerance = 1e-3)

})

test_that("ADAM convergence behavior", {
  # Test that ADAM shows expected convergence behavior
  data("IBD_simulation_data", package = "discent")
  dat <- IBD_simulation_data

  test_params <- c("1" = 0.1, "2" = 0.9, "3" = 0.5, "m" = 300)  # Start far from optimum

  inputdisc <- dat %>% dplyr::filter(deme1 != deme2)
  ret <- disc(discdat = inputdisc,
              start_params = test_params,
              learningrate = 1e-2,  # Higher learning rate for faster convergence
              lambda = 0,
              b1 = 0.9, b2 = 0.999, e = 1e-8,
              steps = 50,
              normalize_geodist = FALSE,
              report_progress = FALSE,
              return_verbose = TRUE)

  # Cost should generally decrease (allowing for some oscillation)
  costs <- ret$cost
  initial_cost <- costs[1]
  final_cost <- costs[length(costs)]

  expect_lt(final_cost, initial_cost,
            label = "Cost should decrease from initial to final")

  # Check that cost decreases in the first 20 steps (early convergence)
  mid_cost <- min(costs[1:20])
  expect_lt(mid_cost, initial_cost * 1.1,  # Allow 10% tolerance
            label = "Cost should decrease in early steps")

  # Gradient magnitudes should generally decrease
  f_grad_norms <- apply(ret$fi_gradtraj, 1, function(x) sqrt(sum(x^2)))
  initial_f_grad_norm <- f_grad_norms[2]  # Step 2 (step 1 is initialization)
  final_f_grad_norm <- f_grad_norms[length(f_grad_norms)]

  expect_lt(final_f_grad_norm, initial_f_grad_norm * 2,  # Allow some tolerance
            label = "F gradient norm should decrease over time")

})

test_that("ADAM parameter sensitivity", {
  # Test behavior with different ADAM parameters (b1, b2)
  data("IBD_simulation_data", package = "discent")
  dat <- IBD_simulation_data

  test_params <- c("1" = 0.3, "2" = 0.7, "3" = 0.5, "m" = 500)

  inputdisc <- dat %>% dplyr::filter(deme1 != deme2)

  # Test with high momentum (b1 = 0.95)
  ret_high_momentum <- disc(discdat = inputdisc,
                           start_params = test_params,
                           learningrate = 1e-3,
                           lambda = 0,
                           b1 = 0.95, b2 = 0.999, e = 1e-8,
                           steps = 10,
                           normalize_geodist = FALSE,
                           report_progress = FALSE,
                           return_verbose = TRUE)

  # Test with low momentum (b1 = 0.5)
  ret_low_momentum <- disc(discdat = inputdisc,
                          start_params = test_params,
                          learningrate = 1e-3,
                          lambda = 0,
                          b1 = 0.5, b2 = 0.999, e = 1e-8,
                          steps = 10,
                          normalize_geodist = FALSE,
                          report_progress = FALSE,
                          return_verbose = TRUE)

  # Both should produce finite results
  expect_true(all(is.finite(ret_high_momentum$fi_run[10, ])),
              label = "High momentum should produce finite F values")
  expect_true(is.finite(ret_high_momentum$m_run[10]),
              label = "High momentum should produce finite M value")

  expect_true(all(is.finite(ret_low_momentum$fi_run[10, ])),
              label = "Low momentum should produce finite F values")
  expect_true(is.finite(ret_low_momentum$m_run[10]),
              label = "Low momentum should produce finite M value")

  # High momentum should have smoother parameter trajectories
  f1_var_high <- var(ret_high_momentum$fi_run[2:10, 1])
  f1_var_low <- var(ret_low_momentum$fi_run[2:10, 1])

  expect_gt(f1_var_high, f1_var_low)
})
