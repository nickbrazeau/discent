test_that("Regularization term isolation", {
  # Test that regularization term 2*lambda*m is correctly added to M gradient
  data("IBD_simulation_data", package = "discent")
  dat <- IBD_simulation_data

  test_params <- c("1" = 0.3, "2" = 0.5, "3" = 0.7, "m" = 400)
  lambda_test <- 0.5

  inputdisc <- dat %>% dplyr::filter(deme1 != deme2)

  # Run with regularization
  ret_with_reg <- disc(discdat = inputdisc,
                       start_params = test_params,
                       learningrate = 1e-3,
                       lambda = lambda_test,
                       steps = 2,
                       normalize_geodist = FALSE,
                       report_progress = FALSE,
                       return_verbose = TRUE)

  # Run without regularization
  ret_no_reg <- disc(discdat = inputdisc,
                     start_params = test_params,
                     learningrate = 1e-3,
                     lambda = 0,
                     steps = 2,
                     normalize_geodist = FALSE,
                     report_progress = FALSE,
                     return_verbose = TRUE)

  # The difference in M gradients should be exactly 2*lambda*m
  m_grad_with_reg <- ret_with_reg$m_gradtraj[2]
  m_grad_no_reg <- ret_no_reg$m_gradtraj[2]
  expected_reg_difference <- 2 * lambda_test * unname(test_params["m"])

  actual_difference <- m_grad_with_reg - m_grad_no_reg

  cat(sprintf("\\nRegularization test - M grad with reg: %.6f, without reg: %.6f\\n",
              m_grad_with_reg, m_grad_no_reg))
  cat(sprintf("Actual difference: %.6f, Expected (2λm): %.6f\\n",
              actual_difference, expected_reg_difference))

  expect_equal(actual_difference, expected_reg_difference, tolerance = 1e-6,
               label = "Regularization difference should equal 2*lambda*m")

  # F gradients should be identical (regularization only affects M)
  for(i in 1:3) {
    expect_equal(ret_with_reg$fi_gradtraj[2, i], ret_no_reg$fi_gradtraj[2, i],
                 tolerance = 1e-8,
                 label = sprintf("F%d gradient should be unaffected by lambda", i))
  }
})

test_that("Regularization scaling with lambda", {
  # Test that regularization effect scales linearly with lambda
  data("IBD_simulation_data", package = "discent")
  dat <- IBD_simulation_data

  test_params <- c("1" = 0.4, "2" = 0.6, "3" = 0.2, "m" = 300)
  lambdas <- c(0.1, 0.2, 0.5, 1.0)

  inputdisc <- dat %>% dplyr::filter(deme1 != deme2)

  m_grads <- numeric(length(lambdas))

  for(i in seq_along(lambdas)) {
    ret <- disc(discdat = inputdisc,
                start_params = test_params,
                learningrate = 1e-3,
                lambda = lambdas[i],
                steps = 2,
                normalize_geodist = FALSE,
                report_progress = FALSE,
                return_verbose = TRUE)
    m_grads[i] <- ret$m_gradtraj[2]
  }

  # Calculate differences from lambda=0 case
  ret_no_reg <- disc(discdat = inputdisc,
                     start_params = test_params,
                     learningrate = 1e-3,
                     lambda = 0,
                     steps = 2,
                     normalize_geodist = FALSE,
                     report_progress = FALSE,
                     return_verbose = TRUE)
  baseline_m_grad <- ret_no_reg$m_gradtraj[2]

  reg_contributions <- m_grads - baseline_m_grad
  expected_contributions <- 2 * lambdas * test_params["m"]

  cat("\\nLambda scaling test:\\n")
  for(i in seq_along(lambdas)) {
    cat(sprintf("λ=%.1f: Actual reg contribution=%.1f, Expected=%.1f\\n",
                lambdas[i], reg_contributions[i], expected_contributions[i]))

    expect_equal(reg_contributions[i], expected_contributions[i], tolerance = 1e-10,
                 label = sprintf("Lambda %.1f scaling", lambdas[i]))
  }
})

test_that("Regularization with different M values", {
  # Test that regularization scales with M value (2*lambda*M)
  data("IBD_simulation_data", package = "discent")
  dat <- IBD_simulation_data

  lambda_test <- 0.3
  m_values <- c(100, 500, 1000)
  f_params <- c("1" = 0.4, "2" = 0.5, "3" = 0.6)

  inputdisc <- dat %>% dplyr::filter(deme1 != deme2)

  reg_contributions <- numeric(length(m_values))

  for(i in seq_along(m_values)) {
    test_params <- c(f_params, "m" = m_values[i])

    # With regularization
    ret_with_reg <- disc(discdat = inputdisc,
                         start_params = test_params,
                         learningrate = 1e-3,
                         lambda = lambda_test,
                         steps = 2,
                         normalize_geodist = FALSE,
                         report_progress = FALSE,
                         return_verbose = TRUE)

    # Without regularization
    ret_no_reg <- disc(discdat = inputdisc,
                       start_params = test_params,
                       learningrate = 1e-3,
                       lambda = 0,
                       steps = 2,
                       normalize_geodist = FALSE,
                       report_progress = FALSE,
                       return_verbose = TRUE)

    reg_contributions[i] <- ret_with_reg$m_gradtraj[2] - ret_no_reg$m_gradtraj[2]
  }

  expected_contributions <- 2 * lambda_test * m_values

  cat("\\nM value scaling test:\\n")
  for(i in seq_along(m_values)) {
    cat(sprintf("M=%d: Actual reg contribution=%.1f, Expected=%.1f\\n",
                m_values[i], reg_contributions[i], expected_contributions[i]))

    expect_equal(reg_contributions[i], expected_contributions[i], tolerance = 1e-10,
                 label = sprintf("M value %d scaling", m_values[i]))
  }
})

test_that("Regularization cost function contribution", {
  # Test that regularization contributes lambda*M^2 to the cost function
  data("IBD_simulation_data", package = "discent")
  dat <- IBD_simulation_data

  test_params <- c("1" = 0.3, "2" = 0.7, "3" = 0.4, "m" = 600)
  lambda_test <- 0.2

  inputdisc <- dat %>% dplyr::filter(deme1 != deme2)

  # With regularization
  ret_with_reg <- disc(discdat = inputdisc,
                       start_params = test_params,
                       learningrate = 0,  # No parameter updates
                       lambda = lambda_test,
                       steps = 2,
                       normalize_geodist = FALSE,
                       report_progress = FALSE,
                       return_verbose = TRUE)

  # Without regularization
  ret_no_reg <- disc(discdat = inputdisc,
                     start_params = test_params,
                     learningrate = 0,
                     lambda = 0,
                     steps = 2,
                     normalize_geodist = FALSE,
                     report_progress = FALSE,
                     return_verbose = TRUE)

  # Cost difference should be proportional to lambda*M^2
  cost_with_reg <- ret_with_reg$cost[2]
  cost_no_reg <- ret_no_reg$cost[2]
  cost_difference <- cost_with_reg - cost_no_reg

  # The regularization term should contribute exactly lambda*M^2 total
  expected_reg_total <- lambda_test * (unname(test_params["m"]))^2
  num_pairs <- sum(inputdisc$deme1 != inputdisc$deme2)

  cat(sprintf("\\nCost regularization test - Cost difference: %.1f\\n", cost_difference))
  cat(sprintf("Expected reg total: %.1f, Number of pairs: %d\\n",
              expected_reg_total, num_pairs))

  # The regularization should increase the cost
  expect_gt(cost_difference, 0,
            label = "Regularization should increase cost")

  # The increase should equal lambda*M^2 (within tolerance)
  expect_equal(cost_difference, expected_reg_total, tolerance = 1e-6,
               label = "Cost increase should equal lambda*M^2")
})

test_that("Zero regularization equivalence", {
  # Test that lambda=0 produces identical results to no regularization
  data("IBD_simulation_data", package = "discent")
  dat <- IBD_simulation_data

  test_params <- c("1" = 0.4, "2" = 0.3, "3" = 0.8, "m" = 700)

  inputdisc <- dat %>% dplyr::filter(deme1 != deme2)

  # With lambda=0
  ret_lambda_zero <- disc(discdat = inputdisc,
                          start_params = test_params,
                          learningrate = 1e-3,
                          lambda = 0,
                          steps = 5,
                          normalize_geodist = FALSE,
                          report_progress = FALSE,
                          return_verbose = TRUE)

  # Explicit lambda=0 (should be identical to above)
  ret_default <- disc(discdat = inputdisc,
                      start_params = test_params,
                      learningrate = 1e-3,
                      lambda = 0,
                      steps = 5,
                      normalize_geodist = FALSE,
                      report_progress = FALSE,
                      return_verbose = TRUE)

  # All results should be identical
  expect_equal(ret_lambda_zero$m_gradtraj, ret_default$m_gradtraj, tolerance = 1e-6,
               label = "M gradients should be identical with lambda=0")
  expect_equal(ret_lambda_zero$fi_gradtraj, ret_default$fi_gradtraj, tolerance = 1e-8,
               label = "F gradients should be identical with lambda=0")
  expect_equal(ret_lambda_zero$cost, ret_default$cost, tolerance = 1e-6,
               label = "Costs should be identical with lambda=0")
  expect_equal(ret_lambda_zero$m_run, ret_default$m_run, tolerance = 1e-6,
               label = "M values should be identical with lambda=0")
  expect_equal(ret_lambda_zero$fi_run, ret_default$fi_run, tolerance = 1e-12,
               label = "F values should be identical with lambda=0")

  cat("\\nZero regularization test - All outputs identical: PASS\\n")
})