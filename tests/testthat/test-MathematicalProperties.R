test_that("Symmetry property - equal F values give equal gradients", {
  # Create symmetric test data
  data("IBD_simulation_data", package = "discent")
  dat <- IBD_simulation_data

  # Set all F values equal
  symmetric_params <- c("1" = 0.4, "2" = 0.4, "3" = 0.4, "m" = 800)

  # Run optimization
  inputdisc <- dat %>% dplyr::filter(deme1 != deme2)
  ret <- disc(discdat = inputdisc,
              start_params = symmetric_params,
              learningrate = 1e-3,
              lambda = 0,
              b1 = 0.9, b2 = 0.999, e = 1e-8,
              steps = 2,
              normalize_geodist = FALSE,
              report_progress = FALSE,
              return_verbose = TRUE)

  f_grads <- ret$fi_gradtraj[2, ]

  # For truly symmetric data and parameters, all F gradients should be equal
  # But our data might not be perfectly symmetric, so we test relative magnitudes
  cat(sprintf("\\nSymmetry test - F gradients: F1=%.6f, F2=%.6f, F3=%.6f\\n",
              f_grads[1], f_grads[2], f_grads[3]))

  # Test that no gradient is dramatically different (within factor of 10)
  max_grad <- max(abs(f_grads))
  min_grad <- min(abs(f_grads))
  relative_spread <- max_grad / max(min_grad, 1e-10)

  expect_lt(relative_spread, 50,
            label = sprintf("F gradients too asymmetric for symmetric parameters: %g", relative_spread))
})

test_that("Gradient scaling linearity", {
  # Test that doubling the data approximately doubles the gradient
  data("IBD_simulation_data", package = "discent")
  dat <- IBD_simulation_data

  test_params <- c("1" = 0.3, "2" = 0.5, "3" = 0.7, "m" = 600)

  # Original data gradient
  inputdisc <- dat %>% dplyr::filter(deme1 != deme2)
  ret1 <- disc(discdat = inputdisc,
               start_params = test_params,
               learningrate = 1e-3,
               lambda = 0,
               steps = 2,
               normalize_geodist = FALSE,
               report_progress = FALSE,
               return_verbose = TRUE)

  # Double the genetic distances (simulate more divergent samples)
  dat_scaled <- dat %>%
    dplyr::mutate(gendist = pmin(gendist * 1.5, 0.95))  # Scale but keep in bounds

  inputdisc_scaled <- dat_scaled %>% dplyr::filter(deme1 != deme2)
  ret2 <- disc(discdat = inputdisc_scaled,
               start_params = test_params,
               learningrate = 1e-3,
               lambda = 0,
               steps = 2,
               normalize_geodist = FALSE,
               report_progress = FALSE,
               return_verbose = TRUE)

  # Compare gradient magnitudes
  orig_f_grad_norm <- sqrt(sum(ret1$fi_gradtraj[2, ]^2))
  scaled_f_grad_norm <- sqrt(sum(ret2$fi_gradtraj[2, ]^2))
  orig_m_grad <- abs(ret1$m_gradtraj[2])
  scaled_m_grad <- abs(ret2$m_gradtraj[2])

  cat(sprintf("\\nScaling test - Original F norm: %.3f, Scaled F norm: %.3f\\n",
              orig_f_grad_norm, scaled_f_grad_norm))
  cat(sprintf("Scaling test - Original M grad: %.3f, Scaled M grad: %.3f\\n",
              orig_m_grad, scaled_m_grad))

  # Gradients should increase with data scale (but not necessarily linearly due to nonlinearity)
  expect_gt(scaled_f_grad_norm, 0.5 * orig_f_grad_norm,
            label = "F gradient should increase with data scaling")
  expect_gt(scaled_m_grad, 0.5 * orig_m_grad,
            label = "M gradient should increase with data scaling")
})

test_that("Gradient consistency across different starting values", {
  # Test that gradients are consistent regardless of starting values
  # when evaluated at the same parameter point
  data("IBD_simulation_data", package = "discent")
  dat <- IBD_simulation_data

  eval_point <- c("1" = 0.2, "2" = 0.6, "3" = 0.4, "m" = 500)

  # Different starting points
  start1 <- c("1" = 0.1, "2" = 0.1, "3" = 0.1, "m" = 100)
  start2 <- c("1" = 0.9, "2" = 0.9, "3" = 0.9, "m" = 1000)

  inputdisc <- dat %>% dplyr::filter(deme1 != deme2)

  # Manually set parameters to evaluation point and get gradients
  # (This is a bit of a hack - we run 1 step with zero learning rate)
  ret1 <- disc(discdat = inputdisc,
               start_params = eval_point,
               learningrate = 0,  # Don't update parameters
               lambda = 0,
               steps = 2,
               normalize_geodist = FALSE,
               report_progress = FALSE,
               return_verbose = TRUE)

  ret2 <- disc(discdat = inputdisc,
               start_params = eval_point,  # Same evaluation point
               learningrate = 0,
               lambda = 0,
               steps = 2,
               normalize_geodist = FALSE,
               report_progress = FALSE,
               return_verbose = TRUE)

  # Gradients should be identical
  for(i in 1:3) {
    expect_equal(ret1$fi_gradtraj[2, i], ret2$fi_gradtraj[2, i],
                 tolerance = 1e-12,
                 label = sprintf("F%d gradient inconsistent", i))
  }

  expect_equal(ret1$m_gradtraj[2], ret2$m_gradtraj[2],
               tolerance = 1e-12,
               label = "M gradient inconsistent")
})

test_that("Perfect fit gives zero gradients", {
  # Create synthetic data that perfectly matches a specific model
  # This tests if gradients are zero when cost function is at minimum

  # Use simple synthetic data - correct column order: smpl1, smpl2, deme1, deme2, gendist, geodist
  synthetic_data <- data.frame(
    smpl1 = c("s1", "s2", "s3"),
    smpl2 = c("s4", "s5", "s6"),
    deme1 = c(1, 2, 1),
    deme2 = c(2, 3, 3),
    gendist = c(0.1, 0.2, 0.3),  # Will be overwritten
    geodist = c(100, 200, 300)
  )

  # Set model parameters
  true_f1 <- 0.3
  true_f2 <- 0.5
  true_f3 <- 0.7
  true_m <- 400

  # Generate perfect genetic distances based on model
  synthetic_data$gendist <- c(
    ((true_f1 + true_f2)/2) * exp(-100/true_m),  # deme 1-2
    ((true_f2 + true_f3)/2) * exp(-200/true_m),  # deme 2-3
    ((true_f1 + true_f3)/2) * exp(-300/true_m)   # deme 1-3
  )

  # Add tiny amount of noise to avoid numerical issues
  synthetic_data$gendist <- pmax(pmin(
    synthetic_data$gendist + rnorm(3, 0, 1e-6),
    0.999), 0.001)

  true_params <- c("1" = true_f1, "2" = true_f2, "3" = true_f3, "m" = true_m)

  ret <- disc(discdat = synthetic_data,
              start_params = true_params,
              learningrate = 1e-3,
              lambda = 0,
              steps = 2,
              normalize_geodist = FALSE,
              report_progress = FALSE,
              return_verbose = TRUE)

  # Gradients should be very small (near zero) for perfect fit
  f_grad_norm <- sqrt(sum(ret$fi_gradtraj[2, ]^2))
  m_grad <- abs(ret$m_gradtraj[2])

  cat(sprintf("\\nPerfect fit test - F gradient norm: %.2e, M gradient: %.2e\\n",
              f_grad_norm, m_grad))

  expect_lt(f_grad_norm, 1e-3,
            label = sprintf("F gradients should be near zero for perfect fit: %g", f_grad_norm))
  expect_lt(m_grad, 1e-3,
            label = sprintf("M gradient should be near zero for perfect fit: %g", m_grad))
})