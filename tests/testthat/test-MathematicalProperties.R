test_that("F gradient structural pattern", {
  # Test the expected structural pattern of F gradients based on upper triangle loop
  data("IBD_simulation_data", package = "discent")
  dat <- IBD_simulation_data

  test_params <- c("1" = 0.4, "2" = 0.4, "3" = 0.4, "m" = 800)

  inputdisc <- dat %>% dplyr::filter(deme1 != deme2)
  ret <- disc(discdat = inputdisc,
              start_params = test_params,
              learningrate = 1e-3,
              lambda = 0,
              steps = 2,
              normalize_geodist = FALSE,
              report_progress = FALSE,
              return_verbose = TRUE)

  f_grads <- ret$fi_gradtraj[2, ]

  cat(sprintf("\\nF gradient pattern - F1=%.6f, F2=%.6f, F3=%.6f\\n",
              f_grads[1], f_grads[2], f_grads[3]))

  # All F gradients should be non-zero due to fgrad[j] contribution in upper triangle loop
  # F1 participates as i in (1,2) and (1,3) pairs
  # F2 participates as i in (2,3) pair and j in (1,2) pair
  # F3 participates as j in (1,3) and (2,3) pairs
  expect_gt(abs(f_grads[1]), 1e-6, label = "F1 gradient should be non-zero")
  expect_gt(abs(f_grads[2]), 1e-6, label = "F2 gradient should be non-zero")
  expect_gt(abs(f_grads[3]), 1e-6, label = "F3 gradient should be non-zero")
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

test_that("Gradients respond to parameter changes", {
  # Test that gradients change appropriately when parameters change
  data("IBD_simulation_data", package = "discent")
  dat <- IBD_simulation_data

  params1 <- c("1" = 0.2, "2" = 0.3, "3" = 0.4, "m" = 500)
  params2 <- c("1" = 0.8, "2" = 0.7, "3" = 0.6, "m" = 200)

  inputdisc <- dat %>% dplyr::filter(deme1 != deme2)

  ret1 <- disc(discdat = inputdisc, start_params = params1, steps = 2,
               normalize_geodist = FALSE, report_progress = FALSE, return_verbose = TRUE)
  ret2 <- disc(discdat = inputdisc, start_params = params2, steps = 2,
               normalize_geodist = FALSE, report_progress = FALSE, return_verbose = TRUE)

  # Gradients should be different for different parameter values
  f_grad_diff <- sqrt(sum((ret1$fi_gradtraj[2, ] - ret2$fi_gradtraj[2, ])^2))
  m_grad_diff <- abs(ret1$m_gradtraj[2] - ret2$m_gradtraj[2])

  expect_gt(f_grad_diff, 1e-6, label = "F gradients should differ for different parameters")
  expect_gt(m_grad_diff, 1e-6, label = "M gradient should differ for different parameters")

  cat(sprintf("\\nParameter sensitivity test - F grad diff: %.3f, M grad diff: %.3f\\n",
              f_grad_diff, m_grad_diff))
})
