test_that("Extreme F values near boundaries", {
  # Test behavior when F values are near 0 or 1 (boundaries of logit space)
  data("IBD_simulation_data", package = "discent")
  dat <- IBD_simulation_data

  # Test near lower boundary
  extreme_low_params <- c("1" = 0.005, "2" = 0.01, "3" = 0.001, "m" = 500)

  inputdisc <- dat %>% dplyr::filter(deme1 != deme2)
  ret_low <- disc(discdat = inputdisc,
                  start_params = extreme_low_params,
                  learningrate = 1e-4,  # Smaller learning rate for stability
                  lambda = 0,
                  steps = 2,
                  normalize_geodist = FALSE,
                  report_progress = FALSE,
                  return_verbose = TRUE)

  # Test near upper boundary
  extreme_high_params <- c("1" = 0.995, "2" = 0.99, "3" = 0.999, "m" = 500)

  ret_high <- disc(discdat = inputdisc,
                   start_params = extreme_high_params,
                   learningrate = 1e-4,
                   lambda = 0,
                   steps = 2,
                   normalize_geodist = FALSE,
                   report_progress = FALSE,
                   return_verbose = TRUE)

  # Check for numerical stability - no NaN/Inf values
  expect_true(all(is.finite(ret_low$fi_gradtraj[2, ])),
              label = "Low boundary F gradients should be finite")
  expect_true(is.finite(ret_low$m_gradtraj[2]),
              label = "Low boundary M gradient should be finite")
  expect_true(all(is.finite(ret_high$fi_gradtraj[2, ])),
              label = "High boundary F gradients should be finite")
  expect_true(is.finite(ret_high$m_gradtraj[2]),
              label = "High boundary M gradient should be finite")

  cat(sprintf("\\nBoundary test - Low F gradients: [%.3f, %.3f, %.3f]\\n",
              ret_low$fi_gradtraj[2,1], ret_low$fi_gradtraj[2,2], ret_low$fi_gradtraj[2,3]))
  cat(sprintf("Boundary test - High F gradients: [%.3f, %.3f, %.3f]\\n",
              ret_high$fi_gradtraj[2,1], ret_high$fi_gradtraj[2,2], ret_high$fi_gradtraj[2,3]))
})

test_that("Extreme M values (migration rate bounds)", {
  # Test behavior with very small and very large migration rates
  data("IBD_simulation_data", package = "discent")
  dat <- IBD_simulation_data

  test_f_params <- c("1" = 0.3, "2" = 0.5, "3" = 0.7)

  inputdisc <- dat %>% dplyr::filter(deme1 != deme2)

  # Very small M (high isolation)
  small_m_params <- c(test_f_params, "m" = 10)
  ret_small_m <- disc(discdat = inputdisc,
                      start_params = small_m_params,
                      learningrate = 1e-3,
                      lambda = 0,
                      steps = 2,
                      normalize_geodist = FALSE,
                      report_progress = FALSE,
                      return_verbose = TRUE)

  # Very large M (high migration)
  large_m_params <- c(test_f_params, "m" = 10000)
  ret_large_m <- disc(discdat = inputdisc,
                      start_params = large_m_params,
                      learningrate = 1e-3,
                      lambda = 0,
                      steps = 2,
                      normalize_geodist = FALSE,
                      report_progress = FALSE,
                      return_verbose = TRUE)

  # Check numerical stability
  expect_true(all(is.finite(ret_small_m$fi_gradtraj[2, ])),
              label = "Small M - F gradients should be finite")
  expect_true(is.finite(ret_small_m$m_gradtraj[2]),
              label = "Small M - M gradient should be finite")
  expect_true(all(is.finite(ret_large_m$fi_gradtraj[2, ])),
              label = "Large M - F gradients should be finite")
  expect_true(is.finite(ret_large_m$m_gradtraj[2]),
              label = "Large M - M gradient should be finite")

  cat(sprintf("\\nM bounds test - Small M (10): M grad = %.6f\\n", ret_small_m$m_gradtraj[2]))
  cat(sprintf("M bounds test - Large M (10000): M grad = %.6f\\n", ret_large_m$m_gradtraj[2]))
})

test_that("Sparse deme pair data", {
  # Test with limited data between some deme pairs (but keeping all demes)
  data("IBD_simulation_data", package = "discent")
  dat <- IBD_simulation_data

  # Keep all demes but filter to reduce data volume
  sparse_data <- dat %>%
    dplyr::filter(deme1 != deme2) %>%
    dplyr::group_by(deme1, deme2) %>%
    dplyr::slice(1:2) %>%  # Keep only first 2 samples per deme pair
    dplyr::ungroup()

  if(nrow(sparse_data) > 0) {
    test_params <- c("1" = 0.3, "2" = 0.5, "3" = 0.7, "m" = 400)

    ret_sparse <- disc(discdat = sparse_data,
                       start_params = test_params,
                       learningrate = 1e-3,
                       lambda = 0,
                       steps = 2,
                       normalize_geodist = FALSE,
                       report_progress = FALSE,
                       return_verbose = TRUE)

    # All gradients should be finite
    expect_true(all(is.finite(ret_sparse$fi_gradtraj[2, ])),
                label = "All F gradients should be finite with sparse data")
    expect_true(is.finite(ret_sparse$m_gradtraj[2]),
                label = "M gradient should be finite with sparse data")

    cat(sprintf("\\nSparse data test - F gradients: [%.6f, %.6f, %.6f]\\n",
                ret_sparse$fi_gradtraj[2,1], ret_sparse$fi_gradtraj[2,2], ret_sparse$fi_gradtraj[2,3]))
  } else {
    skip("No valid sparse data created")
  }
})

test_that("High regularization lambda values", {
  # Test behavior with strong regularization
  data("IBD_simulation_data", package = "discent")
  dat <- IBD_simulation_data

  test_params <- c("1" = 0.3, "2" = 0.5, "3" = 0.7, "m" = 400)

  inputdisc <- dat %>% dplyr::filter(deme1 != deme2)

  # Test with very high lambda
  ret_high_lambda <- disc(discdat = inputdisc,
                         start_params = test_params,
                         learningrate = 1e-3,
                         lambda = 100,  # Very high regularization
                         steps = 2,
                         normalize_geodist = FALSE,
                         report_progress = FALSE,
                         return_verbose = TRUE)

  # M gradient should be dominated by regularization term: 2*lambda*m
  expected_reg_component <- 2 * 100 * 400  # 2 * lambda * m
  actual_m_grad <- ret_high_lambda$m_gradtraj[2]

  cat(sprintf("\\nHigh lambda test - M gradient: %.1f, Expected reg component: %.1f\\n",
              actual_m_grad, expected_reg_component))

  # The regularization component should dominate
  expect_gt(abs(actual_m_grad), expected_reg_component * 0.5,
            label = "High lambda should dominate M gradient")

  # Check numerical stability
  expect_true(is.finite(actual_m_grad),
              label = "M gradient should be finite with high lambda")
  expect_true(all(is.finite(ret_high_lambda$fi_gradtraj[2, ])),
              label = "F gradients should be finite with high lambda")
})

test_that("Low genetic distances (nearly identical samples)", {
  # Test behavior when genetic distance is very low (nearly identical samples)
  data("IBD_simulation_data", package = "discent")
  dat <- IBD_simulation_data

  # Create data with very low genetic distances but preserve geodistance structure
  # Ensure we have representatives from all deme pairs
  low_dist_data <- dat %>%
    dplyr::filter(deme1 != deme2) %>%
    dplyr::group_by(deme1, deme2) %>%
    dplyr::slice(1:3) %>%  # Take first 3 from each deme pair
    dplyr::ungroup() %>%
    dplyr::mutate(gendist = pmax(gendist * 0.001, 0.001))  # Scale down genetic distances

  test_params <- c("1" = 0.3, "2" = 0.5, "3" = 0.7, "m" = 400)

  ret_low_gen <- disc(discdat = low_dist_data,
                      start_params = test_params,
                      learningrate = 1e-3,
                      lambda = 0,
                      steps = 2,
                      normalize_geodist = FALSE,
                      report_progress = FALSE,
                      return_verbose = TRUE)

  # Check for numerical stability with very low genetic distances
  expect_true(all(is.finite(ret_low_gen$fi_gradtraj[2, ])),
              label = "F gradients should be finite with low genetic distances")
  expect_true(is.finite(ret_low_gen$m_gradtraj[2]),
              label = "M gradient should be finite with low genetic distances")

  # Cost should not be NaN or infinite
  expect_true(all(is.finite(ret_low_gen$cost)),
              label = "Cost should be finite with low genetic distances")

  cat(sprintf("\\nLow genetic distance test - Cost: %.3f, M grad: %.3f\\n",
              ret_low_gen$cost[2], ret_low_gen$m_gradtraj[2]))
})