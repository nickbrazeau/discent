# test inspired and in part written by Claude Code (11/2025)

test_that("Extreme F values near boundaries", {
  data("IBD_simulation_data", package = "discent")
  dat <- IBD_simulation_data

  # Test near lower boundary
  extreme_low_params <- c("1" = 0.005, "2" = 0.01, "3" = 0.001, "m" = 50)

  inputdisc <- dat %>% dplyr::filter(deme1 != deme2)
  ret_low <- disc(discdat = inputdisc,
                  start_params = extreme_low_params,
                  learningrate = 1e-4,
                  lambda = 0,
                  steps = 100,
                  normalize_geodist = FALSE,
                  report_progress = FALSE,
                  return_verbose = TRUE)

  # Test near upper boundary
  extreme_high_params <- c("1" = 0.995, "2" = 0.99, "3" = 0.999, "m" = 500)

  ret_high <- disc(discdat = inputdisc,
                   start_params = extreme_high_params,
                   learningrate = 1e-4,
                   lambda = 0,
                   steps = 100,
                   normalize_geodist = FALSE,
                   report_progress = FALSE,
                   return_verbose = TRUE)

  # Check for numerical stability - no NaN/Inf values
  expect_true(all(is.finite(ret_low$fi_gradtraj)),
              label = "Low boundary F gradients should be finite")
  expect_true(all(is.finite(ret_low$m_gradtraj)),
              label = "Low boundary M gradient should be finite")
  expect_true(all(is.finite(ret_high$fi_gradtraj)),
              label = "High boundary F gradients should be finite")
  expect_true(all(is.finite(ret_high$m_gradtraj)),
              label = "High boundary M gradient should be finite")

})

test_that("Extreme M values", {
  # Test behavior with very small and very large migration rates
  data("IBD_simulation_data", package = "discent")
  dat <- IBD_simulation_data

  test_f_params <- c("1" = 0.3, "2" = 0.5, "3" = 0.7)

  inputdisc <- dat %>% dplyr::filter(deme1 != deme2)

  # Very small M (high isolation)
  small_m_params <- c(test_f_params, "m" = 0.01)
  ret_small_m <- disc(discdat = inputdisc,
                      start_params = small_m_params,
                      learningrate = 1e-3,
                      lambda = 0,
                      steps = 100,
                      normalize_geodist = FALSE,
                      report_progress = FALSE,
                      return_verbose = TRUE)

  # Very large M (high migration)
  large_m_params <- c(test_f_params, "m" = 10000)
  ret_large_m <- disc(discdat = inputdisc,
                      start_params = large_m_params,
                      learningrate = 1e-3,
                      lambda = 0,
                      steps = 100,
                      normalize_geodist = FALSE,
                      report_progress = FALSE,
                      return_verbose = TRUE)

  # Check numerical stability
  expect_true(all(is.finite(ret_small_m$fi_gradtraj)),
              label = "Low boundary F gradients should be finite")
  expect_true(all(is.finite(ret_small_m$m_gradtraj)),
              label = "Low boundary M gradient should be finite")
  expect_true(all(is.finite(ret_large_m$fi_gradtraj)),
              label = "High boundary F gradients should be finite")
  expect_true(all(is.finite(ret_large_m$m_gradtraj)),
              label = "High boundary M gradient should be finite")

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

    test_params <- c("1" = 0.3, "2" = 0.5, "3" = 0.7, "m" = 400)

    ret_sparse <- disc(discdat = sparse_data,
                       start_params = test_params,
                       learningrate = 1e-3,
                       lambda = 0,
                       steps = 2,
                       normalize_geodist = FALSE,
                       report_progress = FALSE,
                       return_verbose = TRUE)

    # Check numerical stability
    expect_true(all(is.finite(ret_sparse$fi_gradtraj)),
                label = "F gradients should be finite")
    expect_true(all(is.finite(ret_sparse$m_gradtraj)),
                label = "M gradient should be finite")

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
                         lambda = 1e3,
                         steps = 2,
                         normalize_geodist = FALSE,
                         report_progress = FALSE,
                         return_verbose = TRUE)


  # Check numerical stability
  expect_true(all(is.finite(ret_high_lambda$fi_gradtraj)),
              label = "F gradients should be finite")
  expect_true(all(is.finite(ret_high_lambda$m_gradtraj)),
              label = "M gradient should be finite")
})

test_that("Low genetic distances (nearly identical samples)", {
  # Test behavior when genetic distance is very low (nearly identical samples)
  data("IBD_simulation_data", package = "discent")
  dat <- IBD_simulation_data

  # Create data with very low genetic distances but preserve geodistance structure
  # Ensure we have representatives from all deme pairs
  low_dist_data <- dat %>%
    dplyr::filter(deme1 != deme2) %>%
    dplyr::mutate(gendist = max(gendist * 0.001, 0.001))

  test_params <- c("1" = 0.3, "2" = 0.5, "3" = 0.7, "m" = 40)

  ret_low_gen <- disc(discdat = low_dist_data,
                      start_params = test_params,
                      learningrate = 1e-3,
                      lambda = 0,
                      steps = 2,
                      normalize_geodist = FALSE,
                      report_progress = FALSE,
                      return_verbose = TRUE)

  # Check numerical stability
  expect_true(all(is.finite(ret_low_gen$fi_gradtraj)),
              label = "F gradients should be finite")
  expect_true(all(is.finite(ret_low_gen$m_gradtraj)),
              label = "M gradient should be finite")

  # Cost should not be NaN or infinite
  expect_true(all(is.finite(ret_low_gen$cost)),
              label = "Cost should be finite with low genetic distances")

})
