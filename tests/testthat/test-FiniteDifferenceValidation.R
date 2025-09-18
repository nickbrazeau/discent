test_that("Finite difference gradient validation - F gradients", {
  # Load data
  data("IBD_simulation_data", package = "discent")
  dat <- IBD_simulation_data

  # Test parameters
  test_params <- c("1" = 0.3, "2" = 0.7, "3" = 0.1, "m" = 500)  # Asymmetric to catch bugs
  h <- 1e-6  # finite difference step size
  tolerance <- 1e-5  # numerical tolerance

  # Define cost function for finite differences
  cost_function <- function(params, data) {
    f_vals <- params[1:3]
    m_val <- params[4]

    # Filter data
    clean_data <- data %>%
      dplyr::filter(deme1 != deme2) %>%
      dplyr::mutate(gendist = pmax(pmin(gendist, 0.999), 0.001)) %>%
      dplyr::mutate(gendist = logit(gendist))

    total_cost <- 0
    lambda <- 0  # No regularization for this test

    # Calculate cost using same logic as C++ (upper triangle only)
    for(i in 1:(length(f_vals)-1)) {
      for(j in (i+1):length(f_vals)) {
        pairs_ij <- clean_data %>%
          dplyr::filter((deme1 == i & deme2 == j) | (deme1 == j & deme2 == i))

        if(nrow(pairs_ij) > 0) {
          f_avg <- (f_vals[i] + f_vals[j]) / 2
          r_ij <- f_avg * exp(-pairs_ij$geodist / m_val)

          cost_contributions <- (pairs_ij$gendist - r_ij)^2
          total_cost <- total_cost + sum(cost_contributions, na.rm = TRUE)
        }
      }
    }

    return(total_cost)
  }

  # Calculate analytical gradients using disc
  inputdisc <- dat %>% dplyr::filter(deme1 != deme2)
  ret <- disc(discdat = inputdisc,
              start_params = test_params,
              learningrate = 1e-3,
              lambda = 0,  # No regularization
              b1 = 0.9, b2 = 0.999, e = 1e-8,
              steps = 2,  # Just need first gradient
              normalize_geodist = FALSE,
              report_progress = FALSE,
              return_verbose = TRUE)

  analytical_f_grads <- ret$fi_gradtraj[2, ]  # Second step gradients

  # Calculate numerical gradients using finite differences
  base_params <- c(logit(test_params[1:3]), test_params[4])
  base_cost <- cost_function(base_params, dat)

  numerical_f_grads <- numeric(3)

  for(i in 1:3) {
    # Forward difference
    params_plus <- base_params
    params_plus[i] <- params_plus[i] + h
    cost_plus <- cost_function(params_plus, dat)

    # Backward difference
    params_minus <- base_params
    params_minus[i] <- params_minus[i] - h
    cost_minus <- cost_function(params_minus, dat)

    # Central difference
    numerical_f_grads[i] <- (cost_plus - cost_minus) / (2 * h)
  }

  # Compare analytical vs numerical gradients
  for(i in 1:3) {
    testthat::expect_equal(analytical_f_grads[i], numerical_f_grads[i],  tolerance = tolerance)

  }
})

test_that("Finite difference gradient validation - M gradient", {
  # Load data
  data("IBD_simulation_data", package = "discent")
  dat <- IBD_simulation_data

  # Test parameters
  test_params <- c("1" = 0.3, "2" = 0.7, "3" = 0.1, "m" = 500)
  h <- 1e-6
  tolerance <- 1e-5
  lambda_test <- 0.1  # Test with regularization

  # Define cost function with regularization
  cost_function_with_reg <- function(params, data, lambda) {
    f_vals <- params[1:3]
    m_val <- params[4]

    # Filter data
    clean_data <- data %>%
      dplyr::filter(deme1 != deme2) %>%
      dplyr::mutate(gendist = pmax(pmin(gendist, 0.999), 0.001)) %>%
      dplyr::mutate(gendist = logit(gendist))

    total_cost <- 0

    # Main cost term (without regularization per sample)
    for(i in 1:(length(f_vals)-1)) {
      for(j in (i+1):length(f_vals)) {
        pairs_ij <- clean_data %>%
          dplyr::filter((deme1 == i & deme2 == j) | (deme1 == j & deme2 == i))

        if(nrow(pairs_ij) > 0) {
          f_avg <- (f_vals[i] + f_vals[j]) / 2
          r_ij <- f_avg * exp(-pairs_ij$geodist / m_val)

          cost_contributions <- (pairs_ij$gendist - r_ij)^2  # No regularization here
          total_cost <- total_cost + sum(cost_contributions, na.rm = TRUE)
        }
      }
    }

    # Add regularization term once total (matching C++ implementation)
    total_cost <- total_cost + lambda * m_val^2

    return(total_cost)
  }

  # Calculate analytical M gradient
  inputdisc <- dat %>% dplyr::filter(deme1 != deme2)
  ret <- disc(discdat = inputdisc,
              start_params = test_params,
              learningrate = 1e-3,
              lambda = lambda_test,
              b1 = 0.9, b2 = 0.999, e = 1e-8,
              steps = 2,
              normalize_geodist = FALSE,
              report_progress = FALSE,
              return_verbose = TRUE)

  analytical_m_grad <- as.numeric(ret$m_gradtraj[2])

  # Calculate numerical M gradient
  base_params <- c(logit(test_params[1:3]), test_params[4])

  # Forward difference
  params_plus <- base_params
  params_plus[4] <- params_plus[4] + h
  cost_plus <- cost_function_with_reg(params_plus, dat, lambda_test)

  # Backward difference
  params_minus <- base_params
  params_minus[4] <- params_minus[4] - h
  cost_minus <- cost_function_with_reg(params_minus, dat, lambda_test)

  # Central difference
  numerical_m_grad <- (cost_plus - cost_minus) / (2 * h)

  # Compare
  cat(sprintf("\\nM Gradient - Analytical: %g, Numerical: %g, Diff: %g\\n",
              analytical_m_grad, numerical_m_grad,
              abs(analytical_m_grad - numerical_m_grad)))

  expect_equal(as.numeric(analytical_m_grad), as.numeric(numerical_m_grad), tolerance = tolerance,
               label = "M gradient mismatch")
})
