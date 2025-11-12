## test inspired and partly written by Claude Code

# load common function
cost_function <- function(params, data) {
  # init
  f_vals <- params[c("1", "2", "3")]
  m_val <-  params[c("m")]
  total_cost <- 0

  # Calculate cost using same logic as C++ (upper triangle only)
  for(i in 1:(length(f_vals)-1)) {
    for(j in (i+1):length(f_vals)) {
      pairs_ij <- data %>%
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

test_that("Finite difference gradient validation - F gradients", {
  # Load data
  data("IBD_simulation_data", package = "discent")
  dat <- IBD_simulation_data %>%
    dplyr::filter(deme1 != deme2) %>%
    dplyr::mutate(gendist = ifelse(gendist > 0.9999, 0.9999,
                                   ifelse(gendist < 0.0001, 0.0001,
                                          gendist)))

  test_params <- c("1" = 0.3, "2" = 0.7, "3" = 0.1, "m" = 50)
  h <- 1e-6  # finite difference step size
  tolerance <- 1e-5  # numerical tolerance


  # Calculate gradients using disc
  ret <- disc(discdat = dat,
              start_params = test_params,
              learningrate = 1e-3,
              lambda = 0,
              b1 = 0.9, b2 = 0.999, e = 1e-8,
              steps = 2,
              normalize_geodist = FALSE,
              report_progress = FALSE,
              return_verbose = TRUE)

  # init
  numerical_f_grads <- rep(0, 3)
  for(i in 1:3) {
    # Forward difference
    params_plus <- test_params
    params_plus[i] <- params_plus[i] + h
    cost_plus <- cost_function(params_plus, data = dat)

    # Backward difference
    params_minus <- test_params
    params_minus[i] <- params_minus[i] - h
    cost_minus <- cost_function(params_minus, data  = dat)

    # Central difference
    numerical_f_grads[i] <- (cost_plus - cost_minus) / (2 * h)
  }

  # Compare disc vs finite diff gradients
  disc_grads <- ret$fi_gradtraj[2,]/(ret$fi_run[2,] * (1-ret$fi_run[2,]))

  for(i in 1:3) {
    testthat::expect_equal(disc_grads[i],
                           numerical_f_grads[i],
                           tolerance = tolerance)

  }
})

test_that("Finite difference gradient validation - M gradient", {
  # Load data
  data("IBD_simulation_data", package = "discent")
  dat <- IBD_simulation_data %>%
    dplyr::filter(deme1 != deme2) %>%
    dplyr::mutate(gendist = ifelse(gendist > 0.9999, 0.9999,
                                   ifelse(gendist < 0.0001, 0.0001,
                                          gendist)))
  # Test parameters
  test_params <- c("1" = 0.3, "2" = 0.7, "3" = 0.1, "m" = 500)
  h <- 1e-6
  tolerance <- 1e-5

  # Calculate analytical M gradient
  ret <- disc(discdat = dat,
              start_params = test_params,
              learningrate = 1e-3,
              lambda = 0,
              b1 = 0.9, b2 = 0.999, e = 1e-8,
              steps = 2,
              normalize_geodist = FALSE,
              report_progress = FALSE,
              return_verbose = TRUE)


  # Forward difference
  params_plus <- test_params
  params_plus[4] <- test_params[4] + h
  cost_plus <- cost_function(params_plus, data = dat)

  # Backward difference
  params_minus <- test_params
  params_minus[4] <- params_minus[4] - h
  cost_minus <- cost_function(params_minus, data = dat)

  # Compare disc vs finite diff gradients
  numerical_m_grad <- (cost_plus - cost_minus) / (2 * h)
  disc_grad <- ret$m_gradtraj[2]/ret$m_run[2]


  # Compare
  expect_equal(numerical_m_grad, disc_grad, tolerance = tolerance)

})
