test_that("Fi adam by hand", {
  # sim data
  data("IBD_simulation_data", package = "discent")
  dat <- IBD_simulation_data

  # start params
  our_start_params <- rep(0.2, 3)
  names(our_start_params) <- 1:3
  our_start_params <- c(our_start_params, "m" = 1e3)

  # tidy date and calculate gradient for F1
  input <- dat %>%
    dplyr::filter(deme1 != deme2) %>%
    dplyr::filter(deme1 == 1 | deme2 == 1)


  input <- input %>%
    dplyr::mutate(gendist = ifelse(gendist > 0.999, 0.999,
                                   ifelse(gendist < 0.001, 0.001,
                                          gendist)))

  # now run model
  inputdisc <- dat %>%
    dplyr::filter(deme1 != deme2)
  ret <- disc(discdat = inputdisc,
              start_params = our_start_params,
              learningrate = 1e-3,
              b1 = 0.9,
              b2 = 0.999,
              e = 1e-8,
              steps = 1e3,
              normalize_geodist = F,
              report_progress = T,
              return_verbose = T)
  # back out for Fi adam
  b1 <- 0.9
  b2 <- 0.999
  e <- 1e-8
  learningrate <- 1e-3
  # init case has no prior so that is why b1 * 0
  # perform updates in reparam space
  mt_f1 <- b1 * 0 + (1-b1) *  ret$fi_gradtraj[2,1]
  vt_f1 <- b2 * 0 + (1-b2) *  ret$fi_gradtraj[2,1] *  ret$fi_gradtraj[2,1]
  mt_f1hat <- mt_f1 / (1 - b1^1)
  vt_f1hat <- vt_f1 / (1 - b2^1)
  fnew1 = expit( logit(ret$fi_run[1,1]) - learningrate * (mt_f1hat/(sqrt(vt_f1hat) + e)) )

  # test out
  testthat::expect_equal(fnew1, ret$fi_run[2,1])

})


test_that("M adam by hand", {
  # sim data
  data("IBD_simulation_data", package = "discent")
  dat <- IBD_simulation_data

  # start params
  our_start_params <- rep(0.2, 3)
  names(our_start_params) <- 1:3
  our_start_params <- c(our_start_params, "m" = 5)

  # tidy date and calculate gradient for F1
  input <- dat %>%
    dplyr::filter(deme1 != deme2) %>%
    dplyr::filter(deme1 == 1 | deme2 == 1)


  input <- input %>%
    dplyr::mutate(gendist = ifelse(gendist > 0.999, 0.999,
                                   ifelse(gendist < 0.001, 0.001,
                                          gendist)))

  # now run model
  inputdisc <- dat %>%
    dplyr::filter(deme1 != deme2)
  ret <- disc(discdat = inputdisc,
              start_params = our_start_params,
              learningrate = 1e-3,
              b1 = 0.9,
              b2 = 0.999,
              e = 1e-8,
              steps = 1e3,
              normalize_geodist = T,
              report_progress = T,
              return_verbose = T)
  # back out for m adam
  b1 <- 0.9
  b2 <- 0.999
  e <- 1e-8
  learningrate <- 1e-3
  mt_m <- b1 * 0 + (1-b1) * ret$m_gradtraj[2]
  vt_m <- b2 * 0 + (1-b2) * (ret$m_gradtraj[2]^2)
  mt_mhat <- mt_m / (1 - b1^1)
  vt_mhat <- vt_m / (1 - b2^1)
  mnew1 <- exp( log(ret$m_run[1]) - learningrate * (mt_mhat/(sqrt(vt_mhat) + e)))

  # test out
  testthat::expect_equal(mnew1, ret$m_run[2])

})
