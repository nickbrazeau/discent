test_that("model runs works", {
  #......................
  # simulation some data
  #......................
  set.seed(48)
  K_gendist_geodist <- tibble::tibble(
    smpl1 = paste0("smpl", 1:26),
    smpl2 = paste0("smpl", 27:52),
    locat1 = sample(LETTERS, size = 26, replace = T),
    locat2 = sample(LETTERS, size = 26, replace = T),
    gendist = rbeta(26, 0.5, 0.5),
  )
  geodist <- K_gendist_geodist %>%
    dplyr::select(c("locat1", "locat2")) %>%
    dplyr::filter(!duplicated(.))
  geodist <- geodist %>%
    dplyr::mutate(geodist = rlnorm(dplyr::n(), 3, 1),
                  geodist = ifelse(locat1 == locat2, 0, geodist))
  # locat names
  locatnames <- unique(c(geodist$locat1, geodist$locat2))
  # bring together
  K_gendist_geodist <- dplyr::left_join(K_gendist_geodist, geodist, by = c("locat1", "locat2"))

  start_params = rep(0.1, length(locatnames))
  names(start_params) <- locatnames
  start_params <- c(start_params, "m" = 1e-4)
  #......................
  # run model
  #......................
  mod <- discent::deme_inbreeding_spcoef(K_gendist_geodist = K_gendist_geodist,
                                         start_params = start_params,
                                         m_lowerbound = 0,
                                         m_upperbound = 1,
                                         f_learningrate = 1e-4,
                                         m_learningrate = 1e-10,
                                         full_matrix = FALSE,
                                         steps = 1e2,
                                         report_progress = TRUE)
  testthat::expect_length(mod, 4)
})
