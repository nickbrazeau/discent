test_that("model runs works", {
  #......................
  # simulation some data
  #......................
  set.seed(48)
  locats <- tibble::tibble(
    smpl = paste0("smpl", 1:10),
    locat = sample(LETTERS, size = 10, replace = F)
  )
  # bring together
  K_gendist_geodist <- as.data.frame(t(combn(locats$smpl, 2)))
  colnames(K_gendist_geodist)[1] <- "smpl1"
  colnames(locats)[1] <- "smpl1"
  K_gendist_geodist <- dplyr::left_join(K_gendist_geodist, locats, by = "smpl1")
  colnames(K_gendist_geodist)[2] <- "smpl2"
  colnames(locats)[1] <- "smpl2"
  K_gendist_geodist <- dplyr::left_join(K_gendist_geodist, locats, by = "smpl2")
  colnames(K_gendist_geodist) <- c("smpl1", "smpl2", "locat1", "locat2")

  # now remove selfs (not necessary since we sample w/out replacement above, but for transparency)
  K_gendist_geodist <- K_gendist_geodist %>%
    dplyr::filter(locat1 != locat2)

  # make up some genetic distance
  K_gendist_geodist$gendist <- rbeta(nrow(K_gendist_geodist), 0.5, 0.5)

  geodist <- K_gendist_geodist %>%
    dplyr::select(c("locat1", "locat2")) %>%
    dplyr::filter(!duplicated(.))
  geodist <- geodist %>%
    dplyr::mutate(geodist = rlnorm(dplyr::n(), 3, 1))
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
                                         steps = 1e2,
                                         report_progress = TRUE)
  testthat::expect_length(mod, 4)
})
