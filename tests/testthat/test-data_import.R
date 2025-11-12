test_that("test that genetic data imports correctly", {
  # sim data
  data("IBD_simulation_data", package = "discent")
  dat <- IBD_simulation_data %>%
    dplyr::filter(deme1 != deme2)

  # start params
  our_start_params <- rep(0.2, 3)
  names(our_start_params) <- 1:3
  our_start_params <- c(our_start_params, "m" = 10)
  # locats
  locats <- names(our_start_params)[!grepl("^m$", names(our_start_params))]
  #.................................
  # input into Cpp
  #.................................
  disclist <- wrangle_discentdat(discdat = dat, start_params = our_start_params,
                                 locats = locats, normalize_geodist = T)

  #......................
  # model for export
  #......................
  mod <- disc(discdat = dat,
              start_params = our_start_params,
              learningrate = 1e-5,
              lambda = 1e-1,
              b1 = 0.9,
              b2 = 0.999,
              e = 1e-8,
              steps = 1e3,
              normalize_geodist = T,
              report_progress = F,
              return_verbose = T)
  #............................................................
  # extract
  #...........................................................
  # in
  inputgendat <- disclist$gendist_arr
  inputgeodat <- disclist$geodist_mat
  # out
  outputgendat <- mod$raw_gendist_arr
  outputgeodat <- mod$raw_geodist_mat

  #............................................................
  # COMPARE geodistances
  #...........................................................
  testthat::expect_true(
    identical(as.vector(inputgeodat), unlist(outputgeodat)) # R is column-major, whilst cpp is row-major. This works bc it is a symmetric matrix
  )

  #............................................................
  # COMPARE gendistances
  #...........................................................
  gen_compare <- array(-1.0, dim = c(length(disclist$demes), length(disclist$demes), disclist$n_Kpairmax))

  for (i in seq_along(disclist$demes)) {
    for (j in seq_along(disclist$demes)) {
      for (k in 1:disclist$n_Kpairmax) {
        gen_compare[i,j,k] <- outputgendat[[i]][[j]][k]
      }
    }
  }

  testthat::expect_true(
    identical(inputgendat, gen_compare)
  )


  #............................................................
  # COMPARE raw data
  #...........................................................
  # this is the nesting code from the wrangle function
  demes <- sort(unique(c(dat$deme1, dat$deme2)))
  keyi <- data.frame(deme1 = demes, i = seq_along(demes))
  keyj <- data.frame(deme2 = demes, j = seq_along(demes))

  # this makes an ij key independent of deme names
  gendist <- dat %>%
    expand_pairwise(.) %>% # get all pairwise for full matrix
    dplyr::select(c("deme1", "deme2", "gendist")) %>%
    dplyr::filter(!duplicated(.)) %>%
    dplyr::group_by_at(c("deme1", "deme2")) %>%
    tidyr::nest(.) %>%
    dplyr::left_join(., keyi, by = "deme1") %>%
    dplyr::left_join(., keyj, by = "deme2") %>%
    dplyr::arrange_at(c("i", "j"))

  # the input matrix then becomes [i,j,k] = [deme1, deme2, pair]
  # random checks
  testthat::expect_true(
    identical(as.numeric(gendist$data[[3]][9,1]), disclist$gendist_arr[2,1,9]) )
  testthat::expect_true(
    identical(as.numeric(gendist$data[[2]][12,1]), disclist$gendist_arr[1,3,12]) )

  testthat::expect_true(
    identical(as.numeric(gendist$data[[4]][1,1]), disclist$gendist_arr[2,3,1]) )

  testthat::expect_true(
    identical(as.numeric(gendist$data[[4]][9,1]), disclist$gendist_arr[2,3,9]) )

  testthat::expect_true(
    identical(as.numeric(gendist$data[[6]][5,1]), disclist$gendist_arr[3,2,5]) )


})



