
#' @title Wrapper for Simulated Annealer Optimal Start Parameters
#' @details going to assume one "model" search is suffice
#'
get_SA_wrapper_start <- function(discdat,
                                 initf_learningrate = 1e-5,
                                 initm_learningrate = 1e-10,
                                 momentum = 0.9,
                                 discsteps = 5e2,

                                 initTemp = 1,
                                 annealstep = 1e4,
                                 SLratio = 0.8,
                                 FMSratio = 0.8,
                                 FMLratio = 0.5,
                                 demeSwitchSize = 3,

                                 fstartmin = 0.1,
                                 fstartmax = 0.9,
                                 fstartsteps = 100,

                                 mstartmin = 1e-6,
                                 mstartmax = 1,
                                 mstartsteps = 100,

                                 flearnmin = 1e-15,
                                 flearnmax = 1e-2,
                                 flearnsteps = 100,

                                 mlearnmin = 1e-20,
                                 mlearnmax = 1e-8,
                                 mlearnsteps = 100) {

  #............................................................
  # get modnames
  #...........................................................
  modnms <- unique(discdat$modname)
  # get first one
  SAdat <- lapply(modnms, function(x){
    rw <- min(which(x == discdat$modname))
    return(discdat[rw,])
  }) %>%
    dplyr::bind_rows()

  #............................................................
  # run SA
  #...........................................................
  # remove internal comparisons for DISCent algorithm
  rm_internal_compar <- function(dat){
    dat %>%
      dplyr::filter(deme1 != deme2)
  }
  # start params
  our_start_params <- rep(0.2, 7)
  names(our_start_params) <- c("1", "111", "58", "61", "64", "11", "121")
  our_start_params <- c(our_start_params, "m" = 1e-5)
  # run
  SAdat <- SAdat %>%
    dplyr::mutate(discdat = purrr::map(discdat, rm_internal_compar)) %>%
    dplyr::mutate(SAoptimparms = purrr::map(discdat,
                                                   discent::find_grad_params,
                                                   initstart_params = our_start_params,
                                                   initf_learningrate = initf_learningrate,
                                                   initm_learningrate = initm_learningrate,
                                                   momentum = momentum,
                                                   discsteps = discsteps,
                                                   initTemp = initTemp,
                                                   annealstep = annealstep,
                                                   SLratio = SLratio,
                                                   FMSratio = FMSratio,
                                                   FMLratio = FMLratio,
                                                   demeSwitchSize = demeSwitchSize,

                                                   fstartmin = fstartmin,
                                                   fstartmax = fstartmax,
                                                   fstartsteps = fstartsteps,

                                                   mstartmin = mstartmin,
                                                   mstartmax = mstartmax,
                                                   mstartsteps = mstartsteps,

                                                   flearnmin = flearnmin,
                                                   flearnmax = flearnmax,
                                                   flearnsteps = flearnsteps,

                                                   mlearnmin = mlearnmin,
                                                   mlearnmax = mlearnmax,
                                                   mlearnsteps = mlearnsteps)
    )
  # out
  SAdat <- SAdat %>%
    dplyr::select(c("modname", "SAoptimparms"))
  return(SAdat)
}



get_discentwrapper <- function(discdat = discdat, SAoptimparms = SAoptimparms) {
  # get start params
  discdatfull <- dplyr::left_join(discdat, SAoptimparms, by = "modname")

  # discent wrapper
  run_discent <- function(discdat, SAoptimparms) {
    discent::deme_inbreeding_spcoef(discdat = discdat,
                                    start_params = SAoptimparms$optimalParams$start_params,
                                    f_learningrate = SAoptimparms$optimalParams$f_learn,
                                    m_learningrate = SAoptimparms$optimalParams$m_learn,
                                    momentum = 0.9,
                                    steps = 1e5,
                                    report_progress = FALSE,
                                    return_verbose = FALSE)
  }
  # remove internal comparisons
  rm_internal_compar <- function(dat){
    dat %>%
      dplyr::filter(deme1 != deme2)
  }

  # run disc
  out <- discdatfull %>%
    dplyr::mutate(discdat = purrr::map(discdat, rm_internal_compar)) %>%
    dplyr::mutate(discret = furrr::future_map2(discdat, SAoptimparms,
                                               run_discent))
  return(out)
}
