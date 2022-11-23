

#' @title Wrapper for Search Grid Optimal Start Parameters
#' @details going to assume one "model" search is suffice
#'
get_GS_cost <- function(searchgriddf, discdat) {
  cost <- discent::deme_inbreeding_spcoef(discdat = discdat$discdat,
                                  start_params = searchgriddf$start_params,
                                  f_learningrate = searchgriddf$f_learn,
                                  m_learningrate = searchgriddf$m_learn,
                                  momentum = 0.9,
                                  steps = 1e4,
                                  report_progress = FALSE,
                                  return_verbose = FALSE)$cost[1e4]

  searchgriddf$cost <- cost
  searchgriddf$modname <- discdat$modname
  return(searchgriddf)
}



get_discentwrapper <- function(discdat = discdat, optimstarts = optimstarts) {
  # get start params


  # discent wrapper
  run_discent <- function(discdat, optimstarts) {
    discent::deme_inbreeding_spcoef(discdat = discdat,
                                    start_params = optimstarts$start_params,
                                    f_learningrate = optimstarts$f_learn,
                                    m_learningrate = optimstarts$m_learn,
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
    dplyr::mutate(discret = furrr::future_map2(discdat, optimstarts,
                                               run_discent))
  return(out)
}
