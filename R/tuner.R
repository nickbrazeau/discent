#' @title Find Optimal Start Gradient Parameters using Grid Search
#' @inheritParams deme_inbreeding_spcoef
#' @param discsteps integer; the number of "steps" in the main `DISCent`
#' algorithm" \link{deme_inbreeding_spcoef}
#' @param temptstart_params named double vector; vector of start parameters to serve as template.
#'
#' @param fstartmin double; the minimum value for F inbreeding coefficient in the search grid
#' @param fstartmax double; the maximum value for F inbreeding coefficient in the search grid
#' @param fstartby integer; the number of values between fstartmin and fstartmax to evaluate
#'
#' @param mstartmin double; the minimum value for M global migration rate in the search grid
#' @param mstartmax double; the maximum value for M global migration rate in the search grid
#' @param mstartby integer; the number of values between mstartmin and mstartmax to evaluate
#'
#' @param flearnmin double; the minimum value for F learning rate in the search grid
#' @param flearnmax double; the maximum value for F learning rate in the search grid
#' @param flearnby integer; the number of values between flearnmin and flearnmax to evaluate
#'
#' @param mlearnmin double; the minimum value for M learning rate in the search grid
#' @param mlearnmax double; the maximum value for M learning rate in the search grid
#' @param mlearnby integer; the number of values between mlearnmin and mlearnmax to evaluate
#'
#' @param downsample integer; the subset, or the number of a smaller random set of values, that you want to search in the combination grid
#'
#' @details Using a discrete grid search, identify the best start parameters for
#' inbreeding F values and M global migration parameter as well as the best
#' learning rates. We assume that each deme should have the same F start value for computational tractability.
#'
#' The final cost for each set of start parameters is given by the gradient descent cost after _n_ iterations of
#' the main `DISCent` algorithm" \link{deme_inbreeding_spcoef}.
#' The search grid is defined by the user specified min, max, and steps ("by") of the
#' aforementioned parameters.
#' @export

find_grad_params <- function(discdat,
                             tempstart_params = c(),
                             momentum = 0.9,
                             discsteps = 1e3,

                             fstartmin = 0.1,
                             fstartmax = 0.9,
                             fstartby = 100,

                             mstartmin = 1e-15,
                             mstartmax = 1e-3,
                             mstartby = 100,

                             flearnmin = 1e-15,
                             flearnmax = 1e-2,
                             flearnby = 100,

                             mlearnmin = 1e-15,
                             mlearnmax = 1e-2,
                             mlearnby = 100) {

  #..............................................................
  # Assertions & Catches
  #..............................................................
  if (!all(colnames(discdat) %in% c("smpl1", "smpl2", "deme1", "deme2", "gendist", "geodist"))) {
    stop("The discdat dataframe must contain columns with names: smpl1, smpl2, deme1, deme2, gendist, geodist")
  }
  # make sure correct order
  for (i in 1:6) {
    if (colnames(discdat)[i] != c("smpl1", "smpl2", "deme1", "deme2", "gendist", "geodist")[i]) {
      stop("The discdat dataframe must contain columns with names in the exact order of: smpl1, smpl2, deme1, deme2, gendist, geodist")
    }
  }

  locats <- names(tempstart_params)[!grepl("^m$", names(tempstart_params))]
  if (!all(unique(c(discdat$deme1, discdat$deme2)) %in% locats)) {
    stop("You have cluster names in your discdat dataframe that are not included in your start parameters")
  }
  if (!any(grepl("^m$", names(tempstart_params)))) {
    stop("A m start parameter must be provided (i.e. there must be a vector element name m in your start parameter vector)")
  }
  assert_dataframe(discdat)
  assert_vector(tempstart_params)
  assert_length(tempstart_params, n = (length(unique(c(discdat$deme1, discdat$deme2))) + 1),
                message = "Start params length not correct. You must specificy a start parameter
                           for each deme and the migration parameter, m")
  sapply(tempstart_params[!grepl("^m$", names(tempstart_params))], assert_bounded, left = 0, right = 1, inclusive_left = TRUE, inclusive_right = TRUE)
  assert_single_numeric(momentum)
  assert_single_int(discsteps)

  assert_single_int(downsample)

  assert_single_numeric(fstartmin)
  assert_single_numeric(fstartmax)
  assert_single_int(fstartby)

  assert_single_numeric(mstartmin)
  assert_single_numeric(mstartmax)
  assert_single_int(mstartby)

  assert_single_numeric(flearnmin)
  assert_single_numeric(flearnmax)
  assert_single_int(flearnby)

  assert_single_numeric(mlearnmin)
  assert_single_numeric(mlearnmax)
  assert_single_int(mlearnby)

  #............................................................
  # look up tables
  #...........................................................
  fstartsvec <- seq(from = fstartmin, to = fstartmax, length.out = fstartby)
  mstartsvec <- seq(from = mstartmin, to = mstartmax, length.out = mstartby)
  flearnsvec <- seq(from = flearnmin, to = flearnmax, length.out = flearnby)
  mlearnsvec <- seq(from = mlearnmin, to = mlearnmax, length.out = mlearnby)
  search_grid <- expand.grid(fstartsvec, mstartsvec, flearnsvec, mlearnsvec)
  colnames(search_grid) <- c("fstarts", "mstart", "f_learn", "m_learn")
  # check downsmaple
  assert_bounded(downsample, left = 1, right = nrow(search_grid), inclusive_left = TRUE, inclusive_right = TRUE,
                 message = "Downsample must consider at least 1 search grid parameter, or be less than or equal to the total number of potential search grid combinations")
  # down sample
  search_grid <- search_grid[sample(1:nrow(search_grid), size  = downsample, replace = F), ]

  # liftover to start param format
  liftover_start_params <- function(fstarts, mstarts, start_param_template) {
    out <- start_param_template
    out[names(out) != "m"] <- fstarts
    out[names(out) == "m"] <- mstarts
    return(out)
  }
  search_grid <- search_grid %>%
    dplyr::mutate(start_params = purrr::map2(fstarts, mstarts, liftover_start_params,
                                             start_param_template = tempstart_params)) %>%
    dplyr::select(c("start_params", "f_learn", "m_learn"))

  #......................
  # wrap discent
  #......................
  discent_wrapper <- function(discdat,
                              start_params,
                              f_learningrate,
                              m_learningrate,
                              momentum,
                              steps){

    out <- discent::deme_inbreeding_spcoef(discdat,
                                           start_params,
                                           f_learningrate,
                                           m_learningrate,
                                           momentum,
                                           steps,
                                           report_progress = FALSE,
                                           return_verbose = FALSE)$cost[discsteps]
    return(out)
  }

  search_grid$cost <- furrr::future_pmap_dbl(search_grid, discent_wrapper,
                                             discdat = discdat,
                                             steps = discsteps,
                                             momentum = momentum)


  # out with search grid result
  return(search_grid)
}




