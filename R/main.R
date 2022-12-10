#' @title Identify Deme Inbreeding Spatial Coefficients in Continuous Space
#' @param discdat dataframe; The genetic-geographic data by deme (K)
#' @param start_params named double vector; vector of start parameters.
#' @param f_learningrate double; alpha parameter for how much each "step" is weighted in the gradient descent for inbreeding coefficients
#' @param m_learningrate double; alpha parameter for how much each "step" is weighted in the gradient descent for the migration parameter
#' @param momentum double; gamma parameter for momentum for adaptive learning rate
#' @param steps integer; the number of "steps" as we move down the gradient
#' @param m_lowerbound double; lower limit value for the global "m" parameter; will use a reflected normal within the gradient descent algorithm to adjust any aberrant values
#' @param m_upperbound double; upper limit value for the global "m" parameter; will use a reflected normal within the gradient descent algorithm to adjust any aberrant values
#' @param report_progress boolean; whether or not a progress bar should be shown as you iterate through steps
#' @param return_verbose boolean; whether the inbreeding coefficients and migration rate should be returned for every iteration or
#' only for the final iteration. User will typically not want to store every iteration, which can be memory intensive
#' @details The gen.geo.dist dataframe must be named with the following columns:
#'          "smpl1"; "smpl2"; "deme1"; "deme2"; "gendist"; "geodist"; which corresponds to:
#'          Sample 1 Name; Sample 2 Name; Sample 1 Location; Sample 2 Location;
#'          Pairwise Genetic Distance; Pairwise Geographpic Distance. Note, the order of the
#'          columns do not matter but the names of the columns must match.
#' @details The start_params vector names must match the cluster names (i.e. clusters must be
#'          have a name that we can match on for the starting relatedness paramerts). In addition,
#'          you must provide a start parameter for "m".
#' @description The purpose of this statistic is to identify an inbreeding coefficient, or degree of
#'              relatedness, for a given location in space. We assume that locations in spaces can be
#'              represented as "demes," such that multiple individuals live in the same deme
#'              (i.e. samples are sourced from the same location). The expected pairwise relationship
#'              between two individuals, or samples, is dependent on the each sample's deme's inbreeding
#'              coefficient and the geographic distance between the demes. The program assumes a symmetric distance matrix.
#' @details Note: We have implemented coding decisions to not allow the "f" inbreeding coefficients to be negative by using a
#' logit transformation internally in the code. Similarly, the "m" global migration rate cannot be negative and is bounded to be positive by
#' through taking the absolute value in every iteration.
#' @export
#'

deme_inbreeding_spcoef <- function(discdat,
                                   start_params = c(),
                                   f_learningrate = 1e-5,
                                   m_learningrate = 1e-10,
                                   m_lowerbound = 0,
                                   m_upperbound = Inf,
                                   momentum = 0.9,
                                   steps = 1e3,
                                   report_progress = TRUE,
                                   return_verbose = FALSE){

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

  locats <- names(start_params)[!grepl("^m$", names(start_params))]
  if (!all(unique(c(discdat$deme1, discdat$deme2)) %in% locats)) {
    stop("You have cluster names in your discdat dataframe that are not included in your start parameters")
  }
  if (!any(grepl("^m$", names(start_params)))) {
    stop("A m start parameter must be provided (i.e. there must be a vector element name m in your start parameter vector)")
  }
  assert_dataframe(discdat)
  assert_vector(start_params)
  assert_length(start_params, n = (length(unique(c(discdat$deme1, discdat$deme2))) + 1),
                message = "Start params length not correct. You must specificy a start parameter
                           for each deme and the migration parameter, m")
  sapply(start_params[!grepl("^m$", names(start_params))], assert_bounded, left = 0, right = 1, inclusive_left = TRUE, inclusive_right = TRUE)
  assert_single_numeric(f_learningrate)
  assert_single_numeric(m_learningrate)
  assert_single_numeric(momentum)
  assert_single_numeric(m_lowerbound)
  assert_single_numeric(m_upperbound)
  assert_gr(m_upperbound, m_lowerbound)
  assert_single_int(steps)
  assert_single_logical(report_progress)

  # no missing
  if(sum(is.na(discdat)) != 0) {
    stop("discdat dataframe cannot have missing values")
  }

  #......................
  # check for self comparisons
  #......................
  sapply(discdat$geodist, assert_neq, y = 0,
         message = "No within-deme sample comparisons allowed. Geodistance should not be 0")
  mapply(assert_neq, discdat$deme1, discdat$deme2,
         message = "No within-deme sample comparisons allowed. Locat names should not be the same")


  #..............................................................
  # setup and create progress bars
  #..............................................................
  pb <- utils::txtProgressBar(min = 0, max = steps-1, initial = NA, style = 3)
  args_progress <- list(pb = pb)


  #............................................................
  # R manipulations before C++
  #...........................................................
  # use efficient R functions to group pairs and wrangle data for faster C++ manipulation
  # get deme names and lift over sorted names for i and j
  demes <- sort(unique(c(discdat$deme1, discdat$deme2)))
  keyi <- data.frame(deme1 = demes, i = 1:length(demes))
  keyj <- data.frame(deme2 = demes, j = 1:length(demes))

  # transform data w/ logit
  discdat <- discdat %>%
    dplyr::mutate(gendist = ifelse(gendist > 0.999, 0.999,
                                   ifelse(gendist < 0.001, 0.001,
                                          gendist))) %>% # reasonable bounds on logit
    dplyr::mutate(gendist = logit(gendist))
  # transform start parameters w/ logit
  start_params[names(start_params) != "m"] <- logit(start_params[names(start_params) != "m"])


  # get genetic data by pairs through efficient nest
  gendist <- discdat %>%
    expand_pairwise(.) %>% # get all pairwise for full matrix
    dplyr::select(c("deme1", "deme2", "gendist")) %>%
    dplyr::group_by_at(c("deme1", "deme2")) %>%
    tidyr::nest(.) %>%
    dplyr::left_join(., keyi, by = "deme1") %>%
    dplyr::left_join(., keyj, by = "deme2") %>%
    dplyr::arrange_at(c("i", "j"))


  # put gendist into an array
  # NB we are filling an array with dimension of size:
  #   locations x locations x max K-pairs
  # Will fill in w/ -1 to indicate missing/skip where demes do not
  # have max pairs.
  # This approach wastes memory but allows for a structured array
  # versus a list with varying sizes (and eventually a more efficient for-loop)
  n_Kpairmax <- max(purrr::map_dbl(gendist$data, nrow))
  gendist_arr <- array(data = -1, dim = c(length(locats), length(locats), n_Kpairmax))
  for (i in 1:nrow(gendist)) {
    gendist_arr[gendist$i[i], gendist$j[i], 1:nrow(gendist$data[[i]])] <- unname(unlist(gendist$data[[i]]))
  }

  # put geo information into distance matrix
  geodist <- discdat %>%
    expand_pairwise(.) %>% # get all pairwise for full matrix
    dplyr::select(c("deme1", "deme2", "geodist")) %>%
    dplyr::group_by_at(c("deme1", "deme2")) %>%
    tidyr::nest(.) %>%
    dplyr::left_join(., keyi, by = "deme1") %>%
    dplyr::left_join(., keyj, by = "deme2") %>%
    dplyr::arrange_at(c("i", "j"))

  # simplify geodistance data storage
  geodist$data <- purrr::map_dbl(geodist$data, function(x){
    if (length(unique(unlist(x))) != 1) {
      stop("deme1 and deme2 have different geodistances among P-sample combinations. Distances should all be same among samples")
    }
    return( unique(unlist(x)) ) # all same by unique
  }
  )

  # upper tri
  geodist_mat <- matrix(data = -1, nrow = length(locats), ncol = length(locats))
  for (i in 1:nrow(geodist)) {
    geodist_mat[geodist$i[i], geodist$j[i]] <- geodist$data[i]
  }
  diag(geodist_mat) <- 0

  #..............................................................
  # run efficient C++ function
  #..............................................................

  args <- list(gendist = as.vector(gendist_arr),
               geodist = as.vector(geodist_mat),
               fvec = unname( start_params[!grepl("^m$", names(start_params))] ),
               n_Kpairmax = n_Kpairmax,
               m = unname(start_params["m"]),
               f_learningrate = f_learningrate,
               m_learningrate = m_learningrate,
               m_lowerbound = m_lowerbound,
               m_upperbound = m_upperbound,
               momentum = momentum,
               steps = steps,
               report_progress = report_progress
  )

  # create progress bars
  pb <- txtProgressBar(min = 0, max = steps, initial = NA, style = 3)
  args_progress <- list(pb = pb)
  args_functions <- list(update_progress = update_progress)
  output_raw <- deme_inbreeding_coef_cpp(args, args_functions, args_progress)

  # process output
  colnames(keyi) <- c("Deme", "key")
  if (return_verbose) {
    output <- list(
      deme_key = keyi,
      m_run = output_raw$m_run,
      fi_run = expit(do.call("rbind", output_raw$fi_run)),
      m_update = output_raw$m_update,
      fi_update = do.call("rbind", output_raw$fi_update),
      cost = output_raw$cost,
      Final_Fis = expit(output_raw$Final_Fis),
      Final_m = output_raw$Final_m
    )

  } else {
    output <- list(
      deme_key = keyi,
      cost = output_raw$cost,
      Final_Fis = expit(output_raw$Final_Fis),
      Final_m = output_raw$Final_m)
  }

  # return list
  return(output)

}


#' @title Find Optimal Start Gradient Parameters using Grid Search
#' @inheritParams deme_inbreeding_spcoef
#' @param discsteps integer; the number of "steps" in the main `DISCent`
#' algorithm" \link{deme_inbreeding_spcoef}
#' @param tempstart_params named double vector; vector of start parameters to serve as template.
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
                             downsample = 1e3,

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
  assert_bounded(downsample, left = 1, right = nrow(search_grid), inclusive_left = TRUE, inclusive_right = TRUE)

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

    out <- deme_inbreeding_spcoef(discdat,
                                  start_params,
                                  f_learningrate,
                                  m_learningrate,
                                  momentum,
                                  steps,
                                  report_progress = FALSE,
                                  return_verbose = FALSE)$cost[discsteps]
    return(out)
  }

  search_grid$cost <- purrr::pmap_dbl(search_grid, discent_wrapper,
                                      discdat = discdat,
                                      steps = discsteps,
                                      momentum = momentum)


  # out with search grid result
  return(search_grid)
}
