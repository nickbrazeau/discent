#' @title Identify Deme Inbreeding Spatial Coefficients in Continuous Space
#' @param K_gendist_geodist dataframe; The genetic-geographic data by deme (K)
#' @param start_params named numeric vector; vector of start parameters.
#' @param f_learningrate numeric; alpha parameter for how much each "step" is weighted in the gradient descent for inbreeding coefficients
#' @param m_learningrate numeric; alpha parameter for how much each "step" is weighted in the gradient descent for the migration parameter
#' @param momentum numeric; gamma parameter for momentum for adaptive learning rate
#' @param steps numeric; the number of "steps" as we move down the gradient
#' @param report_progress boolean; whether or not a progress bar should be shown as you iterate through steps
#' @param return_verbose boolean; whether the inbreeding coefficients and migration rate should be returned for every iteration or
#' only for the final iteration. User will typically not want to store every iteration, which can be memory intensive
#' @details The gen.geo.dist dataframe must be named with the following columns:
#'          "smpl1"; "smpl2"; "locat1"; "locat2"; "gendist"; "geodist"; which corresponds to:
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
#' @export
#'

deme_inbreeding_spcoef <- function(K_gendist_geodist,
                                   start_params = c(),
                                   f_learningrate = 1e-5,
                                   m_learningrate = 1e-10,
                                   momentum = 0.9,
                                   steps = 1e3,
                                   report_progress = TRUE,
                                   return_verbose = FALSE){

  #..............................................................
  # Assertions & Catches
  #..............................................................
  if (!all(colnames(K_gendist_geodist) %in% c("smpl1", "smpl2", "locat1", "locat2", "gendist", "geodist"))) {
    stop("The K_gendist_geodist must contain columns with names: smpl1, smpl2, locat1, locat2, gendist, geodist")
  }
  # make sure correct order
  for (i in 1:6) {
    if (colnames(K_gendist_geodist)[i] != c("smpl1", "smpl2", "locat1", "locat2", "gendist", "geodist")[i]) {
      stop("The K_gendist_geodist must contain columns with names in the exact order of: smpl1, smpl2, locat1, locat2, gendist, geodist")
    }
  }

  locats <- names(start_params)[!grepl("m", names(start_params))]
  if (!all(unique(c(K_gendist_geodist$locat1, K_gendist_geodist$locat2)) %in% locats)) {
    stop("You have cluster names in your K_gendist_geodist dataframe that are not included in your start parameters")
  }
  if (!any(grepl("m", names(start_params)))) {
    stop("A m start parameter must be provided (i.e. there must be a vector element name m in your start parameter vector)")
  }
  assert_dataframe(K_gendist_geodist)
  assert_vector(start_params)
  assert_length(start_params, n = (length(unique(c(K_gendist_geodist$locat1, K_gendist_geodist$locat2))) + 1),
                message = "Start params length not correct. You must specificy a start parameter
                           for each deme and the migration parameter, m")
  sapply(start_params[!grepl("m", names(start_params))], assert_bounded, left = 0, right = 1, inclusive_left = TRUE, inclusive_right = TRUE)
  assert_single_numeric(f_learningrate)
  assert_single_numeric(m_learningrate)
  assert_single_numeric(momentum)
  assert_single_int(steps)
  assert_single_logical(report_progress)

  #......................
  # check for self comparisons
  #......................
  sapply(K_gendist_geodist$geodist, assert_neq, y = 0,
         message = "No within-deme sample comparisons allowed. Geodistance should not be 0")
  mapply(assert_neq, K_gendist_geodist$locat1, K_gendist_geodist$locat2,
         message = "No within-deme sample comparisons allowed. Locat names should not be the same")


  #..............................................................
  # setup and create progress bars
  #..............................................................
  # note, samples names there just for sanity
  pb <- txtProgressBar(min = 0, max = steps, initial = NA, style = 3)
  args_progress <- list(pb = pb)


  #............................................................
  # R manipulations before C++
  #...........................................................
  # use efficient R functions to group pairs and wrangle data for faster C++ manipulation
  # get deme names and lift over sorted names for i and j
  demes <- sort(unique(c(K_gendist_geodist$locat1, K_gendist_geodist$locat2)))
  keyi <- data.frame(locat1 = demes, i = 1:length(demes))
  keyj <- data.frame(locat2 = demes, j = 1:length(demes))


  # transform data
  K_gendist_geodist <- K_gendist_geodist %>%
    dplyr::mutate(gendist = discent:::logit(gendist),
                  gendist = ifelse(gendist == Inf, 6, gendist), # reasonable bounds on logit
                  gendist = ifelse(gendist == -Inf, -6, gendist) # reasonable bounds on logit
    )


  # get genetic data by pairs through efficient nest
  gendist <- K_gendist_geodist %>%
    discent:::expand_pairwise(.) %>% # get all pairwise for full matrix
    dplyr::select(c("locat1", "locat2", "gendist")) %>%
    dplyr::group_by_at(c("locat1", "locat2")) %>%
    tidyr::nest(.) %>%
    dplyr::left_join(., keyi, by = "locat1") %>%
    dplyr::left_join(., keyj, by = "locat2") %>%
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
  geodist <- K_gendist_geodist %>%
    discent:::expand_pairwise(.) %>% # get all pairwise for full matrix
    dplyr::select(c("locat1", "locat2", "geodist")) %>%
    dplyr::group_by_at(c("locat1", "locat2")) %>%
    tidyr::nest(.) %>%
    dplyr::left_join(., keyi, by = "locat1") %>%
    dplyr::left_join(., keyj, by = "locat2") %>%
    dplyr::arrange_at(c("i", "j"))

  # simplify geodistance data storage
  geodist$data <- purrr::map_dbl(geodist$data, function(x){
    if (length(unique(unlist(x))) != 1) {
      stop("Locat1 and Locat2 have different geodistances among P-sample combinations. Distances should all be same among samples")
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
               fvec = unname( start_params[!grepl("m", names(start_params))] ),
               n_Kpairmax = n_Kpairmax,
               m = unname(start_params["m"]),
               f_learningrate = f_learningrate,
               m_learningrate = m_learningrate,
               momentum = momentum,
               steps = steps,
               report_progress = report_progress
  )

  # create progress bars
  pb <- txtProgressBar(min = 0, max = steps-1, initial = NA, style = 3)
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
      fi_ada = output_raw$fi_ada_grad,
      m_ada = output_raw$m_ada_grad,
      bGf = output_raw$bGf,
      bGm = output_raw$bGm,
      store_f_learn = output_raw$store_f_learn,
      store_m_learn = output_raw$store_m_learn,
      fi_store_grad = do.call("rbind", output_raw$fi_store_grad),
      m_store_grad = output_raw$m_store_grad,
      cost = output_raw$cost,
      Final_Fis = expit(output_raw$Final_Fis),
      Final_m = output_raw$Final_m)
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
