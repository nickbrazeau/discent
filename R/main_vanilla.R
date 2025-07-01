#' @title Identify Deme Inbreeding Spatial Coefficients in Continuous Space
#' @param discdat dataframe; The genetic-geographic data by deme (K)
#' @param start_params named double vector; vector of start parameters.
#' @param lambda double; A quadratic L2 explicit regularization, or penalty, parameter on "m" parameter. Note, lambda is a scalar such that: \eqn{\lambda m^2}.
#' @param learningrate double; alpha parameter for how much each "step" is weighted in the gradient descent
#' @param b1 double; exponential decay rates for the first moment estimate in the Adam optimization algorithm
#' @param b2 double; exponential decay rates for the second moment estimate in the Adam optimization algorithm
#' @param e double; epsilon (error) for stability in the Adam optimization algorithm
#' @param steps integer; the number of steps as we move down the gradient
#' @param thin integer; the number of steps to keep as part of the output (i.e. if the user specifies 10, every 10th iteration will be kept)
#' @param m_lowerbound double; lower limit value for the global "m" parameter; any "m" value encounter less than the lower bound will be replaced by the lower bound
#' @param m_upperbound double; upper limit value for the global "m" parameter; any "m" value encounter greater than the upper bound will be replaced by the upper bound
#' @param normalize_geodist boolean; whether geographic distances between demes should be normalized (i.e. Min-Max Feature Scaling: \eqn{X' = \frac{X - X_{min}}{X_{max} - X_{min}} }, which places the geodistances on the scale to \eqn{[0-1]}). Helps increase model stability at the expense of complicating the interpretation of the migration rate parameter.
#' @param report_progress boolean; whether or not a progress bar should be shown as you iterate through steps
#' @param return_verbose boolean; whether the inbreeding coefficients and migration rate should be returned for every iteration or
#' only for the final iteration. User will typically not want to store every iteration, which can be memory intensive
#' @description The purpose of this statistic is to identify an inbreeding coefficient, or degree of
#'              relatedness, for a given location in discrete space. We assume that locations in spaces can be
#'              represented as "demes," such that multiple individuals live in the same deme
#'              (i.e. samples are sourced from the same location). The expected pairwise relationship
#'              between two individuals, or samples, is dependent on the each sample's deme's inbreeding
#'              coefficient and the geographic distance between the demes. The program assumes a symmetric distance matrix.
#' @details The gen.geo.dist dataframe must be named with the following columns:
#'          "smpl1"; "smpl2"; "deme1"; "deme2"; "gendist"; "geodist"; which corresponds to:
#'          Sample 1 Name; Sample 2 Name; Sample 1 Location; Sample 2 Location;
#'          Pairwise Genetic Distance; Pairwise Geographpic Distance. Note, the order of the
#'          columns do not matter but the names of the columns must match.
#' @details The start_params vector names must match the cluster names (i.e. clusters must be
#'          have a name that we can match on for the starting relatedness paramerts). In addition,
#'          you must provide a start parameter for "m".
#' @details Note: We have implemented coding decisions to not allow the "f" inbreeding coefficients to be negative by using a
#' logit transformation internally in the code.
#' @details Gradient descent is performed using the Adam (adaptive moment estimation) optimization approach. Default values
#' for moment decay rates, epsilon, and learning rates are taken from \cite{DP Kingma, 2014}.
#' @export

disc <- function(discdat,
                 start_params = NULL,
                 lambda = 0.1,
                 learningrate = 1e-3,
                 m_lowerbound = 0,
                 m_upperbound = Inf,
                 b1 = 0.9,
                 b2 = 0.999,
                 e = 1e-8,
                 steps = 1e3,
                 thin = 1,
                 normalize_geodist = TRUE,
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
  assert_single_numeric(lambda)
  assert_single_numeric(learningrate)
  assert_single_numeric(b1)
  assert_single_numeric(b2)
  assert_single_numeric(e)
  assert_single_numeric(m_lowerbound)
  assert_single_numeric(m_upperbound)
  assert_gr(m_upperbound, m_lowerbound)
  assert_single_int(steps)
  assert_single_int(thin)
  assert_greq(thin, 1, message = "Must be at least 1")
  assert_single_logical(report_progress)
  assert_single_logical(normalize_geodist)

  # no missing
  if(sum(is.na(discdat)) != 0) {
    stop("discdat dataframe cannot have missing values")
  }

  # catch accidental bad F and M start
  if ( any(round(start_params, digits = 1e200) == 0) ) {
    warning("At least one of your start parameters is zero (or essentially zero), which will result in unstable behavior in the Gradient-Descent algorithm. Consider increasing the start parameter.")
  }

  #......................
  # check that disc dat is being used correctly
  #......................
  assert_eq(
    x = choose(length(unique(c(discdat$deme1, discdat$deme2))), 2), # max number of between distances
    y = length(unique(paste(discdat$deme1, discdat$deme2, discdat$geodist, sep = "-"))), # count of unique between distances
    message = "You have pairwise demes with differing geodistances measurements in your data input. Geodistances must be the same between demes (pairwise geodistances should be the same)."
  )


  #......................
  # check for self comparisons
  #......................
  sapply(discdat$geodist, assert_neq, y = 0,
         message = "No within-deme sample comparisons allowed. Geodistance should not be 0")
  mapply(assert_neq, discdat$deme1, discdat$deme2,
         message = "No within-deme sample comparisons allowed. Locat names should not be the same")
  #............................................................
  # R manipulations before C++
  #...........................................................
  disclist <- wrangle_discentdat(discdat, normalize_geodist, start_params, locats)

  #..............................................................
  # run efficient C++ function
  #..............................................................

  args <- list(gendist = as.vector(disclist$gendist_arr),
               geodist = as.vector(disclist$geodist_mat),
               fvec = unname( start_params[!grepl("^m$", names(start_params))] ),
               n_Demes = length(disclist$demes),
               n_Kpairmax = disclist$n_Kpairmax,
               m = unname(start_params["m"]),
               lambda = lambda,
               learningrate = learningrate,
               m_lowerbound = m_lowerbound,
               m_upperbound = m_upperbound,
               b1 = b1,
               b2 = b2,
               e = e,
               steps = steps,
               report_progress = report_progress
  )

  # run
  output_raw <- vanilla_deme_inbreeding_coef_cpp(args)


  # set up thinning
  thin_its <- seq(1, steps, by = thin)
  thin_its <- unique(c(thin_its, steps)) # always include last iteration

  # process output
  colnames(disclist$keyi) <- c("Deme", "key")
  if (return_verbose) {
    output <- list(
      deme_key = disclist$keyi,
      m_run = output_raw$m_run[thin_its],
      fi_run = expit(do.call("rbind", output_raw$fi_run))[thin_its, ],
      m_gradtraj = output_raw$m_gradtraj[thin_its],
      fi_gradtraj = do.call("rbind", output_raw$fi_gradtraj)[thin_its, ],
      m_1moment = output_raw$m_firstmoment[thin_its],
      m_2moment = output_raw$m_secondmoment[thin_its],
      fi_1moment = do.call("rbind", output_raw$fi_firstmoment)[thin_its, ],
      fi_2moment = do.call("rbind", output_raw$fi_secondmoment)[thin_its, ],
      cost = output_raw$cost[thin_its],
      Final_Fis = expit(output_raw$Final_Fis),
      Final_m = output_raw$Final_m,
      raw_geodist_mat = output_raw$raw_geodist_mat,
      raw_gendist_arr = output_raw$raw_gendist_arr
    )

  } else {
    output <- list(
      deme_key = disclist$keyi,
      cost = output_raw$cost[thin_its],
      m_run = output_raw$m_run[thin_its],
      fi_run = expit(do.call("rbind", output_raw$fi_run))[thin_its, ],
      Final_Fis = expit(output_raw$Final_Fis),
      Final_m = output_raw$Final_m)
  }

  # add S3 class structure
  attr(output, "class") <- "DISCresult"
  return(output)
}


#' @title Wrangle DISC data for efficient input into Cpp
#' @inheritParams disc
#' @description Internal wrangling function that expands the DISC genetic data into an array and the geographic
#' distances into a matrix for efficient input into Cpp.
#' @export

wrangle_discentdat <- function(discdat, normalize_geodist, start_params, locats) {
  # use efficient R functions to group pairs and wrangle data for faster C++ manipulation
  # get deme names and lift over sorted names for i and j (note, deme names may be anything, so cannot rely on user indexing)
  demes <- sort(unique(c(discdat$deme1, discdat$deme2)))
  keyi <- data.frame(deme1 = demes, i = seq_len(length(demes)))
  keyj <- data.frame(deme2 = demes, j = seq_len(length(demes)))

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
    dplyr::filter(!duplicated(.)) %>%
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
  for (i in seq_len(nrow(gendist))) {
    gendist_arr[gendist$i[i], gendist$j[i], seq_len(nrow(gendist$data[[i]]))] <- unname(unlist(gendist$data[[i]]))
  }

  # normalize geodistances per user; NB have already removed self comparisons, so no 0s
  if (normalize_geodist) {
    mingeodist <- min(discdat$geodist)
    maxgeodist <- max(discdat$geodist)
    discdat <- discdat %>%
      dplyr::mutate(geodist = (geodist - mingeodist)/(maxgeodist - mingeodist))
  }
  # catch accidental bad M start if user is standardizing distances
  if (normalize_geodist & (start_params[names(start_params) == "m"] > 500) ) {
    warning("You have selected to normalize geographic distances, but your
            migration rate start parameter is large. Please consider placing it on a
            similar scale to your normalized geographic distances for stability.")
  }


  # put geo information into distance matrix
  geodist <- discdat %>%
    expand_pairwise(.) %>% # get all pairwise for full matrix
    dplyr::filter(!duplicated(.)) %>%
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
    # out
    unique(unlist(x))  # all same by unique
  }
  )

  # upper tri
  geodist_mat <- matrix(data = -1, nrow = length(locats), ncol = length(locats))
  for (i in seq_len(nrow(geodist))) {
    geodist_mat[geodist$i[i], geodist$j[i]] <- geodist$data[i]
  }
  diag(geodist_mat) <- 0

  #......................
  # out
  #......................
  ret <- list(
    demes = demes,
    n_Kpairmax = n_Kpairmax,
    keyi = keyi,
    gendist_arr = gendist_arr,
    geodist_mat = geodist_mat
  )
  return(ret)
}
