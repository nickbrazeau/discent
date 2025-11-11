#' @title Identify Deme Inbreeding Spatial Coefficients in Continuous Space
#' @param discdat dataframe; The genetic-geographic data by deme (K). Must contain columns:
#'   \code{smpl1}, \code{smpl2}, \code{deme1}, \code{deme2}, \code{gendist}, \code{geodist}
#' @param start_params named double vector; vector of start parameters. Names must match deme names,
#'   plus one parameter named "m" for migration rate
#' @param lambda double; A quadratic L2 regularization parameter on "m" parameter: \eqn{\lambda m^2}. Default: 0.1
#' @param learningrate double; Learning rate (alpha) for gradient descent optimization. Default: 0.001
#' @param b1 double; Exponential decay rate for first moment estimate in Adam optimizer. Default: 0.9
#' @param b2 double; Exponential decay rate for second moment estimate in Adam optimizer. Default: 0.999
#' @param e double; Small constant for numerical stability in Adam optimizer. Default: 1e-8
#' @param steps integer; Number of optimization steps. Default: 1000
#' @param thin integer; Thinning interval for stored iterations (1 = store all). Default: 1
#' @param normalize_geodist logical; Whether to normalize geographic distances to [0,1] using
#'   min-max scaling: \eqn{X' = \frac{X - X_{min}}{X_{max} - X_{min}}}. Improves numerical
#'   stability but complicates interpretation of migration rate. Default: TRUE
#' @param report_progress logical; Whether to display progress bar during optimization. Default: TRUE
#' @param return_verbose logical; Whether to return full optimization trajectory (TRUE) or just
#'   final results (FALSE). Full trajectory can be memory intensive. Default: FALSE
#' @description This function estimates deme-specific inbreeding coefficients and a global migration
#'   rate from genetic and geographic distance data using an isolation-by-distance model. The model
#'   assumes that genetic similarity between samples decreases exponentially with geographic distance,
#'   modulated by deme-specific inbreeding coefficients.
#'
#'   The gradient calculation has been optimized to ensure all deme pairs contribute symmetrically
#'   to parameter updates, improving convergence and accuracy of estimates.
#' @details The input dataframe must have exactly these column names in order:
#'   \itemize{
#'     \item \code{smpl1}, \code{smpl2}: Sample identifiers
#'     \item \code{deme1}, \code{deme2}: Deme (location) identifiers
#'     \item \code{gendist}: Pairwise genetic distance [0,1]
#'     \item \code{geodist}: Pairwise geographic distance
#'   }
#' @details The \code{start_params} vector must contain:
#'   \itemize{
#'     \item One parameter per unique deme (named with deme identifiers)
#'     \item One parameter named "m" for the migration rate
#'     \item All F parameters must be in [0,1] (inbreeding coefficients)
#'   }
#' @details The model assumes: \eqn{E[r_{ij}] = \frac{F_i + F_j}{2} \exp(-d_{ij}/m)}
#'   where \eqn{r_{ij}} is genetic relatedness, \eqn{F_i} is deme i's inbreeding coefficient,
#'   \eqn{d_{ij}} is geographic distance, and \eqn{m} is the migration rate parameter.
#' @return A list of class "DISCresult" containing:
#'   \itemize{
#'     \item \code{Final_Fis}: Final inbreeding coefficient estimates
#'     \item \code{Final_m}: Final migration rate estimate
#'     \item \code{deme_key}: Mapping of deme names to indices
#'     \item \code{cost}: Final cost function value(s)
#'   }
#'   If \code{return_verbose = TRUE}, additional elements include:
#'   \itemize{
#'     \item \code{fi_run}: F parameter trajectory over iterations
#'     \item \code{m_run}: Migration parameter trajectory
#'     \item \code{fi_gradtraj}: F gradient trajectory
#'     \item \code{m_gradtraj}: Migration gradient trajectory
#'     \item \code{fi_1moment}, \code{fi_2moment}: F parameter Adam moments
#'     \item \code{m_1moment}, \code{m_2moment}: Migration parameter Adam moments
#'   }
#' @export

disc <- function(discdat,
                 start_params = NULL,
                 lambda = 0.1,
                 learningrate = 1e-3,
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
    message = "You have pairwise demes with differing geodistances measurements in your data input. Geodistances must be the same between demes (pairwise geodistances should be the same). Expected %s unique pairs but got %s."
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
  disclist <- discent:::wrangle_discentdat(discdat, normalize_geodist, start_params, locats)

  #..............................................................
  # run efficient C++ function
  #..............................................................

  args <- list(gendist = as.vector(disclist$gendist_arr),
               geodist = as.vector(disclist$geodist_mat),
               fvec = unname( disclist$start_params[!grepl("^m$", names(disclist$start_params))] ),
               n_Demes = length(disclist$demes),
               n_Kpairmax = disclist$n_Kpairmax,
               m = unname(disclist$start_params["m"]),
               learningrate = learningrate,
               lambda = lambda,
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
#' @param locats character vector; names of demes
#' @description Internal wrangling function that expands the DISC genetic data into an array and the geographic
#' distances into a matrix for efficient input into Cpp.
#' @noMd
#' @noRd

wrangle_discentdat <- function(discdat, normalize_geodist, start_params, locats) {
  # use efficient R functions to group pairs and wrangle data for faster C++ manipulation
  # get deme names and lift over sorted names for i and j (note, deme names may be anything, so cannot rely on user indexing)
  demes <- sort(unique(c(discdat$deme1, discdat$deme2)))
  keyi <- data.frame(deme1 = demes, i = seq_len(length(demes)))
  keyj <- data.frame(deme2 = demes, j = seq_len(length(demes)))

  # boundaries on genetic data for logit
  discdat <- discdat %>%
    dplyr::mutate(gendist = ifelse(gendist > 0.9999, 0.9999,
                                   ifelse(gendist < 0.0001, 0.0001,
                                          gendist))) # reasonable bounds for future logit

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
    geodist_mat = geodist_mat,
    start_params = start_params
  )
  return(ret)
}
