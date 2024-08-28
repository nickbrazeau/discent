#' @title Identify Deme Inbreeding Spatial Coefficients in Continuous Space with
#' Particle Swarm Meta-Optimization
#' @inheritParams deme_inbreeding_spcoef_vanilla
#' @param  fi_lowerinit double; The initial deme-inbreeding parameter lower-bound to parameterize
#' each swarm-particle (\eqn{a}), such that starting deme-inbreeding estimates are drawn from
#' a uniform distribution \eqn{f = U_{a,b} }
#' @param  fi_upperinit double; As above, the initial deme-inbreeding parameter upper-bound to parameterize
#' each swarm-particle (\eqn{b}).
#' @param  flearn_lowerinit double; The initial deme-inbreeding learning-rate lower-bound to parameterize
#' each swarm-particle's convergence to the \eqn{f} parameter (\eqn{a} draw from the unifrom \eqn{f = U_{a,b} }).
#' @param  flearn_upperinit double; As above, the initial deme-inbreeding learning-rate upper-bound to parameterize
#' each swarm-particle (\eqn{b}).
#' @param  mlearn_lowerinit double; similar to before, the initial migration learning-rate lower-bound to parameterize
#' each swarm-particle's convergence to the \eqn{m} parameter (\eqn{a} draw from the unifrom \eqn{m = U_{a,b} }).
#' @param  mlearn_upperinit double; As above, the initial migration learning-rate upper-bound to parameterize
#' each swarm-particle (\eqn{b}).
#' @param  c1 double; the "cognitive" coefficient from the PSO algorithm. Essentially, it
#' dictates how strongly the prior particle's positions should be weighted
#' against the entire swarm's historical positions in determining the next step of exploration.
#' @param  c2 double; the "social" coefficient from the PSO. Essentially determines how much weight
#' or influence other particles in the swarm exert on the current particle in determing the next step
#' of exploration.
#' @param  w double; the "inertia" coefficient from the PSO algorithm. Essentially, how
#' strongly the currently velocity (i.e. the direction the particle is headed) should be weighted relative to the prior particle's
#' and swarm's current positions (i.e. prior directions).
#' @param swarmsize integer; the number of particles in the swarm
#' @param  swarmmoves integer; the number of iterations or moves that the particles within the swarm are able
#' to explore before selecting the final particle for the "final" run (note, \emph{moves} is reserved for the swarm's actions, while \emph{steps}
#' is used to describe iterations in the gradient descent algorithm).
#' @param finalsteps integer; the number of "final" steps considered for the "final run" of the gradient descent
#' @param particlesteps integer; the number of steps that a particle takes in the vanilla gradient descent algorithm given its newly initialized start parameters in
#' order to calculate a cost for the new position. Essentially, we consider the vanilla gradient descent model at the current position for a short number of iterations to
#' estimate the "favorability of the positions current footing" or "traction" of the current position being considered.
#' @param report_sd_progress boolean; search chain
#' @param report_fd_progress boolean; final chain
#' @description The Particle Swarm Optimization (PSO) is a meta-optimization (meta-heuristic) approach that attempts to find optimal
#' start parameters for the user to avoid a grid-search approach as would be best practices for fine-tuning the gradient descent
#' approach.
#' @details Default values are based on ***
#' @references Clerc, M., and J. Kennedy. The Particle Swarm — Explosion, Stability, and Convergence in a Multidimensional Complex Space. IEEE Transactions on Evolutionary Computation 6, no. 1 (February 2002): 58–73. Y. H. Shi and R. C. Eberhart, “A modified particle swarm optimizer,” in Proceedings of the IEEE International Conferences on Evolutionary Computation, pp. 69–73, Anchorage, Alaska, USA, May 1998.
#' @export

deme_inbreeding_spcoef_pso <- function(discdat,
                                       m_lowerbound = 1e-10,
                                       m_upperbound = Inf,
                                       fi_lowerinit= 1e-3,
                                       fi_upperinit = 0.3,
                                       flearn_lowerinit = 1e-10,
                                       flearn_upperinit = 1e-2,
                                       mlearn_lowerinit = 1e-15,
                                       mlearn_upperinit = 1e-8,
                                       c1 = 2.0,
                                       c2 = 2.0,
                                       w = 0.73,
                                       b1 = 0.9,
                                       b2 = 0.999,
                                       e = 1e-8,
                                       finalsteps = 1e3,
                                       particlesteps = 1e2,
                                       swarmmoves = 50,
                                       swarmsize = 25,
                                       thin = 1,
                                       normalize_geodist = TRUE,
                                       report_sd_progress = TRUE,
                                       report_fd_progress = TRUE,
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

  assert_dataframe(discdat)
  assert_single_numeric(b1)
  assert_single_numeric(b2)
  assert_single_numeric(e)
  assert_single_numeric(c1)
  assert_single_numeric(c2)
  assert_single_numeric(w)
  assert_single_numeric(m_lowerbound)
  assert_single_numeric(m_upperbound)
  assert_gr(m_upperbound, m_lowerbound)
  assert_single_numeric(fi_lowerinit)
  assert_single_numeric(fi_upperinit)
  assert_gr(fi_upperinit, fi_lowerinit)
  assert_single_numeric(flearn_lowerinit)
  assert_single_numeric(flearn_upperinit)
  assert_gr(flearn_upperinit, flearn_lowerinit)
  assert_single_numeric(mlearn_lowerinit)
  assert_single_numeric(mlearn_upperinit)
  assert_gr(mlearn_upperinit, mlearn_lowerinit)
  assert_single_int(finalsteps)
  assert_single_int(thin)
  assert_greq(thin, 1, message = "Must be at least 1")
  assert_single_logical(report_sd_progress)
  assert_single_logical(report_fd_progress)
  assert_single_logical(normalize_geodist)

  # no missing
  if(sum(is.na(discdat)) != 0) {
    stop("discdat dataframe cannot have missing values")
  }

  # catch accidental bad F and M bound
  if ( any(round(c(fi_lowerinit, m_lowerbound), digits = 1e200) == 0) ) {
    warning("The Fi or M lower-bound is zero (or essentially zero), which will result in unstable behavior in the Gradient-Descent algorithm. Consider increasing the lower-bound limit")
  }

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
  gendist_arr <- array(data = -1, dim = c(length(demes), length(demes), n_Kpairmax))
  for (i in 1:nrow(gendist)) {
    gendist_arr[gendist$i[i], gendist$j[i], 1:nrow(gendist$data[[i]])] <- unname(unlist(gendist$data[[i]]))
  }

  # normalize geodistances per user; NB have already removed self comparisons, so no 0s
  if (normalize_geodist) {
    mingeodist <- min(discdat$geodist)
    maxgeodist <- max(discdat$geodist)
    discdat <- discdat %>%
      dplyr::mutate(geodist = (geodist - mingeodist)/(maxgeodist - mingeodist))
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
  geodist_mat <- matrix(data = -1, nrow = length(demes), ncol = length(demes))
  for (i in 1:nrow(geodist)) {
    geodist_mat[geodist$i[i], geodist$j[i]] <- geodist$data[i]
  }
  diag(geodist_mat) <- 0

  #..............................................................
  # run efficient C++ function
  #..............................................................

  args <- list(gendist = as.vector(gendist_arr),
               geodist = as.vector(geodist_mat),
               n_Demes = length(demes),
               n_Kpairmax = n_Kpairmax,
               m_lowerbound = m_lowerbound,
               m_upperbound = m_upperbound,
               fi_lowerinit = logit( fi_lowerinit ),
               fi_upperinit = logit( fi_upperinit ),
               flearn_lowerinit = flearn_lowerinit,
               flearn_upperinit = flearn_upperinit,
               mlearn_lowerinit = mlearn_lowerinit,
               mlearn_upperinit = mlearn_upperinit,
               b1 = b1,
               b2 = b2,
               e = e,
               c1 = c1,
               c2 = c2,
               w = w,
               swarmmoves = swarmmoves,
               swarmsize = swarmsize,
               particlesteps = particlesteps,
               steps = finalsteps, # the vanilla GD particle just calls this steps, as in the vanilla R file
               report_sd_progress = report_sd_progress,
               report_fd_progress = report_fd_progress,
               return_verbose = return_verbose
  )

  # run
  output_raw <- pso_deme_inbreeding_coef_cpp(args)


  # set up thinning
  thin_its <- seq(1, finalsteps, by = thin)
  thin_its <- unique(c(thin_its, finalsteps)) # always include last iteration

  # process output
  colnames(keyi) <- c("Deme", "key")
  if (return_verbose) {

    #......................
    # tidy up swarm
    #......................
    swarmtidy <- as.data.frame( matrix(NA, nrow = swarmsize * swarmmoves, ncol = 15))
    rowiter <- 1
    for (t in 1:swarmmoves) {
      for (i in 1:swarmsize) {
        swarmtidy[rowiter,1] <- t
        swarmtidy[rowiter,2] <- i
        swarmtidy[rowiter,3:15] <- unlist(output_raw$swarm[[t]][[i]])
        rowiter <- rowiter + 1
      }
    }
    colnames(swarmtidy) <- c("swarm_step_t", "particle_int",
                             "particle_Poscurr_Fstart", "particle_Poscurr_Mstart", "particle_Poscurr_Flearn", "particle_Poscurr_Mlearn",
                             "particle_Posbest_Fstart", "particle_Posbest_Mstart", "particle_Posbest_Flearn", "particle_Posbest_Mlearn", "particle_Posbest_Cost",
                             "particle_Veloccurr_Fstart", "particle_Veloccurr_Mstart", "particle_Veloccurr_Flearn", "particle_Veloccurr_Mlearn")
    #......................
    # rest of output
    #......................
    output <- list(
      deme_key = keyi,
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
      swarm = swarmtidy,
      global_best = output_raw$global_best
    )

  } else {
    output <- list(
      deme_key = keyi,
      cost = output_raw$cost[thin_its],
      m_run = output_raw$m_run[thin_its],
      fi_run = expit(do.call("rbind", output_raw$fi_run))[thin_its, ],
      Final_Fis = expit(output_raw$Final_Fis),
      Final_m = output_raw$Final_m)
  }

  # add S3 class structure
  attr(output, "class") <- "psoDISCresult"
  return(output)
}

