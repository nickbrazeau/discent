#' @title Find Optimal Gradient Parameters
#' @inheritParams deme_inbreeding_spcoef
#' @param discsteps integer; the number of "steps" in the main `DISCent`
#' algorithm" \link{deme_inbreeding_spcoef}
#' @param initstart_params named double vector; vector of start parameters.
#' @param initf_learningrate double; alpha parameter for how much each "step" is weighted in the gradient descent for inbreeding coefficients
#' @param initm_learningrate double; alpha parameter for how much each "step" is weighted in the gradient descent for the migration parameter
#' @param initTemp double; initial temperature for simulated annealer
#' @param annealstep integer; the number of "steps" for the simulated annealer
#' to consider for identifying optimal start and learning parameters
#' @param SLratio double; probability of searching for start parameters versus the
#' learning rate
#' @param FMSratio double; probability of searching for F start parameters versus the
#' M start parameters
#' @param FMLratio double; probability of searching for F learning rate versus the
#' M learning rate
#' @param demeSwitchSize integer; the number of demes to "move" in simulated annealing
#' proposal step
#'
#' @param fstartmin double; the minimum value for F inbreeding coefficient in the search grid
#' @param fstartmax double; the maximum value for F inbreeding coefficient in the search grid
#' @param fstartmin double; the number of values of F inbreeding coefficient to explore in the search grid
#'
#' @param mstartmin double; the minimum value for M global migration rate in the search grid
#' @param mstartmax double; the maximum value for M global migration rate in the search grid
#' @param mstartmin double; the number of values of M global migration rate to explore in the search grid
#'
#' @param flearnmin double; the minimum value for F learning rate in the search grid
#' @param flearnmax double; the maximum value for F learning rate in the search grid
#' @param flearnmin double; the number of values of F learning rate to explore in the search grid
#'
#' @param mlearnmin double; the minimum value for M learning rate in the search grid
#' @param mlearnmax double; the maximum value for M learning rate in the search grid
#' @param mlearnmin double; the number of values of M learning rate to explore in the search grid
#'
#'
#' @details Using simulated annealing, identify the best start parameters for
#' inbreeding F values and M global migration parameter as well as the best
#' learning rates. The simulated annealer assumes a discrete search grid
#' and temperate decreases according to the rate:
#' t = frac{t}{iter}
#' The cost in the simulated annealer is given by the final cost after _n_ iterations of
#' the main `DISCent` algorithm" \link{deme_inbreeding_spcoef}.
#' The search grid is defined by the user specified min, max, and steps of the
#' aforementioned parameters.
#' The various ratios of searching for the "X" parameters versus the "Y" parameters
#' helps to prioritize exploring one set of parameters versus the other. The default
#' value of 0.5 leads to equal probability of exploration.
#' @export

find_grad_params <- function(discdat,
                             initstart_params = c(),
                             initf_learningrate = 1e-5,
                             initm_learningrate = 1e-10,
                             momentum = 0.9,
                             discsteps = 1e3,

                             initTemp = 1,
                             annealsteps = 1e5,
                             SLratio = 0.5,
                             FMSratio = 0.5,
                             FMLratio = 0.5,
                             demeSwitchSize = 1,

                             fstartmin = 0.1,
                             fstartmax = 0.9,
                             fstartsteps = 100,

                             mstartmin = 1e-15,
                             mstartmax = 1e-3,
                             mstartsteps = 100,

                             flearnmin = 1e-15,
                             flearnmax = 1e-2,
                             flearnsteps = 100,

                             mlearnmin = 1e-15,
                             mlearnmax = 1e-2,
                             mlearnsteps = 100) {

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

  locats <- names(initstart_params)[!grepl("^m$", names(initstart_params))]
  if (!all(unique(c(discdat$deme1, discdat$deme2)) %in% locats)) {
    stop("You have cluster names in your discdat dataframe that are not included in your start parameters")
  }
  if (!any(grepl("^m$", names(initstart_params)))) {
    stop("A m start parameter must be provided (i.e. there must be a vector element name m in your start parameter vector)")
  }
  assert_dataframe(discdat)
  assert_vector(initstart_params)
  assert_length(initstart_params, n = (length(unique(c(discdat$deme1, discdat$deme2))) + 1),
                message = "Start params length not correct. You must specificy a start parameter
                           for each deme and the migration parameter, m")
  sapply(initstart_params[!grepl("^m$", names(initstart_params))], assert_bounded, left = 0, right = 1, inclusive_left = TRUE, inclusive_right = TRUE)
  assert_single_numeric(initf_learningrate)
  assert_single_numeric(initm_learningrate)
  assert_single_numeric(momentum)
  assert_single_int(discsteps)

  assert_single_int(annealsteps)
  assert_single_numeric(initTemp)
  assert_single_numeric(SLratio)
  assert_single_numeric(FMSratio)
  assert_single_numeric(FMLratio)
  assert_single_int(demeSwitchSize)

  assert_single_numeric(fstartmin)
  assert_single_numeric(fstartmax)
  assert_single_int(fstartsteps)

  assert_single_numeric(mstartmin)
  assert_single_numeric(mstartmax)
  assert_single_int(mstartsteps)

  assert_single_numeric(flearnmin)
  assert_single_numeric(flearnmax)
  assert_single_int(flearnsteps)

  assert_single_numeric(mlearnmin)
  assert_single_numeric(mlearnmax)
  assert_single_int(mlearnsteps)

  #............................................................
  # look up tables
  #...........................................................
  fstarts <- seq(from = fstartmin, to = fstartmax, length.out = fstartsteps)
  mstarts <- seq(from = mstartmin, to = mstartmax, length.out = mstartsteps)
  flearns <- seq(from = flearnmin, to = flearnmax, length.out = flearnsteps)
  mlearns <- seq(from = mlearnmin, to = mlearnmax, length.out = mlearnsteps)

  #............................................................
  # proposal function
  #...........................................................
  propose <- function(fstarts, mstarts, flearns, mlearns,
                      start_params, f_learn, m_learn,
                      SLratio, FMSratio, FMLratio, demeSwitchSize) {
    #......................
    # store
    #......................
    start_params <- start_params
    f_learn <- f_learn
    m_learn <- m_learn

    #......................
    # determine what is changing and make changes
    #......................
    if (rbinom(1, 1, SLratio)) { # change START params
      if (rbinom(1, 1, FMSratio)) { # change F start
        ds <- sample(names(start_params)[!grepl("^m$", names(start_params))],
                     size = demeSwitchSize,
                     replace = F)
        start_params[names(start_params) %in% ds] <- sample(fstarts,
                                                            size = demeSwitchSize,
                                                            replace = T)
      } else { # change M start
        start_params[names(start_params) == "^m$"] <- sample(mstarts,
                                                             size = 1)
      }

    } else { # change LEARNING rates
      if (rbinom(1, 1, FMLratio)) { # change f learning rate
        f_learn <- sample(flearns, 1)
      } else { # change m learning rate
        m_learn <- sample(mlearns, 1)
      }
    }

    #......................
    # send out
    #......................
    out <- list(
      start_params = start_params,
      f_learn = f_learn,
      m_learn = m_learn
    )
    return(out)
  }


  #............................................................
  # simulated annealer
  #...........................................................
  # STORAGE
  costrun <- rep(NA, annealsteps)
  # INIT
  Temp <- initTemp
  currprop <- list(start_params = initstart_params,
                   f_learn = initf_learningrate,
                   m_learn = initm_learningrate)
  currcost <- discent::deme_inbreeding_spcoef(discdat,
                                              start_params = currprop$start_params,
                                              f_learningrate = currprop$f_learn,
                                              m_learningrate = currprop$m_learn,
                                              momentum = momentum,
                                              steps = discsteps,
                                              report_progress = FALSE,
                                              return_verbose = FALSE)$cost[discsteps]

  # RUN
  for(i in 1:annealsteps) {
    # PROPOSE
    newprop <- propose(start_params = currprop$start_params,
                       f_learn = currprop$f_learn,
                       m_learn = currprop$m_learn, # NB currprop here is only for storage
                       fstarts, mstarts, flearns, mlearns, # R scoping takes care of this redundancy
                       SLratio, FMSratio, FMLratio, demeSwitchSize)
    # CALC COST
    newcost <- discent::deme_inbreeding_spcoef(discdat,
                                               start_params = newprop$start_params,
                                               f_learningrate = newprop$f_learn,
                                               m_learningrate = newprop$m_learn,
                                               momentum = momentum,
                                               steps = discsteps,
                                               report_progress = FALSE,
                                               return_verbose = FALSE)$cost[discsteps]
    # ACCEPT MOVE?
    u <- runif(1)
    p <- min(
      exp(-(newcost-currcost))/Temp,
      1)
    accept <- u <= p

    ## UPDATES
    # update current position and cost
    currprop <- if(accept){newprop}else{currprop}
    currcost <- if(accept){newcost}else{currcost}
    costrun[i] <- currcost
    # update temp
    Temp <- Temp/i
  }

  #............................................................
  # return
  #...........................................................
  out <- list(
    optimalParams = currprop,
    costrun = costrun
  )

  return(out)
}




