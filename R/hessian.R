#' @title Calculate Hessian Matrix from Loss Function
#' @inheritParams disc
#' @param mod class of DISCresult; output from `disc`
#' @description Calculate the Hessian, eigenvalues/vecotrs, and conditional number from the Loss Function using central finite differences
#' @return A list containing:
#'   \itemize{
#'     \item \code{Hessian}: Hessian matrix
#'     \item \code{Eigen}: Eigenvalues and eigenvectors from the Hessian
#'     \item \code{KappaH}: Conditional number from Hessian matrix
#'   }
#' @export

calculate_hessian_eigen <- function(mod, discdat, lambda) {

  #..............................................................
  # Assertions & Catches
  #..............................................................
  assert_dataframe(discdat)
  assert_single_numeric(lambda)
  assert_custom_class(output, "DISCresult")

  #............................................................
  # core
  #............................................................
  # numDeriv expects a function for which the first (vector) argument is used as a parameter vector.
  loss <- function(par_vec, discdat, n_demes, key, lambda) {

    fs <- par_vec[1:n_demes]
    m <- par_vec[n_demes+1]

    # cost[0] = 0.0;
    # for (int i = 0; i < (n_Demes-1); i++) {
    #   for (int j = i+1; j < n_Demes; j++) {
    #     double avg_fvec = (fvec[i] + fvec[j])/2;
    #     double exp_M = exp(-geodist_mat[i][j] / m);
    #     for (int k = 0; k < n_Kpairmax; k++){
    #       if (gendist_arr[i][j][k] != -1) {
    #         cost[0] += pow( (gendist_arr[i][j][k] -  avg_fvec * exp_M), 2);
    #       }
    #     }
    #   }
    # }
    # cost[0] += lambda * m * m; // explicit L2 regularization term at time 0
    cost <- 0
    # cost
    for(i in 1:(n_demes-1)) {
      for(j in (i+1):n_demes) {
        # take care of kij
        k_ij <- discdat %>%
          dplyr::filter((deme1 == key[i] & deme2 == key[j]) | (deme1 == key[j] & deme2 == key[i]))
        if(length(unique(k_ij$geodist)) != 1) {stop("Deme geodist mismatch")}
        # vals
        avg_fvec = (fs[i] + fs[j])/2
        exp_M = exp( unique(k_ij$geodist) / m)
        cost = cost + sum( (k_ij$gendist - avg_fvec * exp_M)^2)
        }
    }

  # Add L2 regularization term
  cost <- cost + lambda * m * m
  # out
  return(cost)
  }

  #......................
  # calculate hessian
  #......................
  par_vec <- c(mod$Final_Fis, mod$Final_m)
  names(par_vec) <- c(mod$deme_key$Deme, "m")
  H <- numDeriv::hessian(loss,
                    x = par_vec,
                    data = dat,
                    n_demes = length(mod$Final_Fis),
                    key = mod$deme_key$Deme,
                    lambda = lambda)

  #......................
  # calculate eigen and conditional number
  #......................
  eig <- eigen(H, symmetric = TRUE)
  kappaH <- kappa(H)


  #......................
  # out
  #......................
  out <- list(
    Hessian = H,
    Eigen = eig,
    KappaH = kappaH
  )

}
