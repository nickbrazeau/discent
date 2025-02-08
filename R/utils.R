#' @title Calculate Hessian Matrix from Cost Chain
#' @param
#' @description
#' @details
#' @returns
#' @export

# Calculate_Hessian_of_Cost <- function(costchain, burn) {
#   # checks
#   # goodegg::assert_vector(costchain)
#   # goodegg::assert_int(burn)
#   # goodegg::assert_le(burn, lenght(costchain), message = "Burn-in exceeds the length of the cost chain")
#
#
#   # core
#   costchain_burned <- costchain[(burn+1):length(costchain)] # burn in
#   nlcost <- length(costchain_burned) # store length
#   H <- matrix(NA, nrow = nlcost, ncol = nlcost) # blank hessian matrix
#   # Compute Hessian using finite differences method
#   for (i in 1:nlcost) {
#     # Finite Differences for On-Diagonal Second Derivative (not mixed)
#     # [f(x+h) - 2f(x) + f(x-h)]/h^2 https://en.wikipedia.org/wiki/Finite_difference#Relation_with_derivatives
#     # NB, little h=1, for 1 step so is dropped
#     if (i == 1) {
#       H[i, i] <- costchain_burned[i+2] - 2 * costchain_burned[i+1] + costchain_burned[i] # forward lead at start
#     } else if (i == nlcost) {
#       H[i, i] <- costchain_burned[i] - 2 * costchain_burned[i-1] + costchain_burned[i-2] # back lag at start
#     } else {
#       H[i, i] <- costchain_burned[i+1] - 2 * costchain_burned[i] + costchain_burned[i-1]
#     }
#     # there are NO mixed partial derivative, off-diagonal since we have one dimensional data
#     # [f(x_i+h, x_j+h) - f(x_i+h, x_j) - f(x_i, x_j+h) + f(x_i, x_j)] / 4h^2 https://math.stackexchange.com/questions/888259/can-someone-explain-in-general-what-a-central-difference-formula-is-and-what-it
#   }
#   # out
#   return(H)
# }


#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL

#------------------------------------------------
# update progress bar.
# pb_list = list of progress bar objects
# name = name of this progress bar
# i = new value of bar
# max_i = max value of bar (close when reach this value)
#' @noRd
update_progress <- function(pb_list, name, i, max_i) {
  utils::setTxtProgressBar(pb_list[[name]], i)
  if (i == max_i) {
    close(pb_list[[name]])
  }
}

#------------------------------------------------
#' @title logit transformation
#' @param p numeric vecotr
#' @description Standard expit formula
#' @noMd
#' @noRd
# no export because simple
logit <- function(p){
  # out
   log(p/(1-p))
}

#------------------------------------------------
#' @title expit transformation
#' @param p numeric vecotr
#' @description Standard expit formula
#' @noMd
#' @noRd
# no export because simple
expit <- function(p){
  # out
  1/(1+exp(-p))
}

#------------------------------------------------
#' @title Simple function for expanding a pairwise matrix
#' @param discdf dataframe; that is specific to the DISCent framework. Must have colnames of: "smpl1", "smpl2", "deme1", "deme2", "gendist", "geodist" in correct order
#' @description This is an internal function that largely lacks generalizability. Do not recommend external use.
#' @noMd
#' @noRd

expand_pairwise <- function(discdf){
  #......................
  # assertions
  #......................
  assert_dataframe(discdf)
  assert_eq(colnames(discdf), c("smpl1", "smpl2", "deme1", "deme2", "gendist", "geodist"))

  #......................
  # expansion
  #......................
  discdfexpand <- discdf
  colnames(discdfexpand) <- c("smpl2", "smpl1", "deme2", "deme1", "gendist", "geodist")
  discdfexpand <- rbind.data.frame(discdf, discdfexpand) # now have all pairwise possibilities
  discdfexpand <- discdfexpand[!duplicated(discdfexpand), ]
  # out
  discdfexpand
}


#............................................................
# DISCresult S3 Class Overloading
#...........................................................
#' @title Check if DISCresult S3 Class
#' @description Overload is: function for determining if object is of class DISCresult
#' @param x DISC result from deme_inbreeding_spcoef function
#' @noMd
#' @export
is.DISCresult <- function(x) {
  inherits(x, "is.DISCresult")
}

#' @title print DISCresult S3 Class
#' @description overload print() function to print summary only
#' @inheritParams is.DISCresult
#' @param ... further arguments passed to or from other methods.
#' @noMd
#' @export
print.DISCresult <- function(x, ...) {

  # print summary only
  cat(crayon::red("Final DISC Range:"),  paste(round(min(x$Final_Fis),2), round(max(x$Final_Fis),2), sep = " - "), "\n")
  cat(crayon::blue("Final Migration Rate:"), round(x$Final_m, 3), "\n")

}

#' @title Summary of DISCresult S3 Class
#' @description overload summary() function.
#' @param object DISCresult Simulation
#' @param ... further arguments passed to or from other methods.
#' @noMd
#' @export
summary.DISCresult <- function(object, ...) {
  # send summary only
  tidyout.DISCresult(object)
}


#............................................................
# organizing disc output
#...........................................................
#' @title Tidy Out Sim Method
#' @description Method assignment
#' @inheritParams is.DISCresult
#' @noMd
#' @export
tidyout <- function(x) {
  UseMethod("tidyout")
}

#' @title Tidy Out Sim
#' @description Function for taking output of SIR NE and lifting it over
#' @inheritParams is.DISCresult
#' @noMd
#' @export

tidyout.DISCresult <- function(x) {
  #......................
  # clean up
  #......................
  fis <- dplyr::bind_cols(x$deme_key, x$Final_Fis) %>%
    dplyr::select(-key) %>%
    magrittr::set_colnames(c("Deme", "DISC"))

  # out
  list(Final_Fis = fis,
       Final_M = x$Final_m)
}
