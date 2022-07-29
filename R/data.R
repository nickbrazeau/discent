#' Malaria Gene Gene Flow up a Mountain
#'
#' A dataset that was simulated using the \code{https://github.com/nickbrazeau/polySimIBD} package:
#' \url{https://github.com/nickbrazeau/polySimIBD}, which is a "cartoon-ish" version for malaria
#' transmission that relies on a discrete-time discrete-loci structured Wright-Fisher model
#'
#' @format A list with the following items:
#' \describe{
#'   \item{locatdat}{Dataframe: Mock latitudes and longitudes for each deme that was simualted on a 30x30 square}
#'   \item{gridmig}{Dataframe: The migration matrix that defines gene flow between demes}
#'   \item{migplot}{ggplot: A rasterized plot of the migration matrix}
#'   \item{discent_dat}{Dataframe: The discent data in the correct format}
#' }

"discent_sim_dat"
