#' Simulated Identity by Descent Data for DISC Analysis
#' @format A dataframe with 45 rows and 6 columns representing pairwise genetic and geographic distances:
#' \describe{
#'   \item{smpl1, smpl2}{Integer sample identifiers (1-10)}
#'   \item{deme1, deme2}{Factor deme (location) identifiers (1-3)}
#'   \item{gendist}{Numeric genetic distances in (0,1) based on simulated identity by descent}
#'   \item{geodist}{Numeric geographic distances (500, 1000) representing spatial separation}
#' }
#' @details This dataset contains pairwise comparisons between samples from 3 demes:
#'   \itemize{
#'     \item Deme 1-2 pairs: Geographic distance = 500
#'     \item Deme 1-3 pairs: Geographic distance = 1000
#'     \item Deme 2-3 pairs: Geographic distance = 1500
#'   }
#'   The data follows an isolation-by-distance model where genetic similarity decreases
#'   exponentially with geographic distance, modulated by deme-specific inbreeding coefficients.
#' @source Simulated toy dataset for testing and demonstration purposes. Generated using
#'   an exponential decay relationship between genetic relatedness and geographic distance.
#'   Not intended for real-world analysis - use for testing DISC functions and understanding
#'   expected input format.
#' @examples
#' # Load the dataset
#' data("IBD_simulation_data")
#'
#' # View structure
#' str(IBD_simulation_data)
#'
#' # Basic usage with disc function
#' \dontrun{
#' start_params <- c("1" = 0.2, "2" = 0.3, "3" = 0.4, "m" = 800)
#' result <- disc(IBD_simulation_data, start_params, steps = 100)
#' }
"IBD_simulation_data"
