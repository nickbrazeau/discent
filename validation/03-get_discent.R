## .................................................................................
## Purpose:
##
## Author: Nick Brazeau
##
## Date: 15 August, 2022
##
## Notes:
## .................................................................................
library(tidyverse)
library(discent)
library(furrr)
library(progressr)
set.seed(48)
setwd("validation") # setwd down one from package for validation work

# read in simulations
simdat <- readRDS("results/sim_data/sim_smpl_hosts_nonlinear_migration.rds") %>%
  dplyr::select(c("name", "discdat"))

#............................................................
#### PART 1: Make parameter search grid ####
#...........................................................
# make search grid to iterate over
searchgrid <- tidyr::expand_grid(f = seq(0.1, 0.9, 0.2),
                                 m = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2), # check distance
                                 m_learningrate = c(1e-14, 1e-12, 1e-10, 1e-8),
                                 f_learningrate = c(1e-9, 1e-7, 1e-5, 1e-3))
# tidy search grid into start params
get_start_params <- function(f, m) {
  our_start_params <- rep(f, 9)
  names(our_start_params) <- as.character(c(1,10,36,54,56,58,76,91,100))
  our_start_params <- c(our_start_params, "m" = f)
  return(our_start_params)
}
searchgrid <- searchgrid %>%
  dplyr::mutate(start_params = purrr::pmap(searchgrid[,c("f", "m")],
                                           get_start_params)) %>%
  dplyr::select(c("start_params", "f_learningrate", "m_learningrate"))

# now nest search grid for every discdat
simdat <- simdat %>%
  dplyr::mutate(sgrid = purrr::map(1:dplyr::n(),
                                   function(x){return(searchgrid)}))



#............................................................
#### PART 2: discent runs ####
#...........................................................
# wrapper function
disc_wrapper <- function(discdat, sgrid, pb) {
  # progress bar
  pb()
  # work
  discdatclean <- discdat %>%
    dplyr::filter(deme1 != deme2)
  # nested loop
  sgrid <- sgrid %>%
    dplyr::mutate(discret = purrr::pmap(sgrid[,c("start_params", "f_learningrate", "m_learningrate")],
                                        discent::deme_inbreeding_spcoef,
                                        discdat = discdatclean, steps = 1e4, momentum = 0.9,
                                        report_progress = F, return_verbose = F)
    )
  return(sgrid)
}



# run with future
future::plan(strategy = multisession, workers = future::availableCores()-1 )
progressr::with_progress({
  pb <- progressr::progressor(steps = nrow(simdat))
  simdat$discret <- furrr::future_pmap(simdat[,c("discdat", "sgrid")],
                                       disc_wrapper,
                                       pb = pb) })

#......................
# save out
#......................
dir.create("results/disc_results/", recursive = T)
saveRDS(simdat, "results/disc_results/simdat_searchgrid_disc_results.rds")
