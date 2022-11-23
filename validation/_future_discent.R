## .................................................................................
## Purpose: Using future to find optimal start params
##
## Author: Nick Brazeau
##
## Date: 13 November, 2022
## Notes: unfortunately had to move to cluster for time
##
## .................................................................................
library(tidyverse)
library(furrr)
library(future)
library(future.batchtools)
library(polySimIBD)
source("Rvalidation/polysim_wrappers.R")
source("Rvalidation/discent_wrappers.R")
source("Rvalidation/utils.R")

#............................................................
# read in results
#...........................................................
beststarts <- readRDS("results/best_starts_for_discdat.RDS")
discdat <- readRDS("results/discdat_from_polySimIBD_maestro.RDS")
# bring home
fulldiscdat <- dplyr::left_join(discdat, beststarts, by = "modname")


#............................................................
# run discent
#...........................................................
ret <- fulldiscdat %>%
  dplyr::select(c("modname", "rep"))
plan(future.batchtools::batchtools_slurm, workers = availableCores(),
     template = "slurm_discent.tmpl")
ret$discret <- furrr::future_pmap(fulldiscdat[,c("discdat", "start_params", "f_learn", "m_learn")],
                                  get_discentwrapper,
                                  .options = furrr_options(seed = TRUE))

# out
dir.create("results", recursive = T)
saveRDS(ret, "results/discresults_for_discdat.RDS")


