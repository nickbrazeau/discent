## .................................................................................
## Purpose: Using Targets to frame out my validation work for discent
##
## Author: Nick Brazeau
##
## Date: 13 November, 2022
## Notes:
##
## .................................................................................
library(tidyverse)
library(furrr)
library(future)
library(future.batchtools)
library(polySimIBD)
source("Rvalidation/polysim_wrappers.R")
source("Rvalidation/utils.R")

#............................................................
# read in and make polysim IBD dataframe
#...........................................................
reps <- 100
maestro <- readRDS("mkdata/simulation_maestro.RDS")
maestro <- lapply(1:reps, function(x){
  maestro <- maestro %>%
    dplyr::mutate(rep = x) %>%
    dplyr::select(c("modname", "rep", dplyr::everything()))
  return(maestro)}
) %>%
  dplyr::bind_rows()
locatcomb <- readRDS("mkdata/simdata/locatcombo.rds")


#............................................................
# run simulations on slurm
#...........................................................
ret <- maestro %>% dplyr::select(c("modname", "rep"))
plan(future.batchtools::batchtools_slurm, workers = availableCores(),
     template = "slurm_discent.tmpl")
ret$discdat <- furrr::future_pmap(maestro[,c("pos", "N", "m", "rho", "mean_coi", "tlim", "migr_mat")],
                                  swfsim_2_discdat_wrapper,
                                  dwnsmplnum = 5,
                                  locatcomb = locatcomb,
                                  .options = furrr_options(seed = TRUE))

dir.create("results")
saveRDS(ret, "results/discdat_from_polySimIBD_maestro.RDS")
