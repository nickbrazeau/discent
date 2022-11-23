## .................................................................................
## Purpose: Using Targets to frame out my validation work for discent
##
## Author: Nick Brazeau
##
## Date: 13 November, 2022
## Notes:
##
# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline # nolint
## .................................................................................
# Load packages required to define the pipeline:
library(targets)
library(tarchetypes)
library(tibble)
library(here)
setwd(here::here())

# Set target options:
tar_option_set(
  packages = c("tidyverse", "discent", "furrr"),
  format = "rds" )
options(clustermq.scheduler = "SLURM",
        clustermq.template = "/nas/longleaf/home/nfb/clustermq_template.slurm")

# Run the R scripts in the R/ folder with your custom functions:
tar_source("validation/Rvalidation")


#............................................................
# read in and make polysim IBD dataframe
#...........................................................
reps <- 100
maestro <- readRDS("validation/mkdata/simulation_maestro.RDS")
maestro <- lapply(1:reps, function(x){return(maestro)}) %>%
  dplyr::bind_rows()
locatcomb <- readRDS("validation/mkdata/simdata/locatcombo.rds")

#............................................................
# make search grid
#...........................................................
# look up tables
fstartsvec <- seq(from = 0.1, to = 0.9, by = 0.05)
mstartsvec <- 10^-seq(1,6)
flearnsvec <- 10^-seq(1,10)
mlearnsvec <- 10^-seq(5,15)
search_grid <- expand.grid(fstartsvec, mstartsvec, flearnsvec, mlearnsvec)
colnames(search_grid) <- c("fstart", "mstart", "f_learn", "m_learn")
# down sample
downsample <- 1e4
search_grid <- search_grid[sample(1:nrow(search_grid), size  = downsample, replace = F), ]
# template start
tempstart_params <- rep(0.1, 100)
names(tempstart_params) <- as.character(1:100)
tempstart_params <- c(tempstart_params, "m" = 1e-3)

# liftover to start param format
liftover_start_params <- function(fstart, mstart, start_param_template) {
  out <- start_param_template
  out[names(out) != "m"] <- fstart
  out[names(out) == "m"] <- mstart
  return(out)
}
search_grid <- search_grid %>%
  dplyr::mutate(start_params = purrr::map2(fstart, mstart, liftover_start_params,
                                           start_param_template = tempstart_params)) %>%
  dplyr::select(c("start_params", "f_learn", "m_learn"))

#............................................................
# WORK
#...........................................................
# run simulations
polysimtargets <- tar_map(value = maestro,
                          names = "sim",
                          tar_target(discdat,
                                     swfsim_2_discdat_wrapper(sim_framework_df = maestro,
                                                              dwnsmplnum = 5,
                                                              locatcomb = locatcomb)))
# bring together discdat initial
rettargets <- tar_combine(combined_disc,
                          polysimtargets,
                          command = dplyr::bind_rows(!!!.x))

# get start params setup
ibd <- tar_target(ibdstartdat, sub_maestro(combined_disc, lvl = "IsoByDist"))
lattice <- tar_target(latticestartdat, sub_maestro(combined_disc, lvl = "lattice"))
torus <- tar_target(torusstartdat, sub_maestro(combined_disc, lvl = "torus"))
nevary <- tar_target(nevarystartdat, sub_maestro(combined_disc, lvl = "NeVary"))

# run start params
ibdstart <- tar_map(value = search_grid,
                    names = "ibd",
                    tar_target(ibdstartdat,
                               get_GS_cost(searchgriddf = search_grid,
                                           discdat = ibd)))
latticestart <- tar_map(value = search_grid,
                        names = "lattice",
                        tar_target(latticestartdat,
                                   get_GS_cost(searchgriddf = search_grid,
                                               discdat = lattice)))
torusstart <- tar_map(value = search_grid,
                      names = "torus",
                      tar_target(torusstartdat,
                                 get_GS_cost(searchgriddf = search_grid,
                                             discdat = torus)))
nevarystart <- tar_map(value = search_grid,
                       names = "nevary",
                       tar_target(nevarystartdat,
                                  get_GS_cost(searchgriddf = search_grid,
                                              discdat = nevary)))


# bring together start param results
starttargets <- tar_combine(combined_start,
                          ibdstart, latticestart, torusstart, nevarystart,
                          command = dplyr::bind_rows(!!!.x))


# add in misspecified dist for cost
misspecified <- tar_target(fulldiscdat, add_misspec_dist(combined_disc))

# bring together
list(polysimtargets, rettargets,
     ibd, lattice, torus, nevary,
     ibdstart, latticestart, torusstart, nevarystart,
     misspecified)

