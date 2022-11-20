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

# Run the R scripts in the R/ folder with your custom functions:
tar_source("validation/Rvalidation")

#............................................................
# Target List
#...........................................................
list(
  tar_target(file1, "validation/mkdata/simulation_maestro.RDS", format = "file"),
  tar_target(sim_framework_df, readRDS(file1)),
  tar_target(file2, "validation/mkdata/simdata/locatcombo.rds", format = "file"),
  tar_target(locatcomb, readRDS(file2)),

  # run polySimIBD simulations and get disc data
  tar_target(discdat, swfsim_2_discdat_wrapper(sim_framework_df = sim_framework_df,
                                       reps = 100,
                                       dwnsmplnum = 5,
                                       locatcomb = locatcomb)),

  # run SA for each simulation type
  tar_target(SAoptimparms, get_SA_wrapper_start(discdat = discdat)),

  # run discent
  tar_target(results, get_discentwrapper(discdat = discdat,
                                         SAoptimparms = SAoptimparms)),

  # make report
  tar_render(report, "validation/validation_onering.Rmd")
)


