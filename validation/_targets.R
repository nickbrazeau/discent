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

# Set target options:
tar_option_set(
  packages = c("tidyverse", "discent", "furrr"),
  format = "rds" )

# Run the R scripts in the R/ folder with your custom functions:
tar_source("validation/R")

#............................................................
# Target List
#...........................................................
list(
  tar_target(file1, "validation/mkdata/simdata/migration_mat_models.RDS", format = "file"),
  tar_target(migmat, readRDS(file1)),
  tar_target(file2, "validation/mkdata/simdata/locatcombo.rds", format = "file"),
  tar_target(locatcomb, readRDS(file2)),

  # run polySimIBD simulations and get disc data
  tar_target(discdat, get_disc_valid_data(migmat = migmat,
                                          reps = 100,
                                          locatcomb = locatcomb,
                                          dwnsmplnum = 5,
                                          Nesize = 25,
                                          nDemes = 7,
                                          demeNames = colnames(migmat$migmat[[1]])))
)


