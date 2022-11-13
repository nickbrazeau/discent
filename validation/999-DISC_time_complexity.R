## .................................................................................
## Purpose: Calculate DISC Time Complexity
##
## Author: Nick Brazeau
##
## Date: 13 August, 2022
##
## Notes:
## .................................................................................
# set seed
set.seed(48)
library(tidyverse)
library(discent)
library(rbenchmark)

#TODO run w/ benchmark (loops better)

#............................................................
# make simulations
#...........................................................
mods <- tibble::tibble(mod = 1:100) %>%
  dplyr::mutate(
    nDemes = sample(5:10, size = 100, replace = T),
    demesize = purrr::map(nDemes,
                          function(x){sapply(1:x, function(x) sample(2:10, size = 1))}),
    distmat = purrr::map(nDemes, discent::drawDistMat, mDist = 100, vDist = 20),
    Ft = sample(seq(0,1,0.001), 100))
# draw simulations
mods$discdat <- purrr::pmap(mods[,c("demesize", "distmat", "Ft")],
                            discent::sim_IBDIBD, rate = 1e-3)
mods <- mods %>%
  dplyr::mutate(discdat = purrr::map(discdat, function(x){
    x %>%
      dplyr::filter(deme1 != deme2)
  }))
#............................................................
# get discent runs
#...........................................................
# start params
our_start_params <- rep(0.2, 2)
names(our_start_params) <- 1:2
our_start_params <- c(our_start_params, "m" = 1e-3)
mods <- mods %>%
  dplyr::mutate(start_params = purrr::map(nDemes, function(x){
                                            ret <- rep(0.2, x)
                                            names(ret) <- 1:x
                                            ret <- c(ret, "m" = 1e-3)
                                            return(ret) }))

mods$disc <- purrr::pmap(mods[,c("discdat", "start_params")],
                         discent::deme_inbreeding_spcoef,
                         momentum = 0.9, steps = 1e4,
                         report_progress = F, return_verbose = F)


