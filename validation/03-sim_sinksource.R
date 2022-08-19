## .................................................................................
## Purpose: Simulation from the spatial DTDLsWF to gauge how
##          the discent model behaves under different "expected"
##          versus realized IBD pairings
##
## Notes:
## .................................................................................
library(tidyverse)
library(raster)
library(polySimIBD)
set.seed(48)
setwd("validation") # setwd down one from package for validation work
#............................................................
#### PART 1: Make distance/migration setup ####
#   same plan of square with 9 samples
#   going to manually set up migration matrix for sink-source dynamics
#...........................................................
ogdistmat <- matrix(NA, nrow = 9, ncol = 9)
rownames(ogdistmat) <- colnames(ogdistmat) <- as.character(c(1,10,36,54,56,58,76,91,100))

#......................
# spring
#......................
springmat <- ogdistmat
springmat[is.na(springmat)] <- 0 # init
springmat[5,] <- 2 # deme 56 has connections out only
diag(springmat) <- 3 # slight preference for staying in home deme than moving


#......................
# blackhole
#......................
blackhole <- ogdistmat
blackhole[is.na(blackhole)] <- 0 # init
blackhole[,5] <- 2 # deme 56 has connections in only
diag(blackhole) <- 3 # slight preference for staying in home deme than moving

#......................
# radiate
#......................
radiate <- ogdistmat
radiate[is.na(radiate)] <- 0 # init
#TODO
diag(radiate) <- 3 # slight preference for staying in home deme than moving



# TODO NFB here



#......................
locatcomb <- readRDS("results/sim_data/locatcombo.rds")
locatcombexpand <- locatcomb
colnames(locatcombexpand) <- c("deme2", "deme1", "distval", "longnum.y", "latnum.y", "longnum.x", "latnum.x")
locatcomb <- rbind.data.frame(locatcomb, locatcombexpand) %>%  # now have all pairwise possibilities
  dplyr::filter(!duplicated(.)) # to remove self comparison duplications
# downsize demes
locats <- locatcomb %>%
  dplyr::filter(deme1 %in% c(1,10,36,54,56,58,76,91,100) &
                  deme2 %in% c(1,10,36,54,56,58,76,91,100)) %>%
  dplyr::filter(!duplicated(.))

#......................
# Location Sampling and migration rate calculations
#   sample lattice locations
#......................
locatcomb <- readRDS("results/sim_data/locatcombo.rds")
locatcombexpand <- locatcomb
colnames(locatcombexpand) <- c("deme2", "deme1", "distval", "longnum.y", "latnum.y", "longnum.x", "latnum.x")
locatcomb <- rbind.data.frame(locatcomb, locatcombexpand) %>%  # now have all pairwise possibilities
  dplyr::filter(!duplicated(.)) # to remove self comparison duplications
# downsize demes
locats <- locatcomb %>%
  dplyr::filter(deme1 %in% c(1,10,36,54,56,58,76,91,100) &
                  deme2 %in% c(1,10,36,54,56,58,76,91,100)) %>%
  dplyr::filter(!duplicated(.))


#### JUST SUM THE DISTANCE MATRIX AND SINK SOURCE MATRIC



#............................................................
#### PART 2: Run Simulations ####
# Will run sWF simulator 100 times for each
# migration scenario
#...........................................................
nreps <- 100
simdat <- lapply(1:nreps, function(x){return(simdat)}) %>%
  dplyr::bind_rows() %>%
  dplyr::arrange(name) %>%
  dplyr::group_by(name) %>%
  dplyr::mutate(name = paste0(name, 1:dplyr::n())) %>%
  dplyr::ungroup()

# magic numbers for simulator
Nesize <- 25
mscale <- 0.5
verity_coi2 <- readRDS("results/optim_coi_lambdas/optim_lambda.RDS")[2]

swf_sim_wrapper <- function(migmat) {
  #......................
  # magic numbers
  #......................
  # Miles et al. 2016 (PMC5052046) & Taylor et al. 2019 (PMC6707449) gives us a recombination rate by 7.4e-7 M/bp
  # Aimee gets this number by taking the inverse of Mile's estimate of the CO recombination rate of 13.5 kb/cM
  rho <- 7.4e-7
  # going to assume we can only detect things 10 generations ago
  tlim <- 10

  # approximate average of Pf3d7 Chromosome Lengths
  pflen <- 1.664e6
  # assuming single chromosome for ease
  # assuming 1e3 loci
  pos <- sort(sample(1.664e6, 1e3))

  # from verity et al coi in the DRC: 2.23 (2.15– 2.31)
  # assuming deme size of 10 for ease
  # tlim at 10 generations as before from verity et al

  #......................
  # run structured WF
  #......................
  swfsim <- polySimIBD::sim_swf(pos =       pos,
                                migr_dist_mat = migmat,
                                N =         rep(Nesize, nrow(migmat)),
                                m =         rep(mscale, nrow(migmat)),
                                rho =       rho,
                                mean_coi =  rep(verity_coi2, nrow(migmat)),
                                tlim =      10)
  return(swfsim)
}


# get simulations
simdat <- simdat %>%
  dplyr::mutate(swfsim = purrr::map(distmat, swf_sim_wrapper))


#............................................................
#### PART 3: Pairwise IBD realizations ####
#...........................................................

# ibd wrapper
get_ibd_wrapper <- function(swfsim, dwnsmplnum, locatcomb) {
  # downsample to "N" individuals per deme
  dwnsmpl <- mapply(function(x,y){sample(x:y, size = dwnsmplnum, replace = F)},
                    x = seq(1, 225, by = 25),
                    y = seq(25, 225, by = 25),
                    SIMPLIFY = F)
  dwnsmpl <- sort(unlist(dwnsmpl))
  # get combinations
  comb_hosts_df <- t(combn(dwnsmpl, 2))
  comb_hosts_list <- split(comb_hosts_df, 1:nrow(comb_hosts_df))
  # get pairwise IBD
  ibd <- purrr::map_dbl(comb_hosts_list, function(hosts, swf) {
    return(polySimIBD::get_pairwise_bv_ibd(swf = swf, host_index = hosts))
  }, swf = swfsim)

  # tidy up for out
  comb_hosts_df <- tibble::as_tibble(comb_hosts_df) %>%
    magrittr::set_colnames(c("smpl1", "smpl2")) %>%
    dplyr::mutate(gendist = as.vector(unlist(ibd)))
  # memberships
  membership_x <- tibble::tibble(smpl1 = 1:225, deme1 = sort(rep(c(1,10,36,54,56,58,76,91,100), 25)))
  membership_y <- tibble::tibble(smpl2 = 1:225, deme2 = sort(rep(c(1,10,36,54,56,58,76,91,100), 25)))
  discdat <- comb_hosts_df %>%
    dplyr::left_join(., membership_x, by = "smpl1") %>%
    dplyr::left_join(., membership_y, by = "smpl2") %>%
    dplyr::left_join(., locatcomb, by = c("deme1", "deme2")) %>%
    dplyr::rename(geodist = distval) %>%
    dplyr::select(c("smpl1", "smpl2", "deme1", "deme2", "gendist", "geodist"))

  # out
  return(discdat)
}

# get pairwise IBD
simdat <- simdat %>%
  dplyr::mutate(discdat = purrr::map(swfsim, get_ibd_wrapper,
                                     locatcomb = locatcomb,
                                     dwnsmplnum = 5))



#......................
# save out
#......................
dir.create("results/sim_data/", recursive = T)
saveRDS(simdat, "results/sim_data/sim_smpl_hosts_sink_source.rds")
