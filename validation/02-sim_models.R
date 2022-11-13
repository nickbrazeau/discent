## .................................................................................
## Purpose: Simulation from the spatial DTDLsWF to gauge how
##          the discent model behaves under different "expected"
##          versus realized IBD pairings
##
## Notes: Migration matrix are from - to format
## .................................................................................
library(tidyverse)
library(raster)
library(polySimIBD)
set.seed(48)

#............................................................
# lattice model and location data
#...........................................................
latticemodel <- readRDS("results/sim_data/lattice_model.rds")
latticemodel <- latticemodel[c(1,11,58,61,64,111,121),]
latticemodel <- latticemodel %>%
  dplyr::arrange(longnum, latnum)

# expand for join
locatcomb <- readRDS("results/sim_data/locatcombo.rds")
locatcombexpand <- locatcomb
colnames(locatcombexpand) <- c("deme2", "deme1", "distval", "longnum.y", "latnum.y", "longnum.x", "latnum.x")
locatcomb <- dplyr::bind_rows(locatcomb, locatcombexpand)
locatcomb <- locatcomb %>%
  dplyr::filter(deme1 %in% latticemodel$deme &
                  deme2 %in% latticemodel$deme) %>%
  dplyr::filter(!duplicated(.)) %>%  #NB, diagonal on distmatrix gets duplicated
  dplyr::mutate(deme1 = as.character(deme1),
                deme2 = as.character(deme2))
#............................................................
#### PART 1: Make distance/migration setup ####
#   same plan of square with 9 samples
#   going to manually set up migration matrix for sink-source dynamics
#...........................................................
ogdistmat <- matrix(0, nrow = 7, ncol = 7)
rownames(ogdistmat) <- colnames(ogdistmat) <- latticemodel$deme

#......................
# spring
#......................
spring <- ogdistmat
spring[4,] <- 2 # middle deme has connections OUT only
diag(spring) <- 1 # if not out, stat in home deme
# liftover to rates
spring <- spring/sum(spring)

#......................
# blackhole
#......................
blackhole <- ogdistmat
blackhole[,4] <- 2 # middle deme has connections IN only
diag(blackhole) <- 1 # if not out, stat in home deme
# liftover to rates
blackhole <- blackhole/sum(blackhole)


#......................
# radiate
#......................
radiate <- ogdistmat
# central radiate
radiate[4,3] <- 2
radiate[4,5] <- 2
# left radiate
radiate[3,1] <- 2
radiate[3,2] <- 2
# right radiate
radiate[5,6] <- 2
radiate[5,7] <- 2
diag(radiate) <- 1 # if not out, stat in home deme
# liftover to rates
radiate <- radiate/sum(radiate)

#......................
# all-in-lover
#......................
lovers <- ogdistmat
diag(lovers) <- 1 # only connect w/ self
# liftover to rates
lovers <- lovers/sum(lovers)

#............................................................
# get locat distance for geodist needed for disc
#...........................................................


#............................................................
#### PART 2: Run Simulations ####
#...........................................................
simdat <- tibble::tibble(
  modname = c("spring", "blackhole", "radiate", "lovers"),
  migmat = list(spring, blackhole, radiate, lovers)
)

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
  nDemes <- nrow(migmat)
  #......................
  # run structured WF
  #......................
  swfsim <- polySimIBD::sim_swf(pos =       pos,
                                migr_mat = migmat,
                                N =         rep(Nesize, nDemes),
                                m =         rep(mscale, nDemes),
                                rho =       rho,
                                mean_coi =  rep(verity_coi2, nDemes),
                                tlim =      10)
  return(swfsim)
}


# get simulations
simdat <- simdat %>%
  dplyr::mutate(swfsim = purrr::map(migmat, swf_sim_wrapper))


#............................................................
#### PART 3: Pairwise IBD realizations ####
#...........................................................
# ibd wrapper
get_ibd_wrapper <- function(swfsim, dwnsmplnum, locatcomb, Nesize, nDemes, demeNames) {
  # downsample to "N" individuals per deme
  dwnsmpl <- mapply(function(x,y){sample(x:y, size = dwnsmplnum, replace = F)},
                    x = seq(1, Nesize * nDemes, by = Nesize),
                    y = seq(Nesize, Nesize * nDemes, by = Nesize),
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
  membership_x <- tibble::tibble(smpl1 = 1:(Nesize * nDemes),
                                 deme1 = unlist(sapply(demeNames, function(x) rep(x,Nesize), simplify = F)))
  membership_y <- tibble::tibble(smpl2 = 1:(Nesize * nDemes),
                                 deme2 = unlist(sapply(demeNames, function(x) rep(x,Nesize), simplify = F)))
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
                                     dwnsmplnum = 5,
                                     Nesize = 25,
                                     nDemes = 7,
                                     demeNames = colnames(ogdistmat)))



#......................
# save out
#......................
dir.create("results/sim_data/", recursive = T)
saveRDS(simdat, "results/sim_data/sim_smpl_hosts_sink_source.rds")
