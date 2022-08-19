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
#### PART 1: Make spatial setup ####
#   plan will be a square always with some "genetic" barrier set of shapes
#   will use 'elevation' to parameterize migration rates
#...........................................................
latticemodel <- readRDS("results/sim_data/lattice_model.rds")
simdat <- tibble::tibble(name = c("ascent", "mtn", "oppcorner"),
                         gridmig = NA,
                         plotObj = NA)


#......................
# ascent
#......................
gridmig <- latticemodel %>%
  dplyr::mutate(migration = purrr::map_dbl(longnum,
                                           function(x){ x * 3 }))

# visualize to confirm
plot(raster::rasterFromXYZ(gridmig[,c(1,2,4)]))
raster::contour(raster::rasterFromXYZ(gridmig[,c(1,2,4)]),
                add = TRUE, drawlabels = FALSE, col = "#969696")

# store, same approx order of mag as distance for migration
simdat$gridmig[1] <- list( gridmig )



#......................
# migration climbing central mountain
#......................
gridmig <- latticemodel %>%
  dplyr::mutate(migration = purrr::map2_dbl(longnum, latnum, function(x, y){
    mvtnorm::dmvnorm(c(x, y),
                     mean = c(50, 50),
                     sigma = matrix(c(0.1, 1e-3, 1e-3, 0.1), ncol = 2),
                     log = T) }))

# visualize to confirm
plot(raster::rasterFromXYZ(gridmig[,c(1,2,4)]))
raster::contour(raster::rasterFromXYZ(gridmig[,c(1,2,4)]),
                add = TRUE, drawlabels = FALSE, col = "#969696")

# store, same approx order of mag as distance for migration
simdat$gridmig[2] <- list( gridmig %>%
                             dplyr::mutate(migration = migration/150) )

#......................
# opposite corners
#......................
gridmig <- latticemodel %>%
  dplyr::mutate(migration = dplyr::case_when(
    longnum <= 50 & latnum <= 50 ~ purrr::map2_dbl(longnum, latnum, function(x, y){
      -mvtnorm::dmvnorm(c(x, y),
                        mean = c(25, 25),
                        sigma = matrix(c(0.1, 1e-3, 1e-3, 0.1), ncol = 2),
                        log = T)}),

    longnum <= 50 & latnum > 50 ~ 0,
    longnum > 50 & latnum > 50 ~ purrr::map2_dbl(longnum, latnum, function(x, y){
      -mvtnorm::dmvnorm(c(x, y),
                        mean = c(75, 75),
                        sigma = matrix(c(0.1, 1e-3, 1e-3, 0.1), ncol = 2),
                        log = T)}),

    longnum > 50 & latnum <= 50 ~ 0

  ))

# make other direction hill
gridmig_min <- min(gridmig$migration, na.rm = T)
gridmig_max <- max(gridmig$migration, na.rm = T)

gridmig <- gridmig %>%
  dplyr::mutate(migration = dplyr::case_when(
    longnum <= 50 & latnum > 50 ~ migration,
    longnum <= 50 & latnum <= 50 ~ migration - gridmig_max,
    longnum > 50 & latnum > 50 ~ migration + gridmig_min,
    longnum > 50 & latnum <= 50 ~ migration
  ))

# visualize to confirm
plot(raster::rasterFromXYZ(gridmig[,c(1,2,4)]))
raster::contour(raster::rasterFromXYZ(gridmig[,c(1,2,4)]),
                add = TRUE, drawlabels = FALSE, col = "#969696")

# store, same approx order of mag as distance for migration
simdat$gridmig[3] <- list( gridmig %>%
                             dplyr::mutate(migration = migration/150) )
# clear this from environment
rm(gridmig)

#......................
# store nice plots
#......................
simdat <- simdat %>%
  dplyr::mutate(plotObj = purrr::map(gridmig, function(x){
    x %>%
      ggplot() +
      geom_tile(aes(x = longnum, y = latnum, fill = migration)) +
      geom_contour(aes(x = longnum, y = latnum, z = migration), color = "black") +
      scale_fill_viridis_c("Migration Topology") +
      coord_fixed() +
      theme_void()
  }))

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


# note already have euclidean now just need to add heights
simdat = simdat %>%
  dplyr::mutate(locats = purrr::map(dplyr::n(), function(x){return(locats)}))



get_dist_matrix  <- function(gridmig, locats, n = 9) {
  #......................
  # liftovers
  #......................
  # storage mat
  mat <- matrix(NA, nrow = n, ncol = n)
  # convert to raster
  rstr <- raster::rasterFromXYZ(gridmig[,c(1,2,4)])

  locats$height <- purrr::pmap_dbl(locats[c("longnum.x", "latnum.x", "longnum.y", "latnum.y")],
                                   function(longnum.x, latnum.x, longnum.y, latnum.y, rstr){
                                     # get height
                                     height1 <- raster::extract(rstr, matrix(c(longnum.x, latnum.x), nrow = 1))
                                     height2 <- raster::extract(rstr, matrix(c(longnum.y, latnum.y), nrow = 1))
                                     height <- sqrt((height2 - height1)^2)
                                     return(height)
                                   }, rstr = rstr)
  # bring together
  locats$distval <- locats$distval + locats$height

  #......................
  # fix values in matrix
  #......................
  # spread out values for matrix
  matlocat_dist <- locats %>%
    dplyr::select(c("deme1", "deme2", "distval")) %>%
    tidyr::pivot_wider(data = .,
                       names_from = "deme2",
                       values_from = "distval")
  # make sure columns are appropriately ordered
  matlocat_dist <- matlocat_dist[,c("1", "10", "36", "54", "56", "58", "76", "91", "100")]

  rownames(matlocat_dist) <- colnames(matlocat_dist)
  # convert to matrix
  matlocat_dist <- as.matrix(matlocat_dist)
  # scale
  matlocat_dist <- 1/matlocat_dist
  # now add in diagonal and offset for stay more often
  diag(matlocat_dist) <- 0
  diag(matlocat_dist) <- apply(matlocat_dist, 1, max, na.rm = T) + apply(matlocat_dist, 1, max, na.rm = T)/10
  # now re-make to probability matrix
  matlocat_dist <- matlocat_dist/rowSums(matlocat_dist)
  return(matlocat_dist)
}

#......................
# find distance matrices
#......................
simdat <- simdat %>%
  dplyr::mutate(distmat = purrr::pmap(simdat[,c("gridmig", "locats")], get_dist_matrix, n = 9))


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
saveRDS(simdat, "results/sim_data/sim_smpl_hosts_nonlinear_migration.rds")
