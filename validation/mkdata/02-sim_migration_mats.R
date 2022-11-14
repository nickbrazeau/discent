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
latticemodel <- readRDS("validation/mkdata/simdata/lattice_model.rds")
latticemodel <- latticemodel[c(1,11,58,61,64,111,121),]
latticemodel <- latticemodel %>%
  dplyr::arrange(longnum, latnum)

# expand for join
locatcomb <- readRDS("validation/mkdata/simdata/locatcombo.rds")
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
# outs
#...........................................................
migmats <- tibble::tibble(
  modname = c("source", "sink", "radiate", "lovers"),
  migmat = list(spring, blackhole, radiate, lovers)
)

saveRDS(migmats, "validation/mkdata/simdata/migration_mat_models.RDS")

