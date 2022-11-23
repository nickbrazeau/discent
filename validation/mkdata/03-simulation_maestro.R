## .................................................................................
## Purpose: Make the simulation df needed to run polysimibd and convert to discdat
##
## Author: Nick Brazeau
##
## Date: 19 November, 2022
##
## Notes:
## .................................................................................
library(tidyverse)
simulation_df <- readRDS("validation/mkdata/simdata/migmat_framework.RDS")
simulation_df <- dplyr::bind_rows(simulation_df, simulation_df[1,])
simulation_df$modname[4] <- "NeVary"

#......................
# magic numbers
#......................
lambdaCOI <-  readRDS("validation/mkdata/simdata/optim_lambda.RDS")[2] #Default from verity et al for a coi 2 = 1.593624 given that the COI in the DRC: 2.23 (2.15– 2.31)
tlim <- 10 # assume IBD to 10 generations as before from verity et al
mscale <- 0.5 # assume that our mix of superinfection vs coinfection is split
nDemes <- nrow(simulation_df$migmat[[1]])

# Miles et al. 2016 (PMC5052046) & Taylor et al. 2019 (PMC6707449) gives us a recombination rate by 7.4e-7 M/bp
# Aimee gets this number by taking the inverse of Mile's estimate of the CO recombination rate of 13.5 kb/cM
rho <- 7.4e-7

# approximate average of Pf3d7 Chromosome Lengths
pflen <- 1.664e6
# assuming single chromosome for ease
# assuming 1e3 loci
pos <- sort(sample(pflen, 1e3))


#............................................................
# bring together
#...........................................................
# magic numbers common to all sims
simulation_df$rho <- list(rho)
simulation_df$pos <- list(pos)
simulation_df$m <- list(rep(mscale, nDemes))
simulation_df$tlim <- 10
simulation_df$mean_coi = list(rep(lambdaCOI, nDemes))
simulation_df$demeNames <- list(colnames(simulation_df$migmat[[1]]))

# varying params
Nbasic <- rep(25, nDemes)
Nvary <- rep(c(5,10,18,25,50,75,100,150,225,300), sqrt(nDemes))
simulation_df$N <- NA
simulation_df$N[1:3] <- list(Nbasic)
simulation_df$N[4] <- list(Nvary)

# rename
simulation_df <- simulation_df %>%
  dplyr::rename(migr_mat = migmat)

#............................................................
# save out
#...........................................................
saveRDS(simulation_df, "validation/mkdata/simulation_maestro.RDS")
