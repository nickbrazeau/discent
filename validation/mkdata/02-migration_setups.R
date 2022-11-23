## .................................................................................
## Purpose: Make various migration matrices to put DISCent through the wringer
##
## Author: Nick Brazeau
##
## Date: 17 November, 2022
##
## Notes:
## .................................................................................
library(tidyverse)

#............................................................
# Items I will make
#...........................................................
migmatdf <- tibble::tibble(modname = c("IsoByDist",
                                       "lattice",
                                       "torus"))


#............................................................
##### PART 0: Make a Square Matrix #####
#...........................................................
nCell <- 100
coords <- round(seq(1, nCell, by = 10))
latticemodel <- expand.grid(coords, coords)
plot(latticemodel)
colnames(latticemodel) <- c("longnum", "latnum")
demeNames <- 1:nrow(latticemodel)

latticemodel <- latticemodel %>%
  dplyr::mutate(deme = demeNames)
#......................
# cartesian distance matrix
#......................
# get combinations I need
locatcomb <- t(combn(sort(latticemodel$deme), 2)) %>%
  tibble::as_tibble(., .name_repair = "minimal") %>%
  magrittr::set_colnames(c("deme1", "deme2"))
# get selfs
selves <- tibble::tibble(deme1 = unique(c(locatcomb$deme1, locatcomb$deme2)),
                         deme2 = unique(c(locatcomb$deme1, locatcomb$deme2)))
locatcomb <- dplyr::bind_rows(locatcomb, selves)
# euc
locatcomb <- locatcomb %>%
  dplyr::mutate(geodist = furrr::future_pmap_dbl(locatcomb, function(deme1, deme2){
    # get long lat
    xy1 <- latticemodel[latticemodel$deme == deme1, c("longnum", "latnum")]
    xy2 <- latticemodel[latticemodel$deme == deme2, c("longnum", "latnum")]
    # euclidean distance
    euc <- dist(rbind(xy1, xy2), method = "euclidean")
    return(euc)})
  )
# make symmetrical
locatcomb_expand <- locatcomb
colnames(locatcomb_expand) <- c("deme2", "deme1", "geodist")
locatcomb <- dplyr::bind_rows(locatcomb, locatcomb_expand)
# remove duplicated diagonals
locatcomb <- locatcomb %>%
  dplyr::filter(!duplicated(locatcomb))
# now tidy up
latticemodel_x <- latticemodel %>%
  dplyr::rename(deme1 = deme,
                deme1longnum = longnum,
                deme1latnum = latnum)
latticemodel_y <- latticemodel %>%
  dplyr::rename(deme2 = deme,
                deme2longnum = longnum,
                deme2latnum = latnum)
locatcomb <- locatcomb %>%
  dplyr::left_join(., latticemodel_x, by = "deme1") %>%
  dplyr::left_join(., latticemodel_y, by = "deme2")

# expect this to be lower tri + upper tri + diagonals
goodegg:::assert_eq(nrow(locatcomb), choose(100,2)*2 + 100)
# save out for downstream
dir.create("validation/mkdata/simdata/", recursive = T)
saveRDS(locatcomb, "validation/mkdata/simdata/locatcombo.rds")


#......................
# util functions
#......................
make_wide_dist_mat <- function(locatcomb) {
  goodegg:::assert_eq(colnames(locatcomb),
                      c("deme1", "deme2", "geodist", "deme1longnum", "deme1latnum", "deme2longnum", "deme2latnum"))

  locatcombdistmat <- locatcomb %>%
    dplyr::select(!dplyr::contains("longnum")) %>%
    dplyr::select(!dplyr::contains("latnum")) %>%
    dplyr::arrange(deme1, deme2) %>%
    tidyr::pivot_wider(.,
                       names_from = deme2, values_from = geodist)

  locatcombdistmat <- as.matrix(locatcombdistmat)
  rownames(locatcombdistmat) <- locatcombdistmat[,1]
  locatcombdistmat <- locatcombdistmat[,2:ncol(locatcombdistmat)]
  return(locatcombdistmat)

}

#............................................................
##### PART 1: Isolation by Distance Matrix #####
#...........................................................
# assume rates are 1/dist
migmatdf$migmat <- NA
migmatdf$migmat[migmatdf$modname == "IsoByDist"] <- list(1/make_wide_dist_mat(locatcomb))

#............................................................
##### PART 2: Lattice Matrix #####
#...........................................................
nInds <- 100
nMov <- sqrt(nInds)
latticemigmat <- matrix(0, nInds, nInds)
# deme looks up = row + sqrt(nInds) in square
for (i in 1:nrow(latticemigmat)) {
  if (i + nMov < ncol(latticemigmat)) {
    latticemigmat[i, i + nMov] <- 1
  }
}

# deme looking down is just the reflection
latticemigmat[lower.tri(latticemigmat)] <- t(latticemigmat)[lower.tri(latticemigmat)]

# deme looks right  = row + 1 unless divisible by sqrt(nInds) (wall)
for (i in 1:nrow(latticemigmat)) {
  if (i %% nMov != 0) {
    latticemigmat[i, i + 1] <- 1
  }
}
# assume stay home 25% as much as move
diag(latticemigmat) <- 0.25
# store
migmatdf$migmat[migmatdf$modname == "lattice"] <- list(latticemigmat)


#............................................................
##### PART 3: Torus Matrix #####
#...........................................................
# adding to the latticemigmat a wrap around feature to be a torus
torusmigmat <- latticemigmat
wrparnd <- seq(nMov, nInds, by = nMov)
wrparnd_back <- seq(1, nInds, by = nMov)
# add in wrap around
torusmigmat[cbind(wrparnd, wrparnd_back)] <- 1

# store
migmatdf$migmat[migmatdf$modname == "torus"] <- list(torusmigmat)

#............................................................
# save out
#...........................................................
dir.create("validation/mkdata/simdata/", recursive = T)
saveRDS(migmatdf, "validation/mkdata/simdata/migmat_framework.RDS")
