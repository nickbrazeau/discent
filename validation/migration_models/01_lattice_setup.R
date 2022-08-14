## .................................................................................
## Purpose: Make a lattice model to use for the polySimIBD modeling framework
##
## Notes:
## .................................................................................
library(tidyverse)
library(raster)
set.seed(48)

#............................................................
# Make a spatial setup to ask two question (01; 02 scripts)
#     Will use a Lattice Model
#...........................................................
#......................
# Lattice
#......................
nCell <- 100
coords <- round(seq(1, nCell, by = 10))
latticemodel <- expand.grid(coords, coords)
plot(latticemodel)
colnames(latticemodel) <- c("longnum", "latnum")

# store, same approx order of mag as distance for migration
latticemodel <- latticemodel %>%
  dplyr::mutate(deme = 1:dplyr::n())

#............................................................
# cartesian distance matrix
#...........................................................
#......................
# get combinations I need
#......................
locatcomb <- t(combn(sort(latticemodel$deme), 2)) %>%
  tibble::as_tibble(., .name_repair = "minimal") %>%
  magrittr::set_colnames(c("deme1", "deme2"))
# get selfs
selves <- tibble::tibble(deme1 = unique(c(locatcomb$deme1, locatcomb$deme2)),
                         deme2 = unique(c(locatcomb$deme1, locatcomb$deme2)))
locatcomb <- dplyr::bind_rows(locatcomb, selves)
#......................
# calculate for euclidean and height
#......................
# euc and height
locatcomb$distval <- furrr::future_pmap_dbl(locatcomb, function(deme1, deme2){
  # get long lat
  xy1 <- latticemodel[latticemodel$deme == deme1, c("longnum", "latnum")]
  xy2 <- latticemodel[latticemodel$deme == deme2, c("longnum", "latnum")]

  # euclidean distance
  euc <- dist(rbind(xy1, xy2), method = "euclidean")
  return(euc)}
  )


# now combine
latticemodel_x <- latticemodel %>%
  dplyr::rename(deme1 = deme)
latticemodel_y <- latticemodel %>%
  dplyr::rename(deme2 = deme)
locatcomb <- locatcomb %>%
  dplyr::left_join(., latticemodel_x, by = "deme1") %>%
  dplyr::left_join(., latticemodel_y, by = "deme2")

# save out
saveRDS(latticemodel, file = "data/sim_data/lattice_model.rds")
saveRDS(locatcomb, file = "data/sim_data/locatcombo.rds")
