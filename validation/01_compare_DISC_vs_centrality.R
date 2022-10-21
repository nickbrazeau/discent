## .................................................................................
## Purpose: Compare 100 DISC models with mean deme centrality calculated
## from network models
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
library(tidygraph)
library(ggraph)

#............................................................
# make simulations
#...........................................................
nmods <- 100
mods <- tibble::tibble(mod = 1:nmods) %>%
  dplyr::mutate(
    nDemes = sample(5:10, size = 100, replace = T),
    demesize = purrr::map(nDemes,
                          function(x){sapply(1:x, function(x) sample(2:10, size = 1))}),
    distmat = purrr::map(nDemes, discent::drawDistMat, rateDist = 1e-3),
    Ft = sample(seq(0,1,0.001), nmods))
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


mods <- mods %>%
  dplyr::mutate(Final_Fis = purrr::map(disc, "Final_Fis"),
                Final_m = purrr::map_dbl(disc, "Final_m"))

#............................................................
#### Compare w/ Centrality ####
#...........................................................
calc_centrality_by_smpl <- function(discdat) {
  # get node deme membership
  discdatexpand <- discdat
  colnames(discdat) <- c("smpl2", "smpl1", "deme2", "deme1", "gendist", "geodist")
  expand <- dplyr::bind_rows(discdat, discdatexpand)
  membership <- expand %>%
    dplyr::rename(deme = deme1,
                  name = smpl1) %>%
    dplyr::mutate(name = as.character(name)) %>%
    dplyr::select(c("name", "deme")) %>%
    dplyr::filter(!duplicated(.))

  # now get centrality
  centralret <- discdat %>%
    tidygraph::as_tbl_graph(., directed = F) %>%
    dplyr::mutate(host_importance = tidygraph::centrality_pagerank(weights = gendist)) %>%
    tidygraph::activate("nodes") %>%
    tibble::as_tibble() %>%
    dplyr::left_join(., membership, by = "name") %>%
    dplyr::rename(Deme = deme)
  return(centralret)
}


corr_cent_disc <- function(disc, centralret) {
  discret <- dplyr::bind_cols(disc$deme_key, disc = disc$Final_Fis)
  cordat <- dplyr::left_join(centralret, discret, by = "Deme")
  return(cor(cordat$host_importance, cordat$disc))
}

#......................
# bring together
#......................
mods <- mods %>%
  dplyr::mutate(centralret = purrr::map(discdat, calc_centrality_by_smpl),
                corr = purrr::map2_dbl(disc, centralret, corr_cent_disc))

mods %>%
  dplyr::filter(Final_m < 1) %>%
  ggplot() +
  geom_boxplot(aes(x = corr)) +
  coord_flip() +
  theme_linedraw()



mods %>%
  dplyr::filter(Final_m < 1) %>%
  ggplot() +
  geom_point(aes(x = Ft, y = corr)) +
  theme_linedraw()


summary(unlist(mods$Final_Fis))

