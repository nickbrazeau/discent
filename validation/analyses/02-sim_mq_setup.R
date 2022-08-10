## .................................................................................
## Purpose: Simulation from the spatial DTDLsWF to gauge "power" of discent model
##
## Notes:
## .................................................................................
library(tidyverse)
library(raster)
library(polySimIBD)
set.seed(44)


#............................................................
# QUESTION 1: Can we use Discent to do model comparisons/capture differences in geodistance
#...........................................................
#......................
# read in previous spatial setups
#......................
locatcomb <- readRDS("data/sim_data/locatcombo.rds")
locatcombexpand <- locatcomb
colnames(locatcombexpand) <- c("deme2", "deme1", "distval", "longnum.y", "latnum.y", "longnum.x", "latnum.x")
locatcomb <- rbind.data.frame(locatcomb, locatcombexpand) # now have all pairwise possibilities

#............................................................
# Make iso by distance matrix
#...........................................................
gridmigmat <- locatcomb %>%
  dplyr::filter(deme1 %in% c(1,10,56,91,100) &
                  deme2 %in% c(1,10,56,91,100)) %>%
  dplyr::filter(!duplicated(.)) %>%
  dplyr::select(c("deme1", "deme2", "distval")) %>%
  tidyr::pivot_wider(data = .,
                     names_from = "deme2",
                     values_from = "distval")
# make sure columns are appropriately ordered
gridmigmat <- gridmigmat[,c("1", "10", "56", "91", "100")]
# convert to matrix
gridmigmat <- as.matrix(gridmigmat)
# now add in diagonal and offset for stay more often -- 50% time at home
diag(gridmigmat) <- apply(gridmigmat, 1, max, na.rm = T) + apply(gridmigmat, 1, max, na.rm = T)/2
# now re-make to probability matrix
gridmigmat <- gridmigmat/rowSums(gridmigmat)


#............................................................
# run sWF simulator
#...........................................................
# magic numbers outside
Nesize <- 25
mscale <- 0.5
verity_coi2 <- readRDS("results/optim_coi_lambdas/optim_lambda.RDS")[2]

swf_sim_wrapper <- function(migmat, Nesize, mscale) {
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

  # from verity et al coi in the DRC: 2.23 (2.15– 2.31), going to use
  # a lambda of 2 as it is close
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
                                tlim =      tlim)
  return(swfsim)
}

#......................
# run simulations
#......................
swfsim <- swf_sim_wrapper(gridmigmat, Nesize = Nesize, mscale = mscale)


#............................................................
# Pairwise IBD realizations
#...........................................................
# downsample to 5 individuals per deme
dwnsmpl <- mapply(function(x,y){sample(x:y, size = 5, replace = F)},
       x = seq(1, 125, by = 25),
       y = seq(25, 125, by = 25),
       SIMPLIFY = F)

dwnsmpl <- sort(unlist(dwnsmpl))

# get pairwise IBD
comb_hosts_df <- t(combn(dwnsmpl, 2))
comb_hosts_list <- split(comb_hosts_df, 1:nrow(comb_hosts_df))
ibd <- purrr::map_dbl(comb_hosts_list, function(hosts, swf) {
  return(polySimIBD::get_pairwise_bv_ibd(swf = swf, host_index = hosts))
}, swf = swfsim)


#......................
# tidy up for out
#......................
comb_hosts_df <- tibble::as_tibble(comb_hosts_df) %>%
  magrittr::set_colnames(c("p1", "p2")) %>%
  dplyr::mutate(ibd = as.vector(unlist(ibd)))
# memberships
membership_x <- tibble::tibble(p1 = 1:125, deme1 = sort(rep(c(1,10,56,91,100), 25)))
membership_y <- tibble::tibble(p2 = 1:125, deme2 = sort(rep(c(1,10,56,91,100), 25)))
comb_hosts_df <- comb_hosts_df %>%
  dplyr::left_join(., membership_x) %>%
  dplyr::left_join(., membership_y) %>%
  dplyr::left_join(., locatcomb)



#............................................................
#### viz connections ####
#...........................................................
# see if connections are as expected
library(tidygraph)
library(ggraph)
adj_graph <- comb_hosts_df %>%
  tidygraph::as_tbl_graph(., directed = F)


comb_hosts_dfexpand <- comb_hosts_df
colnames(comb_hosts_dfexpand) <- c("p2", "p1", "ibd", "deme2", "deme1", "distval", "longnum.y", "latnum.y", "longum.x", "latnum.x")
expand <- dplyr::bind_rows(comb_hosts_df, comb_hosts_dfexpand)
membership <- expand %>%
  dplyr::rename(deme = deme1,
                name = p1) %>%
  dplyr::mutate(name = as.character(name)) %>%
  dplyr::select(c("name", "deme")) %>%
  dplyr::filter(!duplicated(.))


adj_graph %>%
  tidygraph::activate(., "nodes") %>%
  dplyr::left_join(., membership) %>%
  dplyr::mutate(community = as.factor(tidygraph::group_louvain(weights = ibd))) %>%
  tidygraph::activate("edges") %>%
  dplyr::filter(ibd > 0.25) %>%
  ggraph::ggraph(layout = 'kk') +
  ggraph::geom_edge_link(aes(width = ibd,
                             color = ibd)) +
  #ggraph::geom_node_point(aes(color = community), size = 3) +
  ggraph::geom_node_point(aes(color = factor(deme)), size = 3) +
  ggraph::scale_edge_width_continuous(range = c(0, 1), guide = "none") +
  scale_edge_color_viridis("IBD", values = c(0,1), option = "plasma") +
  scale_color_brewer(palette = "Set1") +
  ggraph::theme_graph() +
  theme(legend.position = "bottom")

#......................
# centrality
#......................
adj_graph %>%
  activate(nodes) %>%
  dplyr::left_join(., membership) %>%
  mutate(host_importance = tidygraph::centrality_pagerank(weights = ibd)) %>%
  tibble::as_tibble() %>%
  ggplot() +
  geom_boxplot(aes(x = factor(deme), y = host_importance, group = deme))

#............................................................
# save out
#...........................................................
saveRDS(swfsim, "data/sim_data/mq_swf_simulations.rds")
saveRDS(comb_hosts_df, "data/sim_data/sim_smpl_hosts_mq.rds")
