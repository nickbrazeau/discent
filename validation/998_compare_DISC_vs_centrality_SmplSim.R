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


#............................................................
#### viz connections ####
#...........................................................
# see if connections are as expected

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
