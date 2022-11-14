
#' @title Run polySimIBD Wrapper
swf_sim_wrapper <- function(migmat) {
  #......................
  # magic numbers
  #......................
  # misc
  Nesize <- 25
  mscale <- 0.5
  verity_coi2 <- 1.593624

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


#' @title Convert swfsim to DISC data
#' @detials Pairwise IBD is from Verity et al 2020. Downsample to allow Ne
#' size to be large but not to overwhelm calculations and to simulate
#' real-life of subset sampling
#
get_swfsim_2_discdat <- function(swfsim, dwnsmplnum, locatcomb, Nesize, nDemes, demeNames) {
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
                                 deme1 = as.character(unlist(sapply(demeNames, function(x) rep(x,Nesize), simplify = F))))
  membership_y <- tibble::tibble(smpl2 = 1:(Nesize * nDemes),
                                 deme2 = as.character(unlist(sapply(demeNames, function(x) rep(x,Nesize), simplify = F))))
  # confirm compat
  locatcomb <- locatcomb %>%
    dplyr::mutate(deme1 = as.character(deme1),
                  deme2 = as.character(deme2))
  # bring together
  discdat <- comb_hosts_df %>%
    dplyr::left_join(., membership_x, by = "smpl1") %>%
    dplyr::left_join(., membership_y, by = "smpl2") %>%
    dplyr::left_join(., locatcomb, by = c("deme1", "deme2")) %>%
    dplyr::rename(geodist = distval) %>%
    dplyr::select(c("smpl1", "smpl2", "deme1", "deme2", "gendist", "geodist"))

  # out
  return(discdat)
}

#' @title Wrap Simulation to DISC data spectrum
#'

get_disc_valid_data <- function(migmat, reps,
                                dwnsmplnum, locatcomb,
                                Nesize, nDemes, demeNames) {
  # get full migmat
  migmatfull <- lapply(1:reps, function(x){return(migmat)}) %>%
    dplyr::bind_rows()
  # run polysimIBD and get IBD for discdat
  ret <- migmatfull %>%
    dplyr::mutate(swfsim = furrr::future_map(migmat, swf_sim_wrapper,
                                             seed = T),
                  discdat = furrr::future_map(swfsim, get_swfsim_2_discdat,
                                       dwnsmplnum, locatcomb,
                                       Nesize, nDemes, demeNames,
                                       seed = T)
    )
  return(ret)

}
