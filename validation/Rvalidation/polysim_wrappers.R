# expand out reps
expand_sim_framework_df <- function(sim_framework_df, reps) {
  out <- dplyr::bind_rows(lapply(1:reps, function(x){return(sim_framework_df)}))
  return(out)
}


swfsim_2_discdat_wrapper <- function(pos, N, m,
                                  rho, mean_coi, tlim,
                                  migr_mat,
                                  dwnsmplnum,
                                  demeNames, locatcomb) {
    # run swfsim
    swfsim <- polySimIBD::sim_swf(pos, N, m, rho, mean_coi, tlim, migr_mat)
    # # get ibd dat
    # ibddat <- get_swfsim_2_ibd(swfsim, N, dwnsmplnum)
    # # tidy to discent
    # ret <- IBDdat_tidy_out_2_discdat(N = N, demeNames = demeNames,
    #                                  locatcomb = locatcomb,
    #                                  IBDdat = ibddat)
    # # out
    # ret <- ret %>%
    #   dplyr::select(c("modname", "discdat"))
    # return(ret)
    return(0)
}




#' @title Convert swfsim to IBD
#' @detials Pairwise IBD is from Verity et al 2020. Downsample to allow Ne
#' size to be large but not to overwhelm calculations and to simulate
#' real-life of subset sampling
#
get_swfsim_2_ibd <- function(swfsim, N, dwnsmplnum){

    # get start and end ind counts
    inds <- lapply(N, function(x){seq(1, x, by = 1)}) # list of inds by deme
    end <- cumsum(sapply(inds, max))
    start <- end + 1
    start <- c(1, start[1:(length(start)-1)])
    # downsample to "N" individuals per deme
    dwnsmpl <- mapply(function(x,y){sample(x:y, size = dwnsmplnum, replace = F)},
                      x = start, y = end, SIMPLIFY = F)
    dwnsmpl <- sort(unlist(dwnsmpl))
    # get combinations
    comb_hosts_df <- t(combn(dwnsmpl, 2))
    comb_hosts_list <- split(comb_hosts_df, 1:nrow(comb_hosts_df))
    # get pairwise IBD
    ibd <- purrr::map_dbl(comb_hosts_list, function(hosts, swf) {
      return(polySimIBD::get_pairwise_bv_ibd(swf = swf, host_index = hosts))
    }, swf = swfsim)
    # tidy up for out
    ret <- tibble::as_tibble(comb_hosts_df) %>%
      magrittr::set_colnames(c("smpl1", "smpl2")) %>%
      dplyr::mutate(gendist = as.vector(unlist(ibd)))
    return(ret)
  }


IBDdat_tidy_out_2_discdat <- function(N, demeNames, locatcomb, IBDdat){
    #......................
    # memberships
    #......................
    # individuals and deme count
    membership_x <- tibble::tibble(smpl1 = 1:sum(N),
                                   deme1 = as.character(rep(demeNames, N)))
    membership_y <- tibble::tibble(smpl2 = 1:sum(N),
                                   deme2 = as.character(rep(demeNames, N)))
    #......................
    # expand locat comb
    #......................
    # confirm compat
    locatcomb <- locatcomb %>%
      dplyr::mutate(deme1 = as.character(deme1),
                    deme2 = as.character(deme2)) %>%
      dplyr::select(c("deme1", "deme2", "geodist"))

    locatcombexpand <- locatcomb
    colnames(locatcombexpand) <- c("deme2", "deme1", "geodist")
    locatcomb <- dplyr::bind_rows(locatcomb, locatcombexpand)

    #......................
    # bring together
    #......................
    discdat <- IBDdat %>%
      dplyr::left_join(., membership_x, by = "smpl1") %>%
      dplyr::left_join(., membership_y, by = "smpl2") %>%
      dplyr::left_join(., locatcomb, by = c("deme1", "deme2")) %>%
      dplyr::select(c("smpl1", "smpl2", "deme1", "deme2", "gendist", "geodist"))
    return(discdat)
  }

