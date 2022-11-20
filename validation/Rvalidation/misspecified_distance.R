add_misspec_dist <- function(prediscdat, seed = 48) {
  set.seed(seed)
  # extract out for liftover
  base <- prediscdat %>%
    dplyr::filter(modname == "IsoByDist")
  dists <- base$discdat$geodists[[1]]
  wrongdist <- rexp(n = length(dists), rate = mean(dists))
  # liftover
  base <- base %>%
    dplyr::mutate(discdat =
                    purrr::map(discdat,
                               function(x){ x$geodists <- wrongdist
                               return(x) }))
  # out
  return(dplyr::bind_rows(prediscdat, base))
}
