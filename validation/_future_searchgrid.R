#............................................................
# make search grid
#...........................................................
# look up tables
fstartsvec <- seq(from = 0.1, to = 0.9, by = 0.05)
mstartsvec <- 10^-seq(1,6)
flearnsvec <- 10^-seq(1,10)
mlearnsvec <- 10^-seq(5,15)
search_grid <- expand.grid(fstartsvec, mstartsvec, flearnsvec, mlearnsvec)
colnames(search_grid) <- c("fstart", "mstart", "f_learn", "m_learn")
# down sample
downsample <- 1e4
search_grid <- search_grid[sample(1:nrow(search_grid), size  = downsample, replace = F), ]
# template start
tempstart_params <- rep(0.1, 100)
names(tempstart_params) <- as.character(1:100)
tempstart_params <- c(tempstart_params, "m" = 1e-3)

# liftover to start param format
liftover_start_params <- function(fstart, mstart, start_param_template) {
  out <- start_param_template
  out[names(out) != "m"] <- fstart
  out[names(out) == "m"] <- mstart
  return(out)
}
search_grid <- search_grid %>%
  dplyr::mutate(start_params = purrr::map2(fstart, mstart, liftover_start_params,
                                           start_param_template = tempstart_params)) %>%
  dplyr::select(c("start_params", "f_learn", "m_learn"))
