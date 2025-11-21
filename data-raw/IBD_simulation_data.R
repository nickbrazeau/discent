## .................................................................................
## Purpose: Forward simulator for DISC
## .................................................................................
# Starting Parameters
ndemes <- 5
true_F <- round((rbeta(n = ndemes, shape1 = 2, shape2 = 2) + 0.05), 2)
names(true_F) <- LETTERS[1:ndemes]
true_M <- 500
names(true_M) <- "m"
# Make up some Geodistances, ignoring triangularity of geospatial data:
geodistmat <- matrix(NA, nrow = ndemes, ncol = ndemes)
geodistmat[upper.tri(geodistmat, diag = F)] <- rexp(n = sum(upper.tri(geodistmat, diag = F)), rate = 1/500)
geodistmat[lower.tri(geodistmat, diag = FALSE)] <- t(geodistmat)[lower.tri(geodistmat, diag = FALSE)]
diag(geodistmat) <- 0
rownames(geodistmat) <- colnames(geodistmat) <- LETTERS[1:ndemes]
geodist_matrix <- geodistmat

# run simulation
set.seed(48)
IBD_simulation_data <- run_forward_disc(true_params = c(true_F, true_M),
                                        geodist_matrix = geodistmat,
                                        samples_per_deme = c(3,3,3,3,3),
                                        overdispersion = 200)



#............................................................
# out
#...........................................................

usethis::use_data(IBD_simulation_data, overwrite = TRUE)
