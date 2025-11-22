#------------------------------------------------
#' @title Draw from Beta-binomial distribution
#'
#' @description Draw from Beta-binomial distribution.
#'
#' @param n number of draws.
#' @param k number of binomial trials.
#' @param alpha first shape parameter of beta distribution.
#' @param beta second shape parameter of beta distribution.
#'
#' @export

rbetabinom <- function(n = 1, k = 10, alpha = 1, beta = 1) {
  p <- rbeta(n = n, shape1 = alpha, shape2 = beta)
  ret <- rbinom(n = n, size = k, prob = p)
  return(ret)
}

#' @title DISCent Forward Simulator
#' @param true_params named vector; True F_i and M values to simulate. At least one element must be named "m".
#' @param geodist_matrix named matrix of geodistances; must be a nxn geodistance matrix, with n corresponding to number of demes and deme names matching the true_param names
#' @param samples_per_deme integer vector; number of pairs within each deme
#' @param overdispersion numeric; overdispersion parameter in the beta-binomial model; Default: 200
#' @description Uses the DISCent equation, \eqn{r_ij = (f_i + f_j/2)e^(-d_{ij}/m)} to simulate
#' data expected genetic distances between sample pairs within demes. We assume a beta-binomial
#' model to create realistic "noise" amongst the pairs, where \eqn{Y|p ~ binom(n,p)} such that
#' \eqn{p ~ beta(\mu * \phi, (1 - \mu) * \phi)}, where \eqn{\mu} is average relatedness (defined by \eqn{\rij}) and \eqn{\phi}
#' is a concentration parameter of that relatedness. In our instance the beta-binomial returns counts (ie number of related segments),
#' we fix the segment denominator at 100. Overdispersion can be adjusted to produce more or less correlation
#' of relatedness
#'
#'
#' @return A dataframe of class "DISCsim" containing columns:
#'   \itemize{
#'     \item \code{Final_Fis}: Final inbreeding coefficient estimates
#'     \item \code{Final_m}: Final migration rate estimate
#'     \item \code{deme_key}: Mapping of deme names to indices
#'     \item \code{cost}: Final cost function value(s)
#'   }
#' @export

run_forward_disc <- function(true_params, geodist_matrix, samples_per_deme, overdispersion = 200) {
  #............................................................
  # checks
  #............................................................
  assert_vector(true_params)
  assert_symmetric_matrix(geodist_matrix)
  assert_eq(any(names(true_params) %in% "m"), TRUE, message = "One element in your true_params vector must be named m")
  assert_leq(length(true_params[names(true_params) != "m"]), nrow(geodist_matrix),
             message = "Your geodistance matrix and true_params vectors are not aligned")
  assert_eq(all(names(true_params[names(true_params) != "m"]) %in% rownames(geodist_matrix)), TRUE,
            message = "Your geodistance matrix and true_params vector are not named properly")
  assert_eq(all(names(true_params[names(true_params) != "m"]) %in% colnames(geodist_matrix)), TRUE,
            message = "Your geodistance matrix and true_params vector are not named properly")
  assert_int(samples_per_deme)
  assert_vector(samples_per_deme)
  assert_leq(length(true_params[names(true_params) != "m"]), length(samples_per_deme),
             message = "Samples per deme must be same length as number of demes")
  assert_numeric(overdispersion)

  #............................................................
  # setup (func, storage, etc)
  #............................................................
  rijcalc <- function(from, to, geodist, fi, fj, m){
    ((fi + fj)/2) * exp(-geodist/m)
  }

  #............................................................
  # core
  #............................................................
  #......................
  # arrange data
  #......................
  geodist_matrix_long <- geodist_matrix %>%
    as.data.frame() %>%
    dplyr::mutate(from = rownames(.)) %>%
    tidyr::pivot_longer(., cols = -from, names_to = "to", values_to = "geodist")

  # upper triangle only
  geodist_matrix_long <- geodist_matrix_long %>%
    mutate(utkey = paste0(pmin(from, to), pmax(from, to))) %>%
    dplyr::filter(!duplicated(utkey)) %>%
    dplyr::select(-c("utkey"))

  # tidy
  keyfi <- tibble::tibble(from = names(true_params[names(true_params) != "m"]),
                          fi = true_params[names(true_params) != "m"])
  keyfj <- tibble::tibble(to = names(true_params[names(true_params) != "m"]),
                          fj = true_params[names(true_params) != "m"])

  demekey <- expand_grid(from = names(true_params[names(true_params) != "m"]),
                         to = names(true_params[names(true_params) != "m"])) %>%
    dplyr::left_join(., keyfi, by = "from") %>%
    dplyr::left_join(., keyfj, by = "to")

  # now we can create sim dat
  simdat <- dplyr::full_join(geodist_matrix_long, demekey, by = c("from", "to")) %>%
    dplyr::filter(from != to) %>%
    dplyr::mutate(m = true_params[names(true_params) == "m"]) %>%
    dplyr::mutate(expected_gendist = purrr::pmap_dbl(., rijcalc)) %>%
    dplyr::rename(deme1 = from,
                  deme2 = to)

  #......................
  # pair-level
  #......................
  demekey1 <- tibble::tibble(
    deme1 = rep(names(true_params[names(true_params) != "m"]), times = samples_per_deme),
    smpl1 = 1:sum(samples_per_deme)
  )
  demekey2 <- tibble::tibble(
    deme2 = rep(names(true_params[names(true_params) != "m"]), times = samples_per_deme),
    smpl2 = 1:sum(samples_per_deme)
  )
  pairs <- tidyr::expand_grid(smpl1 = demekey1$smpl1, smpl2 = demekey2$smpl2)
  # bring it all together
  fulldat <- pairs %>%
    dplyr::left_join(., demekey1, by = "smpl1") %>%
    dplyr::left_join(., demekey2, by = "smpl2") %>%
    dplyr::left_join(., simdat, by = c("deme1", "deme2")) %>%
    dplyr::filter(deme1 != deme2) %>%
    dplyr::filter(!is.na(geodist)) # upper triangle from above

  #......................
  # beta-binomial model
  #......................
  demepairs <- fulldat %>%
    dplyr::group_by(deme1, deme2) %>%
    tidyr::nest()
  # wrapper for beta binomial
  rbbwrap <- function(data, overdispersion){
    # get mu
    mu <- max(data$expected_gendist) # single value
    # run beta-binomial
    data$gendist <- rbetabinom(n = nrow(data),
                               k = overdispersion,
                               alpha = mu*overdispersion,
                               beta =  (1-mu)*overdispersion) / 100
    return(data)
  }

  fulldat <- demepairs %>%
    dplyr::mutate(data = purrr::map(data, rbbwrap, overdispersion = overdispersion)) %>%
    tidyr::unnest(cols = data) %>%
    dplyr::ungroup()


  #............................................................
  # tidy up out
  #............................................................
  discdat <- fulldat %>%
    dplyr::select(c("smpl1", "smpl2", "deme1", "deme2", "gendist", "geodist"))
  return(list(
    discdat = discdat,
    fulldat = fulldat
  ))
}

