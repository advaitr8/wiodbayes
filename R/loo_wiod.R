#' Leave one out Cross Validation for WIOD Bayesian modeling
#'
#' This function allows you to compare Bayesian models using leave one out cross validation
#' @param stanfit_object this is the stanfit object that is the argument for the loo-cv function
#' @export
#' @examples
#' sample examples here
#'

loo_wiod <- function(stanfit_object){
  library("loo")
  # Extract pointwise log-likelihood and compute LOO
  log_lik <- extract_log_lik(stanfit_object,
                             merge_chains = F)

  # as of loo v2.0.0 we can optionally provide relative effective sample sizes
  # when calling loo, which allows for better estimates of the PSIS effective
  # sample sizes and Monte Carlo error
  r_eff <- relative_eff(exp(log_lik))

  loo_object <- loo(log_lik,
                    r_eff = r_eff,
                    cores = parallel::detectCores())
  return(loo_object)
}