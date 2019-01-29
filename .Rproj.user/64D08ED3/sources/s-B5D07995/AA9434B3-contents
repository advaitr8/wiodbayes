#' WIOD Bayesian modeling
#'
#' This function allows you to run a variety of Bayesian (hierarchical) models on WIOD data and key statistics to unearth persistent statistical patterns across time and countries
#' @param statistic this is the data
#' @param model is the type of model you want to run, choose from normal, lognormal, gamma, Weibull and exponential
#' @param pooling is the level of pooling desired, choose from complete, partial or none
#' @keywords wiodbayes
#' @export
#' @examples
#' sample examples here
#'

fit_wiod <- function(statistic,
                     levels_id_vector,
                     model,
                     pooling){
  library("rstan")
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = TRUE)
  #statistic is the vector of data
  #model specifies the type of model to run
  #pooling is the nature of pooling required across individual subsamples, the default is complete pooling
  N <- length(statistic)
  n_levels <- length(unique(levels_id_vector))

  #####################
  # Normal - Complete
  ####################

  if(model == "normal"){
    if(pooling == "complete"){
      normal_complete_model <- "
      data{
        int N;
        real statistic[N];
    }
      parameters{
        real mu;
        real <lower = 0> sigma;
      }
      model{
        for(i in 1:N){
          statistic[i] ~ normal(mu,sigma);
      }
        mu ~ normal(0,10);
        sigma ~ cauchy(0,20);
      }
      "
      stanfit_object <- stan(model_code = normal_complete_model,
                              data = list("statistic",
                                          "N"),
                              iter = 1000,
                              chains = 3)
    }
  }
  #####################
  # Normal - Partial
  ####################
  if(model == "normal"){
    if(pooling == "partial"){
      normal_partial_model <- "
      data{
        int N;
        int n_levels;
        int levels_id_vector[N];
        real statistic[N];
    }
      parameters{
        real mu[n_levels];
        real mu_all_levels;
        real <lower = 0> tau_all_levels;
        real <lower = 0> sigma[n_levels];
      }
      model{
        for(i in 1:N){
          statistic[i] ~ normal(mu[levels_id_vector[i]],sigma[levels_id_vector[i]]);
      }
        mu ~ normal(mu_all_levels,tau_all_levels);
        mu_all_levels ~ normal(0,5);
        tau_all_levels ~ cauchy(0,20);
        sigma ~ cauchy(0,20);
      }
      "
      stanfit_object <- stan(model_code = normal_partial_model,
                             data = list("statistic",
                                         "levels_id_vector",
                                         "n_levels",
                                         "N"),
                             iter = 1000,
                             chains = 3)

    }
  }

  #####################
  # Normal - None
  ####################
  if(model == "normal"){
    if(pooling == "none"){
      normal_none_model <- "
      data{
        int N;
        int n_levels;
        int levels_id_vector[N];
        real statistic[N];
    }
      parameters{
        real mu[n_levels];
        real <lower = 0> sigma[n_levels];
      }
      model{
        for(i in 1:N){
          statistic[i] ~ normal(mu[levels_id_vector[i]],sigma[levels_id_vector[i]]);
      }
        mu ~ normal(0,10);
        sigma ~ cauchy(0,20);
      }
      "
      stanfit_object <- stan(model_code = normal_none_model,
                             data = list("statistic",
                                         "levels_id_vector",
                                         "n_levels",
                                         "N"),
                             iter = 1000,
                             chains = 3)

    }
  }
  if(model == "lognormal"){
    levels_id_vector <- levels_id_vector[statistic != 0]
    statistic <- statistic[statistic != 0]
    N <- length(statistic)
  }
  #####################
  # Lognormal - Complete
  ####################
  if(model == "lognormal"){
    if(pooling == "complete"){
      lognormal_complete_model <- "
      data{
        int N;
        real statistic[N];
    }
      parameters{
        real mu;
        real <lower = 0> sigma;
      }
      model{
        for(i in 1:N){
          statistic[i] ~ lognormal(mu,sigma);
      }
        mu ~ normal(0,10);
        sigma ~ cauchy(0,20);
      }
      "
      stanfit_object <- stan(model_code = lognormal_complete_model,
                             data = list("statistic",
                                         "levels_id_vector",
                                         "n_levels",
                                         "N"),
                             iter = 1000,
                             chains = 3)

    }
  }
  #####################
  # Lognormal - Partial
  ####################
  if(model == "lognormal"){
    if(pooling == "partial"){
      lognormal_partial_model <- "
      data{
        int N;
        int n_levels;
        int levels_id_vector[N];
        real statistic[N];
    }
      parameters{
        real mu[n_levels];
        real mu_all_levels;
        real <lower = 0> tau_all_levels;
        real <lower = 0> sigma[n_levels];
      }
      model{
        for(i in 1:N){
          statistic[i] ~ lognormal(mu[levels_id_vector[i]],sigma[levels_id_vector[i]]);
      }
        mu ~ normal(mu_all_levels,tau_all_levels);
        mu_all_levels ~ normal(0,5);
        tau_all_levels ~ cauchy(0,20);
        sigma ~ cauchy(0,20);
      }
      "
      stanfit_object <- stan(model_code = lognormal_partial_model,
                             data = list("statistic",
                                         "levels_id_vector",
                                         "n_levels",
                                         "N"),
                             iter = 1000,
                             chains = 3)
    }
  }
  #####################
  # Lognormal - None
  ####################
  if(model == "lognormal"){
    if(pooling == "none"){
      lognormal_none_model <- "
      data{
        int N;
        int n_levels;
        int levels_id_vector[N];
        real statistic[N];
    }
      parameters{
        real mu[n_levels];
        real <lower = 0> sigma[n_levels];
      }
      model{
        for(i in 1:N){
          statistic[i] ~ lognormal(mu[levels_id_vector[i]],sigma[levels_id_vector[i]]);
      }
        mu ~ normal(0,10);
        sigma ~ cauchy(0,20);
      }
      "
      stanfit_object <- stan(model_code = lognormal_none_model,
                             data = list("statistic",
                                         "levels_id_vector",
                                         "n_levels",
                                         "N"),
                             iter = 1000,
                             chains = 3)

  }
  }
  return(stanfit_object)
}


