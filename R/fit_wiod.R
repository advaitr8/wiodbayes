#' WIOD Bayesian modeling
#'
#' This function allows you to run a variety of Bayesian (hierarchical) models on WIOD data and key statistics to unearth persistent statistical patterns across time and countries
#' @param statistic this is the data
#' @param model is the type of model you want to run, choose from normal, lognormal, gamma, Weibull and exponential
#' @param levels_id_vector is the vector that shows which subsample a particular observation belongs to
#' @param pooling is the level of pooling desired, choose from complete, partial or none
#' @keywords wiodbayes
#' @export
#' @examples
#' sample examples here
#'

fit_wiod <- function(statistic,
                     levels_id_vector = 0,
                     model,
                     pooling){
  library("rstan")
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = TRUE)
  library("loo")
  #statistic is the vector of data
  #levels_id_vector is the vector that shows which subsample a particular observation belongs to
  #model specifies the type of model to run
  #pooling is the nature of pooling required across individual subsamples, the default is complete pooling
  N <- length(statistic)
  n_levels <- length(unique(levels_id_vector))
  if(model == "lognormal" || model == "weibull"){
    levels_id_vector <- levels_id_vector[statistic != 0]
    statistic <- statistic[statistic != 0]
    N <- length(statistic)
  }

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
      generated quantities{
        # real statistic_pred[N];
        vector[N] log_lik;
        for(i in 1:N){
          # statistic_pred[i] = normal_rng(mu,sigma);
          log_lik[i] = normal_lpdf(statistic[i] | mu, sigma);
        }
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
        real <lower = 0> tau_mu_all_levels;
        real <lower = 0> sigma[n_levels];
        real <lower = 0> sigma_all_levels;
        real <lower = 0> tau_sigma_all_levels;
      }
      model{
        for(i in 1:N){
          statistic[i] ~ normal(mu[levels_id_vector[i]],sigma[levels_id_vector[i]]);
      }
        mu ~ normal(mu_all_levels,tau_mu_all_levels);
        mu_all_levels ~ normal(0,5);
        tau_mu_all_levels ~ cauchy(0,20);
        sigma ~ cauchy(sigma_all_levels,tau_sigma_all_levels);
        sigma_all_levels ~ normal(0,0.1);
        tau_sigma_all_levels ~ normal(10,0.1);
      }
      generated quantities{
        # real statistic_pred[N];
        vector[N] log_lik;
        for(i in 1:N){
          # statistic_pred[i] = normal_rng(mu[levels_id_vector[i]],sigma[levels_id_vector[i]]);
          log_lik[i] = normal_lpdf(statistic[i] | mu[levels_id_vector[i]], sigma[levels_id_vector[i]]);
        }
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
          mu[levels_id_vector[i]] ~ normal(0,10);
          sigma[levels_id_vector[i]] ~ cauchy(0,20);
        }
      }
      generated quantities{
        real statistic_pred[N];
        vector[N] log_lik;
        for(i in 1:N){
          statistic_pred[i] = normal_rng(mu[levels_id_vector[i]],sigma[levels_id_vector[i]]);
          log_lik[i] = normal_lpdf(statistic[i] | mu[levels_id_vector[i]], sigma[levels_id_vector[i]]);
        }
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
      generated quantities{
        real statistic_pred[N];
        vector[N] log_lik;
        for(i in 1:N){
          statistic_pred[i] = lognormal_rng(mu,sigma);
          log_lik[i] = lognormal_lpdf(statistic[i] | mu, sigma);
        }
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
        real <lower = 0> sigma_all_levels;
        real <lower = 0> tau_sigma_all_levels;
      }
      model{
        for(i in 1:N){
          statistic[i] ~ lognormal(mu[levels_id_vector[i]],sigma[levels_id_vector[i]]);
      }
        mu ~ normal(mu_all_levels,tau_all_levels);
        mu_all_levels ~ normal(0,5);
        tau_all_levels ~ cauchy(0,20);
        sigma ~ cauchy(sigma_all_levels,tau_sigma_all_levels);
        sigma_all_levels ~ normal(0,0.1);
        tau_sigma_all_levels ~ normal(10,0.1);
      }
      generated quantities{
        # real statistic_pred[N];
        vector[N] log_lik;
        for(i in 1:N){
          # statistic_pred[i] = lognormal_rng(mu[levels_id_vector[i]],sigma[levels_id_vector[i]]);
          log_lik[i] = lognormal_lpdf(statistic[i] | mu[levels_id_vector[i]], sigma[levels_id_vector[i]]);
        }
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
          mu[levels_id_vector[i]] ~ normal(0,10);
          sigma[levels_id_vector[i]] ~ cauchy(0,20);
        }
      }
      generated quantities{
        real statistic_pred[N];
        vector[N] log_lik;
        for(i in 1:N){
          statistic_pred[i] = lognormal_rng(mu[levels_id_vector[i]],sigma[levels_id_vector[i]]);
          log_lik[i] = lognormal_lpdf(statistic[i] | mu[levels_id_vector[i]], sigma[levels_id_vector[i]]);
        }
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

  #####################
  # Skewnormal - Complete
  ####################
  if(model == "skew_normal"){
    if(pooling == "complete"){
      skew_normal_complete_model <- "
      data{
        int N;
        real statistic[N];
      }
      parameters{
        real mu;
        real <lower = 0> omega;
        real alpha;
      }
      model{
        for(i in 1:N){
          statistic[i] ~ skew_normal(mu, omega, alpha);
      }
        mu ~ normal(0,10);
        omega ~ normal(0,2);
        alpha ~ normal(0,0.5);
      }
      generated quantities{
        real statistic_pred[N];
        vector[N] log_lik;
        for(i in 1:N){
          statistic_pred[i] = skew_normal_rng(mu,sigma);
          log_lik[i] = skew_normal_lpdf(statistic[i] | mu, omega, alpha);
        }
      }
    "
    stanfit_object <- stan(model_code = skew_normal_complete_model,
                           data = list("statistic",
                                       "N"),
                           iter = 1000,
                           chains = 3)
    }
}


  #####################
  # Skewnormal - Partial
  ####################
  if(model == "skew_normal"){
    if(pooling == "partial"){
      skew_normal_partial_model <- "
      data{
        int N;
        int n_levels;
        int levels_id_vector[N];
        real statistic[N];
      }

      parameters{
        real mu[n_levels];
        # real mu_all_levels;
        # real <lower = 0> tau_mu_all_levels;
        real <lower = 0.1> omega[n_levels];
        # real <lower = 0.1> omega_all_levels;
        # real <lower = 0> tau_omega_all_levels;
        real alpha[n_levels];
        # real alpha_all_levels;
        # real <lower = 0> tau_alpha_all_levels;
      }

      model{
        for(i in 1:N){
          statistic[i] ~ skew_normal(mu[levels_id_vector[i]],
                                     omega[levels_id_vector[i]],
                                     alpha[levels_id_vector[i]]);
        }
        # mu ~ normal(mu_all_levels,tau_mu_all_levels);
        # omega ~ normal(omega_all_levels, tau_omega_all_levels);
        # alpha ~ normal(alpha_all_levels,tau_alpha_all_levels);
        # mu_all_levels ~ normal(0,5);
        # omega_all_levels ~ gamma(2,2);
        # alpha_all_levels ~ normal(0,5);
        # tau_mu_all_levels ~ cauchy(0,10);
        # tau_omega_all_levels ~ cauchy(0,10);
        # tau_alpha_all_levels ~ cauchy(0,10);
      }
      "
      stanfit_object <- stan(model_code = skew_normal_partial_model,
                             data = list("statistic",
                                         "levels_id_vector",
                                         "n_levels",
                                         "N"),
                             iter = 1000,
                             chains = 3)
    }
  }

  #####################
  # Skewnormal - None
  ####################
  if(model == "skew_normal"){
    if(pooling == "none"){
      skew_normal_none_model <- "
      data{
        int N;
        int n_levels;
        int levels_id_vector[N];
        real statistic[N];
      }
      parameters{
        real mu[n_levels];
        real <lower = 0> omega[n_levels];
        real alpha[n_levels];
      }
      model{
        for(i in 1:N){
          statistic[i] ~ skew_normal(mu[levels_id_vector[i]],
                                     omega[levels_id_vector[i]],
                                     alpha[levels_id_vector[i]]);
          mu[levels_id_vector[i]] ~ normal(0,10);
          omega[levels_id_vector[i]] ~ normal(0,5);
          alpha[levels_id_vector[i]] ~ normal(0,1);
        }
      }
      generated quantities{
        real statistic_pred[N];
        vector[N] log_lik;
        for(i in 1:N){
          statistic_pred[i] = skew_normal_rng(mu[levels_id_vector[i]],sigma[levels_id_vector[i]]);
          log_lik[i] = skew_normal_lpdf(statistic[i] | mu[levels_id_vector[i]], sigma[levels_id_vector[i]]);
      }
    }
      "
      stanfit_object <- stan(model_code = skew_normal_none_model,
                             data = list("statistic",
                                         "levels_id_vector",
                                         "n_levels",
                                         "N"),
                             iter = 1000,
                             chains = 3)
    }
}

  #####################
  # Gamma - Complete
  ####################

  if(model == "gamma"){
    if(pooling == "complete"){
      gamma_complete_model <- "
      data{
        int N;
        real statistic[N];
      }
      parameters{
        real <lower = 0> alpha;
        real <lower =0> beta;
      }
      model{
        for(i in 1:N){
          statistic[i] ~ gamma(alpha,beta);
        }
        alpha ~ exponential(3);
        beta ~ exponential(3);
      }
      generated quantities{
        real statistic_pred[N];
        vector[N] log_lik;
        for(i in 1:N){
          statistic_pred[i] = gamma_rng(alpha, beta);
          log_lik[i] = gamma_lpdf(statistic[i] | alpha, beta);
        }
      }
      "
      stanfit_object <- stan(model_code = gamma_complete_model,
                             data = list("statistic",
                                         "N"),
                             iter = 1000,
                             chains = 3)
    }
  }
  #####################
  # Gamma - Partial
  ####################

  if(model == "gamma"){
    if(pooling == "partial"){
      gamma_partial_model <- "
      data{
        int N;
        int n_levels;
        int levels_id_vector[N];
        real statistic[N];
      }
      parameters{
        real <lower = 0.1> alpha[n_levels];
        real <lower = 0.1> lambda_alpha;
        real <lower = 0.1> beta[n_levels];
        real <lower = 0.1> lambda_beta;
      }
      model{
        for(i in 1:N){
          statistic[i] ~ gamma(alpha[levels_id_vector[i]], beta[levels_id_vector[i]]);
        }
        alpha ~ exponential(lambda_alpha);
        lambda_alpha ~ normal(0,5);
        beta ~ exponential(lambda_beta);
        lambda_beta ~ normal(0,5);
      }
      generated quantities{
        # real statistic_pred[N];
        vector[N] log_lik;
        for(i in 1:N){
          # statistic_pred[i] = gamma_rng(alpha[levels_id_vector[i]], beta[levels_id_vector[i]]);
          log_lik[i] = gamma_lpdf(statistic[i] | alpha[levels_id_vector[i]], beta[levels_id_vector[i]]);
        }
      }
      "
      stanfit_object <- stan(model_code = gamma_partial_model,
                             data = list("statistic",
                                         "levels_id_vector",
                                         "n_levels",
                                         "N"),
                             iter = 1000,
                             chains = 3)
    }
}

  #####################
  # Gamma - None
  ####################
  if(model == "gamma"){
    if(pooling == "none"){
      gamma_none_model <- "
      data{
        int N;
        int n_levels;
        int levels_id_vector[N];
        real statistic[N];
      }
      parameters{
        real <lower = 0.1, upper = 2> alpha[n_levels];
        real <lower = 0.1, upper = 2> beta[n_levels];
      }
      model{
        for(i in 1:N){
          statistic[i] ~ gamma(alpha[levels_id_vector[i]], beta[levels_id_vector[i]]);
          alpha[levels_id_vector[i]] ~ exponential(3);
          beta[levels_id_vector[i]] ~ exponential(3);
        }
      }
      generated quantities{
        real statistic_pred[N];
        vector[N] log_lik;
        for(i in 1:N){
          statistic_pred[i] = gamma_rng(alpha[levels_id_vector[i]], beta[levels_id_vector[i]]);
          log_lik[i] = gamma_lpdf(statistic[i] | alpha[levels_id_vector[i]], beta[levels_id_vector[i]]);
      }
    }"
      stanfit_object <- stan(model_code = gamma_none_model,
                             data = list("statistic",
                                         "levels_id_vector",
                                         "n_levels",
                                         "N"),
                             iter = 1000,
                             chains = 3)
    }
}

  #####################
  # Exponential - Complete
  ####################
  if(model == "exponential"){
    if(pooling == "complete"){
      exponential_complete_model <- "
       data{
        int N;
        real statistic[N];
      }
      parameters{
        real lambda;
      }
      model{
        for(i in 1:N){
          statistic[i] ~ exponential(lambda);
        }
        lambda ~ exponential(2);
      }
      generated quantities{
        real statistic_pred[N];
        vector[N] log_lik;
        for(i in 1:N){
          statistic_pred[i] = exponential_rng(lambda);
          log_lik[i] = exponential_lpdf(statistic[i] | lambda);
        }
      }
      "
      stanfit_object <- stan(model_code = exponential_complete_model,
                             data = list("statistic",
                                         "N"),
                             iter = 1000,
                             chains = 3)
    }
  }

  #####################
  # Exponential - Partial
  ####################
  if(model == "exponential"){
    if(pooling == "partial"){
      exponential_partial_model <- "
      data{
        int N;
        int n_levels;
        int levels_id_vector[N];
        real statistic[N];
    }
      parameters{
        real <lower = 0.1> lambda[n_levels];
        real <lower = 0.1> theta;
      }
      model{
        for(i in 1:N){
          statistic[i] ~ exponential(lambda[levels_id_vector[i]]);
        }
        lambda ~ exponential(theta);
        theta ~ exponential(2);
      }
      generated quantities{
        # real statistic_pred[N];
        vector[N] log_lik;
        for(i in 1:N){
          # statistic_pred[i] = exponential_rng(lambda[levels_id_vector[i]]);
          log_lik[i] = exponential_lpdf(statistic[i] | lambda[levels_id_vector[i]]);
        }
      }
      "
      stanfit_object <- stan(model_code = exponential_partial_model,
                             data = list("statistic",
                                         "levels_id_vector",
                                         "n_levels",
                                         "N"),
                             iter = 1000,
                             chains = 3)
    }
  }

  #####################
  # Exponential - None
  ####################
  if(model == "exponential"){
    if(pooling == "none"){
      exponential_none_model <- "
      data{
        int N;
        int n_levels;
        int levels_id_vector[N];
        real statistic[N];
      }
      parameters{
        real <lower = 0.1> lambda[n_levels];
      }
      model{
        for(i in 1:N){
          statistic[i] ~ exponential(lambda[levels_id_vector[i]]);
          lambda[levels_id_vector[i]] ~ exponential(2);
        }
      }
      generated quantities{
        real statistic_pred[N];
        vector[N] log_lik;
        for(i in 1:N){
          statistic_pred[i] = exponential_rng(lambda[levels_id_vector[i]]);
          log_lik[i] = exponential_lpdf(statistic[i] | lambda[levels_id_vector[i]]);
        }
      }

      "
      stanfit_object <- stan(model_code = exponential_none_model,
                             data = list("statistic",
                                         "levels_id_vector",
                                         "n_levels",
                                         "N"),
                             iter = 1000,
                             chains = 3)
    }
}

  #####################
  # Weibull - Complete
  ####################
  if(model == "weibull"){
    if(pooling == "complete"){
      weibull_complete_model <- "
      data{
        int N;
        real statistic[N];
      }
      parameters{
        real <lower = 0.1> alpha;
        real <lower = 0.1> sigma;
      }
      model{
        for(i in 1:N){
          statistic[i] ~ weibull(alpha,sigma);
      }
        alpha ~ exponential(2);
        sigma ~ exponential(2);
      }
      generated quantities{
        real statistic_pred[N];
        vector[N] log_lik;
        for(i in 1:N){
          statistic_pred[i] = weibull_rng(alpha, sigma);
          log_lik[i] = weibull_lpdf(statistic[i] | alpha, sigma);
        }
      }
      "
      stanfit_object <- stan(model_code = weibull_complete_model,
                             data = list("statistic",
                                         "N"),
                             iter = 1000,
                             chains = 3)
    }
  }

  #####################
  # Weibull - Partial
  ####################
  if(model == "weibull"){
    if(pooling == "partial"){
      weibull_partial_model <- "
      data{
        int N;
        int n_levels;
        int levels_id_vector[N];
        real statistic[N];
      }
      parameters{
        real <lower = 0.1> alpha[n_levels];
        real <lower = 0.1> lambda_alpha;
        real <lower = 0.1> sigma[n_levels];
        real <lower = 0.1> lambda_sigma;
      }
      model{
        for(i in 1:N){
          statistic[i] ~ weibull(alpha[levels_id_vector[i]], sigma[levels_id_vector[i]]);
      }
        alpha ~ exponential(lambda_alpha);
        lambda_alpha ~ exponential(0.5);
        sigma ~ exponential(lambda_sigma);
        lambda_sigma ~ exponential(0.5);
      }
      generated quantities{
        # real statistic_pred[N];
        vector[N] log_lik;
        for(i in 1:N){
          # statistic_pred[i] = weibull_rng(alpha[levels_id_vector[i]], sigma[levels_id_vector[i]]);
          log_lik[i] = weibull_lpdf(statistic[i] | alpha[levels_id_vector[i]], sigma[levels_id_vector[i]]);
        }
      }
      "
      stanfit_object <- stan(model_code = weibull_partial_model,
                             data = list("statistic",
                                         "levels_id_vector",
                                         "n_levels",
                                         "N"),
                             iter = 1000,
                             chains = 3)
    }
}

  #####################
  # Weibull - None
  ####################
  if(model == "weibull"){
    if(pooling == "none"){
      weibull_none_model <- "
      data{
        int N;
        int n_levels;
        int levels_id_vector[N];
        real statistic[N];
      }
      parameters{
        real <lower = 0.1> alpha[n_levels];
        real <lower = 0.1> sigma[n_levels];
      }
      model{
        for(i in 1:N){
          statistic[i] ~ weibull(alpha[levels_id_vector[i]], sigma[levels_id_vector[i]]);
          alpha[levels_id_vector[i]] ~ exponential(2);
          sigma[levels_id_vector[i]] ~ exponential(2);
        }
      }
      generated quantities{
        real statistic_pred[N];
        vector[N] log_lik;
        for(i in 1:N){
          statistic_pred[i] = weibull_rng(alpha[levels_id_vector[i]], sigma[levels_id_vector[i]]);
          log_lik[i] = weibull_lpdf(statistic[i] | alpha[levels_id_vector[i]], sigma[levels_id_vector[i]]);
      }
    }
      "
      stanfit_object <- stan(model_code = weibull_none_model,
                             data = list("statistic",
                                         "levels_id_vector",
                                         "n_levels",
                                         "N"),
                             iter = 1000,
                             chains = 3)
    }
  }

  ######
  # returning the Stanfit object
  ######
  return(stanfit_object)
}


