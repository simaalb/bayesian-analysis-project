data {
  int<lower=1> N;
  int<lower=1> J; #species 
  vector[N] y;
  vector[N] bill_depth_s;
  array[N] int<lower=0, upper=1> sex_male;
  array[N] int<lower=1, upper=J> species;
  // prior hyperparameters
  real prior_alpha_mean;
  real<lower=0> prior_alpha_sd;
  real<lower=0> prior_beta_sd;
  real<lower=0> prior_sigma_sd;
}
parameters {
  vector[J] alpha;     // species-specific intercepts
  real beta_depth;     // shared slope
  real beta_sex;       // shared sex effect
  real<lower=0> sigma;
}
model {
  alpha ~ normal(prior_alpha_mean, prior_alpha_sd);
  beta_depth ~ normal(0, prior_beta_sd);
  beta_sex ~ normal(0, prior_beta_sd);
  sigma ~ normal(0, prior_sigma_sd);
  for (i in 1:N)
    y[i] ~ normal(alpha[species[i]]
                  + beta_depth * bill_depth_s[i]
                  + beta_sex   * sex_male[i], sigma);
}
generated quantities {
  vector[N] log_lik;
  array[N] real y_rep;
  for (i in 1:N) {
    real mu_i = alpha[species[i]]
                + beta_depth * bill_depth_s[i]
                + beta_sex   * sex_male[i];
    log_lik[i] = normal_lpdf(y[i] | mu_i, sigma);
    y_rep[i] = normal_rng(mu_i, sigma);
  }
}
