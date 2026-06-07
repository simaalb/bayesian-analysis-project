data {
  int<lower=1> N;
  int<lower=1> J;
  vector[N] y;
  vector[N] bill_depth_s;
  array[N] int<lower=0, upper=1> sex_male;
  array[N] int<lower=1, upper=J> species;
  real prior_alpha_mean;
  real<lower=0> prior_alpha_sd;
  real<lower=0> prior_group_sd;
  real<lower=0> prior_beta_sd;
  real<lower=0> prior_sigma_sd;
}
parameters {
  real mu_alpha;
  real<lower=0> sigma_alpha;
  real mu_depth;
  real<lower=0> sigma_depth;

  // raw parameters — always Normal(0,1), no funnel
  vector[J] alpha_raw;
  vector[J] beta_depth_raw;

  real beta_sex;
  real<lower=0> sigma;
}
transformed parameters {
  // reconstruct the actual species-level parameters
  vector[J] alpha = mu_alpha + sigma_alpha * alpha_raw;
  vector[J] beta_depth = mu_depth + sigma_depth * beta_depth_raw;
}
model {
  mu_alpha ~ normal(prior_alpha_mean, prior_alpha_sd);
  sigma_alpha ~ normal(0, prior_group_sd);
  mu_depth ~ normal(0, prior_beta_sd);
  sigma_depth ~ normal(0, prior_group_sd);

  // raw priors — simple Normal(0,1), no dependence on sigma_alpha/sigma_depth
  alpha_raw ~ normal(0, 1);
  beta_depth_raw ~ normal(0, 1);

  beta_sex ~ normal(0, prior_beta_sd);
  sigma ~ normal(0, prior_sigma_sd);

  for (i in 1:N)
    y[i] ~ normal(alpha[species[i]]
                  + beta_depth[species[i]] * bill_depth_s[i]
                  + beta_sex * sex_male[i], sigma);
}
generated quantities {
  vector[N] log_lik;
  array[N] real y_rep;
  for (i in 1:N) {
    real mu_i = alpha[species[i]]
                + beta_depth[species[i]] * bill_depth_s[i]
                + beta_sex * sex_male[i];
    log_lik[i] = normal_lpdf(y[i] | mu_i, sigma);
    y_rep[i] = normal_rng(mu_i, sigma);
  }
}
