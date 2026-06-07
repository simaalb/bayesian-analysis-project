data {
  int<lower=1> N; # observations
  vector[N] y; # response var bill length
  vector[N] bill_depth_s; 
  array[N] int<lower=0, upper=1> sex_male; 
  real prior_alpha_mean; #prior mean for intercept
  real<lower=0> prior_alpha_sd; #prior sd for intercept
  real<lower=0> prior_beta_sd; #prior sd for regression coef
  real<lower=0> prior_sigma_sd; #prior sd for residual std dev
}
parameters {
  real alpha;
  real beta_depth;
  real beta_sex;
  real<lower=0> sigma;
}
model {
  alpha ~ normal(prior_alpha_mean, prior_alpha_sd);
  beta_depth ~ normal(0, prior_beta_sd);
  beta_sex  ~ normal(0, prior_beta_sd);
  sigma ~ normal(0, prior_sigma_sd);
  y ~ normal(alpha + beta_depth * bill_depth_s + beta_sex * to_vector(sex_male), sigma);
}
generated quantities {
  vector[N] log_lik;
  array[N] real y_rep;
  for (i in 1:N) {
    real mu_i = alpha + beta_depth * bill_depth_s[i]
                + beta_sex * sex_male[i];
    log_lik[i] = normal_lpdf(y[i] | mu_i, sigma);
    y_rep[i] = normal_rng(mu_i, sigma);
  }
}
