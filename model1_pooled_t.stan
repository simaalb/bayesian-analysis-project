data {
  int<lower=1> N;
  vector[N] y;
  vector[N] bill_depth_s;
  array[N] int<lower=0, upper=1> sex_male;

  real prior_alpha_mean;
  real<lower=0> prior_alpha_sd;
  real<lower=0> prior_beta_sd;
  real<lower=0> prior_sigma_sd;
}

parameters {
  real alpha;
  real beta_depth;
  real beta_sex;
  real<lower=0> sigma;

  // Degrees of freedom for Student-t likelihood
  real<lower=2> nu;
}

model {
  alpha ~ normal(prior_alpha_mean, prior_alpha_sd);
  beta_depth ~ normal(0, prior_beta_sd);
  beta_sex ~ normal(0, prior_beta_sd);
  sigma ~ normal(0, prior_sigma_sd);

  // Prior for degrees of freedom
  nu ~ gamma(2, 0.1);

  y ~ student_t(
    nu,
    alpha + beta_depth * bill_depth_s + beta_sex * to_vector(sex_male),
    sigma
  );
}

generated quantities {
  vector[N] log_lik;
  array[N] real y_rep;

  for (i in 1:N) {
    real mu_i = alpha + beta_depth * bill_depth_s[i]
                + beta_sex * sex_male[i];

    log_lik[i] = student_t_lpdf(y[i] | nu, mu_i, sigma);
    y_rep[i] = student_t_rng(nu, mu_i, sigma);
  }
}
