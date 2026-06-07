rm(list = ls())

library(rstan)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(bayesplot)
library(loo)
library(posterior)

theme_set(theme_minimal())
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

#--------1- prepare data---------

penguins <- readRDS("penguins.RDS")

#keep only the variables used in the models and remove incomplete rows
penguins <- penguins %>%
  select(species, sex, bill_length, bill_depth) %>%
  filter(!is.na(species), !is.na(sex), !is.na(bill_length), !is.na(bill_depth))

penguins$species <- factor(penguins$species,levels = c("Adelie", "Chinstrap", "Gentoo"))
penguins$sex <- factor(penguins$sex)

#exploratory summaries
print(head(penguins))
print(summary(penguins))

summary_by_group <- penguins %>%
  group_by(species, sex) %>%
  summarise(
    n = n(),
    mean_bill_length = mean(bill_length),
    sd_bill_length= sd(bill_length),
    mean_bill_depth = mean(bill_depth),
    sd_bill_depth= sd(bill_depth),
    .groups = "drop"
  )
print(summary_by_group)

#plots
ggplot(penguins, aes(x = bill_length, fill = species)) +
  geom_histogram(binwidth = 1, color = "white") +
  facet_wrap(~species) +
  labs(title = "Bill length by species", x = "Bill length (mm)", y = "Count") +
  theme(legend.position = "none")

ggplot(penguins, aes(x = bill_depth, y = bill_length, color = species)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Bill length vs bill depth", x = "Bill depth (mm)", y = "Bill length (mm)")

ggplot(penguins, aes(x = species, y = bill_length, fill = sex)) +
  geom_boxplot() +
  labs(title = "Bill Length by Species and Sex")

#variables for Stan
penguins$species_id <- as.integer(penguins$species)
penguins$sex_male <- as.integer(penguins$sex == "male")
penguins$bill_depth_s <- as.numeric(scale(penguins$bill_depth))

N <- nrow(penguins)
J <- nlevels(penguins$species)
y <- penguins$bill_length
bill_depth_s <- penguins$bill_depth_s
sex_male <- penguins$sex_male
species <- penguins$species_id

cat("N =", N, " J =", J, "\n")

#----------2-data lists for Stan----------

stan_data_m1 <- list(
  N = N, y = y, bill_depth_s = bill_depth_s, sex_male = sex_male,
  prior_alpha_mean = 40,
  prior_alpha_sd   = 10,
  prior_beta_sd    = 5,
  prior_sigma_sd   = 5
)

stan_data_m2 <- list(
  N = N, J = J, y = y, bill_depth_s = bill_depth_s, sex_male = sex_male, species = species,
  prior_alpha_mean = 40,
  prior_alpha_sd   = 10,
  prior_beta_sd    = 5,
  prior_sigma_sd   = 5
)

stan_data_m3 <- list(
  N = N, J = J, y = y, bill_depth_s = bill_depth_s, sex_male = sex_male, species = species,
  prior_alpha_mean = 40,
  prior_alpha_sd   = 10,
  prior_group_sd   = 5,
  prior_beta_sd    = 5,
  prior_sigma_sd   = 5
)


#----------3- fit models---------

model1 = stan_model("model1_pooled.stan")
fit1 <- sampling(model1, data=stan_data_m1, iter=2000, chains=4, seed=123)

model2 = stan_model("model2_separate.stan")
fit2 <- sampling(model2, data=stan_data_m2, iter=2000, chains=4, seed=123)

model3 = stan_model("model3_hierarchical.stan")
fit3 <- sampling(model3, data=stan_data_m3, iter=2000, chains=4, seed=123, control = list(adapt_delta = 0.99,max_treedepth=15))

#----------4- convergence diagnostics-----------

print(fit1, pars = c("alpha", "beta_depth", "beta_sex", "sigma"))
print(fit2, pars = c("alpha", "beta_depth", "beta_sex", "sigma"))
print(fit3, pars = c("mu_alpha", "sigma_alpha", "mu_depth", "sigma_depth","alpha", "beta_depth", "beta_sex", "sigma"))

rstan::traceplot(fit1, pars = c("alpha", "beta_depth","beta_sex", "sigma", "lp__"), inc_warmup = FALSE)
rstan::traceplot(fit2, pars = c("alpha[1]", "alpha[2]", "alpha[3]", "beta_depth", "beta_sex", "sigma", "lp__"), inc_warmup = FALSE)
rstan::traceplot(fit3, pars = c("mu_alpha", "sigma_alpha", "mu_depth", "sigma_depth", "beta_sex", "sigma", "lp__"), inc_warmup = FALSE)

check_hmc_diagnostics(fit1)
check_hmc_diagnostics(fit2)
check_hmc_diagnostics(fit3)

count_divergences <- function(fit) {
  sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
  sapply(sampler_params, function(x) sum(x[, "divergent__"]))
}

cat("Model 1 divergences per chain:\n")
print(count_divergences(fit1))
cat("Model 2 divergences per chain:\n")
print(count_divergences(fit2))
cat("Model 3 divergences per chain:\n")
print(count_divergences(fit3))

summarise_selected <- function(fit, pars) {
  summarise_draws(as_draws_df(fit), ess_bulk, ess_tail, rhat) %>%
    filter(variable %in% pars)
}

print(summarise_selected(fit1, c("alpha", "beta_depth", "beta_sex", "sigma")))
print(summarise_selected(fit2, c("alpha[1]", "alpha[2]", "alpha[3]", "beta_depth", "beta_sex", "sigma")))
print(summarise_selected(fit3, c("mu_alpha", "sigma_alpha", "alpha[1]", "alpha[2]", "alpha[3]",
                                 "mu_depth", "sigma_depth", "beta_sex", "sigma")))

#-------------5- posterior predictive checks-----------

draws1 <- rstan::extract(fit1)
draws2 <- rstan::extract(fit2)
draws3 <- rstan::extract(fit3)
y_rep1 <- draws1$y_rep
y_rep2 <- draws2$y_rep
y_rep3 <- draws3$y_rep

ppc_dens_overlay(y, y_rep1[1:100, ]) + ggtitle("Model 1: pooled")
ppc_dens_overlay(y, y_rep2[1:100, ]) + ggtitle("Model 2: separate")
ppc_dens_overlay(y, y_rep3[1:100, ]) + ggtitle("Model 3: hierarchical")

ppc_stat(y, y_rep1, stat = "min") + ggtitle("Model 1: T(y) = min")
ppc_stat(y, y_rep2, stat = "min") + ggtitle("Model 2: T(y) = min")
ppc_stat(y, y_rep3, stat = "min") + ggtitle("Model 3: T(y) = min")

ppc_stat(y, y_rep1, stat = "max") + ggtitle("Model 1: T(y) = max")
ppc_stat(y, y_rep2, stat = "max") + ggtitle("Model 2: T(y) = max")
ppc_stat(y, y_rep3, stat = "max") + ggtitle("Model 3: T(y) = max")

color_scheme_set("brewer-Paired")
ppc_stat_2d(y, y_rep1, stat = c("mean", "sd")) + ggtitle("Model 1")
ppc_stat_2d(y, y_rep2, stat = c("mean", "sd")) + ggtitle("Model 2")
ppc_stat_2d(y, y_rep3, stat = c("mean", "sd")) + ggtitle("Model 3")
color_scheme_set()

#combined summary grid

color_scheme_set("brewer-Paired")

grid.arrange(
  ppc_dens_overlay(y, y_rep1[1:100, ]) + ggtitle("Density overlay") + theme(legend.position = "none"),
  ppc_stat(y, y_rep1, stat = "min") + ggtitle("T(y) = min") + theme(legend.position = "none"),
  ppc_stat(y, y_rep1, stat = "max") + ggtitle("T(y) = max") + theme(legend.position = "none"),
  ppc_stat_2d(y, y_rep1, stat = c("mean", "sd")) + ggtitle("T(y) = (mean, sd)") + theme(legend.position = "none"),
  ncol = 2, top = "PPC summary — Model 1: pooled"
)

grid.arrange(
  ppc_dens_overlay(y, y_rep2[1:100, ]) + ggtitle("Density overlay")+ theme(legend.position = "none"),
  ppc_stat(y, y_rep2, stat = "min") + ggtitle("T(y) = min") + theme(legend.position = "none"),
  ppc_stat(y, y_rep2, stat = "max") + ggtitle("T(y) = max")  + theme(legend.position = "none"),
  ppc_stat_2d(y, y_rep2, stat = c("mean", "sd")) + ggtitle("T(y) = (mean, sd)") + theme(legend.position = "none"),
  ncol = 2, top = "PPC summary — Model 2: separate"
)

grid.arrange(
  ppc_dens_overlay(y, y_rep3[1:100, ]) + ggtitle("Density overlay") + theme(legend.position = "none"),
  ppc_stat(y, y_rep3, stat = "min") + ggtitle("T(y) = min")+ theme(legend.position = "none"),
  ppc_stat(y, y_rep3, stat = "max") + ggtitle("T(y) = max") + theme(legend.position = "none"),
  ppc_stat_2d(y, y_rep3, stat = c("mean", "sd")) + ggtitle("T(y) = (mean, sd)") + theme(legend.position = "none"),
  ncol = 2, top = "PPC summary — Model 3: hierarchical"
)

color_scheme_set()

#----------6- LOO model comparison----------
log_lik1 <- extract_log_lik(fit1, merge_chains = FALSE)
log_lik2 <- extract_log_lik(fit2, merge_chains = FALSE)
log_lik3 <- extract_log_lik(fit3, merge_chains = FALSE)

loo1 <- loo(log_lik1)
loo2 <- loo(log_lik2)
loo3 <- loo(log_lik3)

print(loo1)
print(loo2)
print(loo3)
print(loo_compare(loo1, loo2, loo3))


#--------7- predictive accuracy summary--------

pred_summary <- function(y, y_rep, model_name) {
  y_pred <- colMeans(y_rep)
  data.frame(
    model = model_name,
    RMSE = sqrt(mean((y - y_pred)^2)),
    MAE = mean(abs(y - y_pred))
  )
}

pred_results <- bind_rows(
  pred_summary(y, y_rep1, "Model 1: pooled"),
  pred_summary(y, y_rep2, "Model 2: separate"),
  pred_summary(y, y_rep3, "Model 3: hierarchical")
)
print(pred_results)

pred_df <- data.frame(
  observed = rep(y, 3),
  predicted = c(colMeans(y_rep1),colMeans(y_rep2),colMeans(y_rep3)),
  species = rep(penguins$species, 3),
  sex = rep(penguins$sex, 3),
  model = rep(c("Model 1: pooled","Model 2: separate","Model 3: hierarchical"),
  each = N)
)

ggplot(pred_df, aes(x = observed,y = predicted, color = species, shape = sex)) +
  geom_point(alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  facet_wrap(~model) +
  labs(title = "Predicted vs observed bill length", x = "Observed bill length (mm)",
    y = "Predicted bill length (mm)",
    shape = "Sex"
  )

#--------8- sensitivity analysis--------

# 8.1 fit models

model1_t = stan_model("model1_pooled_t.stan")
fit1_t <- sampling(model1_t, data=stan_data_m1, iter=2000, chains=4, seed=123)

model2_t = stan_model("model2_separate_t.stan")
fit2_t <- sampling(model2_t, data=stan_data_m2, iter=2000, chains=4, seed=123)

model3_t = stan_model("model3_hierarchical_t.stan")
fit3_t <- sampling(model3_t, data=stan_data_m3, iter=2000, chains=4, seed=123, control=list(adapt_delta = 0.99, max_treedepth = 15))


# 8.2 convergence diagnostics

check_hmc_diagnostics(fit1_t)
check_hmc_diagnostics(fit2_t)
check_hmc_diagnostics(fit3_t)

print(fit1_t, pars = c("alpha", "beta_depth", "beta_sex", "sigma", "nu"))
print(fit2_t, pars = c("alpha", "beta_depth", "beta_sex", "sigma", "nu"))
print(fit3_t, pars = c("mu_alpha", "sigma_alpha", "mu_depth", "sigma_depth", "alpha", "beta_depth", "beta_sex", "sigma", "nu"))

# 8.3 LOO comparison

log_lik1_t <- extract_log_lik(fit1_t, merge_chains = FALSE)
log_lik2_t <- extract_log_lik(fit2_t, merge_chains = FALSE)
log_lik3_t <- extract_log_lik(fit3_t, merge_chains = FALSE)

loo1_t <- loo(log_lik1_t)
loo2_t <- loo(log_lik2_t)
loo3_t <- loo(log_lik3_t)

cat("Student-t Model 1 LOO:\n")
print(loo1_t)

cat("Student-t Model 2 LOO:\n")
print(loo2_t)

cat("Student-t Model 3 LOO:\n")
print(loo3_t)

cat("Normal vs Student-t: Model 1\n")
print(loo_compare(loo1, loo1_t))

cat("Normal vs Student-t: Model 2\n")
print(loo_compare(loo2, loo2_t))

cat("Normal vs Student-t: Model 3\n")
print(loo_compare(loo3, loo3_t))

cat("All six models together:\n")
print(loo_compare(loo1, loo2, loo3, loo1_t, loo2_t, loo3_t))

# 8.4 posterior comparison: key parameters

compare_posteriors <- function(fit_normal, fit_t, pars, model_name) {
  
  normal_sum <- as.data.frame(summary(fit_normal, pars = pars)$summary)
  normal_sum$parameter <- rownames(normal_sum)
  normal_sum$likelihood <- "Normal"
  
  t_sum <- as.data.frame(summary(fit_t, pars = pars)$summary)
  t_sum$parameter <- rownames(t_sum)
  t_sum$likelihood <- "Student-t"
  
  out <- bind_rows(normal_sum, t_sum) %>%
    select(parameter, likelihood, mean, sd, `2.5%`, `97.5%`)
  
  cat("\nPosterior comparison for", model_name, "\n")
  print(out)
  
  return(out)
}

post_comp_m1 <- compare_posteriors(
  fit1, fit1_t, pars = c("alpha", "beta_depth", "beta_sex", "sigma"), model_name = "Model 1"
)

post_comp_m2 <- compare_posteriors(
  fit2, fit2_t, pars = c("alpha[1]", "alpha[2]", "alpha[3]","beta_depth", "beta_sex", "sigma"), model_name = "Model 2"
)

post_comp_m3 <- compare_posteriors(
  fit3, fit3_t,
  pars = c("mu_alpha", "sigma_alpha","mu_depth", "sigma_depth","alpha[1]", "alpha[2]", "alpha[3]","beta_sex", "sigma"),
  model_name = "Model 3"
)

# 8.5 posterior predictive checks

draws1_t <- rstan::extract(fit1_t)
draws2_t <- rstan::extract(fit2_t)
draws3_t <- rstan::extract(fit3_t)

y_rep1_t <- draws1_t$y_rep
y_rep2_t <- draws2_t$y_rep
y_rep3_t <- draws3_t$y_rep

color_scheme_set("brewer-Paired")

grid.arrange(
  ppc_dens_overlay(y, y_rep1[1:100, ]) + ggtitle("Normal likelihood") +theme(legend.position = "none"),
  ppc_dens_overlay(y, y_rep1_t[1:100, ]) + ggtitle("Student-t likelihood") + theme(legend.position = "none"),
  ncol = 2,
  top = "Model 1 sensitivity: Normal vs Student-t likelihood"
)

grid.arrange(
  ppc_dens_overlay(y, y_rep2[1:100, ]) +ggtitle("Normal likelihood") + theme(legend.position = "none"),
  ppc_dens_overlay(y, y_rep2_t[1:100, ]) + ggtitle("Student-t likelihood") + theme(legend.position = "none"),
  ncol = 2,
  top = "Model 2 sensitivity: Normal vs Student-t likelihood"
)

grid.arrange(
  ppc_dens_overlay(y, y_rep3[1:100, ]) +ggtitle("Normal likelihood") + theme(legend.position = "none"),
  ppc_dens_overlay(y, y_rep3_t[1:100, ]) + ggtitle("Student-t likelihood") + theme(legend.position = "none"),
  ncol = 2,
  top = "Model 3 sensitivity: Normal vs Student-t likelihood"
)

color_scheme_set()
