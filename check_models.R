# Notes -------------------------------------------------------------------

# Master file to do:
# - Prior predictive checks
# - Fake data simulation (from priors)
# - Fit fake data (to see if we can recover parameters)
# For different models:
# - RW: random walk model
# - AR: autoregressive model (order 1, fixed effects)
# - MixedAR: mixed effect autoregressive model (order 1)
# - SSM: Hidden Markov Model (Gaussian measurement error and mixed autoregressive model for the latent dynamic)
# - SSMX: Hidden Markov Model with eXogeneous variables, following a horseshoe prior (two parametrisations for the horseshoe)

# NB: for AR/SSM, identifiability issue when alpha -> 1 (we can't determine S_inf); make sure draw is OK

# Initialisation ----------------------------------------------------------

rm(list = ls()) # Clear Workspace (but better to restart session)

set.seed(1744834965) # Reproducibility (Stan use a different seed)

library(tidyverse)
library(cowplot)
library(HuraultMisc) # Functions shared across projects
library(rstan)
rstan_options(auto_write = TRUE) # Save compiled model
options(mc.cores = parallel::detectCores()) # Parallel computing
source("functions.R") # Additional functions

#### OPTIONS
score <- "EASI"
model_name <- "SSM"
run_prior <- FALSE # prior distribution
run_fake <- FALSE # fit fake data
n_it <- 2000
n_chains <- 4
####

score <- match.arg(score, c("EASI", "oSCORAD", "SCORAD", "POEM"))
model_name <- match.arg(model_name, c("RW", "AR", "MixedAR", "SSM", "SSMX"))

# Results files
stan_code <- file.path("Models", paste0(model_name, ".stan"))
prior_file <- file.path("Results", paste0("prior_", score, "_", model_name, ".rds"))
par0_file <- file.path("Results", paste0("par0_", score, "_", model_name, ".rds"))
fake_file <- file.path("Results", paste0("fake_", score, "_", model_name, ".rds"))

# Model parameters
if (model_name == "RW") {
  param_pop <- c("sigma")
  param_ind <- c()
}
if (model_name == "AR") {
  param_pop <- c("sigma", "alpha", "S_inf", "b")
  param_ind <- c()
}
if (model_name == "MixedAR") {
  param_pop <- c("sigma",
                 "mu_alpha", "phi_alpha",
                 "mu_inf", "sigma_inf")
  param_ind <- c("alpha", "S_inf", "b")
}
if (model_name == "SSM") {
  param_pop <- c("sigma_tot", "rho2", "sigma_lat", "sigma_meas", "MDC",
                 "mu_alpha", "phi_alpha",
                 "mu_inf", "sigma_inf")
  param_ind <- c("alpha", "S_inf", "b")
}
if (model_name == "SSMX") {
  param_pop <- c("sigma_tot", "rho2", "sigma_lat", "sigma_meas", "MDC",
                 "mu_alpha", "phi_alpha",
                 "mu_inf", "sigma_inf",
                 "beta")
  param_ind <- c("alpha", "S_inf", "b", "f")
}

if (any(run_prior, run_fake)) {
  compiled_model <- stan_model(stan_code)
}

score_char <- data.frame(Score = c("SCORAD", "oSCORAD", "EASI", "POEM"),
                         Range = c(103, 83, 72, 28),
                         MCID = c(8.7, 8.2, 6.6, 3.4)) %>%
  filter(Score == score)

# Data --------------------------------------------------------------------
# Load data to generate fake data with similar characteristics
# Also to extract the matrix of biomarkers to avoid generating one

l <- load_dataset()
dp <- l$patient_data
dt <- l$severity_data
pt <- unique(dt[["Patient"]])
bio <- as.matrix(dp[, colnames(dp) != "Patient"]) # matrix of biomarkers (including treatment, age, sex...)

n_pt <- length(pt)
n_dur <- 24 / 2 + 1

# Prior predictive check ------------------------------------------------------

param_obs <- c("S_rep")
param <- c(param_pop, param_ind, param_obs)

data_prior <- list(
  N_obs = n_pt,
  N_mis = n_pt * (n_dur - 1),
  N_pt = n_pt,
  max_score = score_char$Range,
  idx_obs = ((1:n_pt) - 1) * (n_dur) + 1,
  idx_mis = setdiff(1:(n_pt * n_dur), ((1:n_pt) - 1) * (n_dur) + 1),
  S_obs = runif(n_pt, 0, score_char$Range),
  
  N_test = 0,
  idx_test = vector(),
  S_test = vector(),
  
  # For horseshoe
  p0 = 5,
  slab_scale = 1,
  slab_df = 5,
  N_cov = ncol(bio),
  X_cov = bio,
  parametrisation = 1,
  
  run = 0,
  rep = 1
)

if (run_prior) {
  fit_prior <- sampling(compiled_model,
                        data = data_prior,
                        pars = param,
                        iter = n_it,
                        chains = n_chains,
                        control = list(adapt_delta = 0.9))
  saveRDS(fit_prior, file = prior_file)
  par0 <- extract_parameters(fit_prior, param, param_ind, param_obs, pt, data_prior) # save for comparing prior to posterior
  saveRDS(par0, file = par0_file)
} else {
  fit_prior <- readRDS(prior_file)
  par0 <- readRDS(par0_file)
}

# If divergent transitions
# Probably caused by the fact that prior but probability on regions of high curvature (e.g. hierarchical sigma close to 0)
# In that case shouldn't be a problem
# But check anyway if chains are not getting stuck...

if (FALSE) {
  check_hmc_diagnostics(fit_prior)
  # pairs(fit_prior, pars = param_pop)
  # plot(fit_prior, pars = param_pop, plotfun = "trace")
  
  # Distribution of parameters
  plot(fit_prior, pars = setdiff(param_pop, "beta"), plotfun = "hist")
  if (length(param_ind) > 0) {plot(fit_prior, pars = paste0(param_ind, "[1]"), plotfun = "hist")}
  if (model_name == "SSMX") {plot(fit_prior, pars = "beta[1]", plotfun = "hist")}
  
  # Predictive distribution
  lapply(sample(pt, 4),
         function(pid) {
           library(ggplot2)
           ggplot(data = subset(par0, Variable == "S_rep" & Patient == pid),
                  aes(x = Week, ymin = pmax(`5%`, 0), ymax = pmin(`95%`, score_char$Range))) +
             geom_ribbon(alpha = 0.5) +
             scale_y_continuous(limits = c(0, score_char$Range)) +
             theme_bw(base_size = 15) +
             theme(panel.grid.minor.y = element_blank())
         }) %>%
    plot_grid(plotlist = ., ncol = 2)
  
  # score-50 (e.g. EASI-50)
  apply(rstan::extract(fit_prior, pars = "S_rep")[[1]],
        1,
        function(x) {
          k <- with(data_prior, 1:N_pt)
          mult <- with(data_prior, (N_obs + N_mis) / N_pt)
          start_k <- (k - 1) * mult + 1
          end_k <- k * mult
          mean(x[end_k] < 0.5 * x[start_k])
        }) %>%
    hist(., breaks = with(data_prior, (0:N_pt) / N_pt), probability = TRUE, col = "#B2001D", main = "", xlab = "EASI-50")
  
}

# Generate fake data ---------------------------------------------------------------

# Take one draw (different draws corresponds to different a priori pattern in the data)
draw <- 35

# True parameters
true_param <- extract_parameters_from_draw(fit_prior, c(param_pop, param_ind), draw)
true_param[["Patient"]] <- pt[true_param[["Index"]]]

# Dataframe
sim <- rstan::extract(fit_prior, pars = "S_rep")[[1]][draw, ]
# sim <- round(sim) # Round for observed severity
fd <- data.frame(Patient = rep(pt, each = n_dur),
                 Week = 2 * rep(0:(n_dur - 1), n_pt),
                 S = sim)
fd$S[!(fd$Week %in% c(0, 2, 4, 8, 12, 24))] <- NA

# Plot trajectories
lapply(sample(pt, 4),
       function(pid) {
         library(ggplot2)
         ggplot(data = subset(fd, Patient == pid),
                aes(x = Week)) +
           geom_point(aes(y = S)) +
           scale_y_continuous(limits = c(0, score_char$Range)) +
           labs(title = paste("Patient", pid),
                subtitle = paste("alpha = ", signif(true_param %>% filter(Parameter == "alpha" & Patient == pid) %>% pull(Value), 2))) +
           theme_bw(base_size = 20) +
           theme(panel.grid.minor.y = element_blank())
       }) %>%
  plot_grid(plotlist = ., ncol = 2)

if (model_name == "SSMX") {
  true_param %>%
    filter(Parameter == "beta") %>%
    pull(Value) %>%
    sort() %>%
    barplot()
}

# Fit fake data ---------------------------------------------------------

data_fake <- with(fd,
                  list(
                    N_obs = sum(!is.na(S)),
                    N_mis = sum(is.na(S)),
                    N_pt = length(unique(Patient)),
                    max_score = score_char$Range,
                    idx_obs = which(!is.na(S)),
                    idx_mis = which(is.na(S)),
                    S_obs = na.omit(S),
                    
                    N_test = 0,
                    idx_test = vector(),
                    S_test = vector(),
                    
                    # For horseshoe
                    p0 = 5,
                    slab_scale = 1,
                    slab_df = 5,
                    N_cov = ncol(bio),
                    X_cov = bio,
                    parametrisation = 0,
                    
                    run = 1,
                    rep = 1
                  ))

if (run_fake) {
  fit_fake <- sampling(compiled_model,
                       data = data_fake,
                       pars = param,
                       iter = n_it,
                       chains = n_chains,
                       control = list(adapt_delta = case_when(model_name %in% c("SSM", "SSMX") ~ 0.99,
                                                              TRUE ~ 0.9)))
  saveRDS(fit_fake, file = fake_file)
} else {
  fit_fake <- readRDS(fake_file)
}

# Fake data check ----------------------------------------------------------

if (FALSE) {
  
  check_hmc_diagnostics(fit_fake)
  
  pairs(fit_fake, pars = setdiff(param_pop, "beta"))
  # pairs(fit_fake, pars = paste0("beta[", 1:5, "]"))
  # print(fit_fake, pars = param_pop)
  
  # Sensitivity to prior
  par <- extract_parameters(fit_fake, param, param_ind, param_obs, pt, data_fake)
  HuraultMisc::check_model_sensitivity(par0, par, param)
  
  ## Can we recover known parameters?
  tmp <- HuraultMisc::summary_statistics(fit_fake, param) %>%
    inner_join(true_param, by = c("Variable" = "Parameter", "Index")) %>%
    rename(True = Value)
  tmp$Patient <- factor(tmp$Patient, levels = pt)
  
  # Population parameters
  tmp %>%
    filter(Variable %in% setdiff(param_pop, "beta")) %>%
    ggplot(aes(x = Variable)) +
    geom_pointrange(aes(y = Mean, ymin = `5%`, ymax = `95%`)) +
    geom_point(aes(y = True), col = "#E69F00", size = 2) +
    coord_flip() +
    labs(x = "", y = "Estimate") +
    theme_bw(base_size = 20)
  
  # Patient parameters
  lapply(param_ind,
         function(par_name) {
           tmp %>%
             filter(Variable == par_name) %>%
             mutate(Patient = fct_reorder(Patient, True)) %>%
             ggplot(aes(x = Patient)) +
             geom_pointrange(aes(y = Mean, ymin = `5%`, ymax = `95%`)) +
             geom_point(aes(y = True), col = "#E69F00", size = 2) +
             coord_flip() +
             labs(y = par_name) +
             theme_bw(base_size = 15)
         }) %>%
    plot_grid(plotlist = ., nrow = 1)
  
  if (model_name == "SSMX") {
    # beta
    tmp %>%
      filter(Variable == "beta") %>%
      mutate(Index = fct_reorder(factor(Index), True)) %>%
      ggplot(aes(x = Index)) +
      geom_pointrange(aes(y = Mean, ymin = `5%`, ymax = `95%`)) +
      geom_point(aes(y = True), col = "#E69F00", size = 2) +
      coord_flip() +
      labs(x = "", y = "Estimate") +
      theme_bw(base_size = 20)
    
    # Coverage of beta
    HuraultMisc::plot_coverage(extract(fit_fake, pars = "beta")[[1]],
                               subset(true_param, Parameter == "beta")$Value)
  }
  
  # Coverage of the posterior predictive distribution
  yrep_fake <- rstan::extract(fit_fake, pars = "S_rep")[[1]]
  HuraultMisc::plot_coverage(yrep_fake, fd[["S"]])
  
  # Coverage of patient-dependent parameters
  lapply(param_ind,
         function(x) {
           HuraultMisc::plot_coverage(rstan::extract(fit_fake, pars = x)[[1]],
                                      subset(true_param, Parameter == x)$Value)
         }) %>%
    plot_grid(plotlist = .)
  
  ## Posterior predictive checks
  ssi <- full_join(HuraultMisc::extract_distribution(fit_fake, "S_rep", type = "hdi", CI_level = seq(0.1, 0.9, 0.1)),
                   get_index(pt, data_fake),
                   by = "Index")
  pl <- lapply(sample(pt, 4),
               function(pid) {
                 PPC_fanchart(ssi, fd %>% rename(y = "S"), pid, score_char$Range) +
                   labs(y = score, title = paste("Patient", pid))
               })
  plot_grid(get_legend(pl[[1]] + theme(legend.position = "top")),
            plot_grid(plotlist = lapply(pl,
                                        function(p) {
                                          p + theme(legend.position = "none")
                                        }),
                      nrow = 2, labels = "AUTO"),
            nrow = 2, rel_heights = c(.1, .9))
  
}
