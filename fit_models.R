# Notes -------------------------------------------------------------------

# Master file to fit models for the different scores:
# - RW: random walk model
# - AR: autoregressive model (order 1, fixed effects)
# - MixedAR: mixed effect autoregressive model (order 1)
# - SSM: Hidden Markov Model (Gaussian measurement error and mixed autoregressive model for the latent dynamic)
# - SSMX: Hidden Markov Model with eXogeneous variables, following a horseshoe prior (two parametrisations for the horseshoe)

# Initialisation ----------------------------------------------------------

rm(list = ls()) # Clear Workspace (but better to restart session)

library(HuraultMisc) # Functions shared across projects
library(tidyverse)
library(cowplot)
library(rstan)
rstan_options(auto_write = TRUE) # Save compiled model
options(mc.cores = parallel::detectCores()) # Parallel computing
source("functions.R") # Additional functions

seed <- 462528635 # seed also used for stan
set.seed(seed)

#### OPTIONS
score <- "EASI"
model_name <- "SSM"
run <- FALSE
n_it <- 2000
n_chains <- 4
####

score <- match.arg(score, c("EASI", "SCORAD", "oSCORAD", "POEM"))
model_name <- match.arg(model_name, c("RW", "AR", "MixedAR", "SSM", "SSMX"))

stan_code <- file.path("Models", paste0(model_name, ".stan"))
res_file <- file.path("Results", paste0("fit_", score, "_", model_name, ".rds"))
par_file <- file.path("Results", paste0("par_", score, "_", model_name, ".rds"))
par0_file <- file.path("Results", paste0("par0_", score, "_", model_name, ".rds"))

if (model_name == "RW") {
  param_pop <- c("sigma")
  param_ind <- c()
  param_obs <- c("S_rep") # "S_mis"
}
if (model_name == "AR") {
  param_pop <- c("sigma", "alpha", "S_inf", "b")
  param_ind <- c()
  param_obs <- c("S_rep") # "S_mis"
}
if (model_name == "MixedAR") {
  param_pop <- c("sigma",
                 "mu_alpha", "phi_alpha",
                 "mu_inf", "sigma_inf")
  param_ind <- c("alpha", "S_inf", "b")
  param_obs <- c("S_rep") # "S_mis"
}
if (model_name == "SSM") {
  param_pop <- c("sigma_tot", "rho2", "sigma_lat", "sigma_meas", "MDC",
                 "mu_alpha", "phi_alpha",
                 "mu_inf", "sigma_inf")
  param_ind <- c("alpha", "S_inf", "b")
  param_obs <- c("S_lat", "S_rep")
}
if (model_name == "SSMX") {
  param_pop <- c("sigma_tot", "rho2", "sigma_lat", "sigma_meas", "MDC",
                 "mu_alpha", "phi_alpha",
                 "mu_inf", "sigma_inf",
                 "beta")
  param_ind <- c("alpha", "S_inf", "b", "f")
  param_obs <- c("S_lat", "S_rep")
}
param <- c(param_pop, param_ind, param_obs)

score_char <- data.frame(Score = c("SCORAD", "oSCORAD", "EASI", "POEM"),
                         Range = c(103, 83, 72, 28),
                         MCID = c(8.7, 8.2, 6.6, 3.4)) %>%
  filter(Score == score)

# Data --------------------------------------------------------------------

l <- load_dataset()
dp <- l$patient_data
dt <- l$severity_data
pt <- unique(dt[["Patient"]])
bio <- as.matrix(dp[, colnames(dp) != "Patient"]) # matrix of biomarkers (including treatment, age, sex...)

# Model -------------------------------------------------------------------

format_data <- function(df, score) {
  list(
    N_obs = sum(!is.na(df[, score])),
    N_mis = sum(is.na(df[, score])),
    N_pt = length(unique(df$Patient)),
    max_score = score_char$Range,
    idx_obs = which(!is.na(df[, score])),
    idx_mis = which(is.na(df[, score])),
    S_obs = na.omit(df[, score]),
    
    N_test = 0,
    idx_test = vector(),
    S_test = vector(),
    
    # For horsehoe
    p0 = 5,
    slab_scale = 1,
    slab_df = 5,
    N_cov = ncol(bio),
    X_cov = bio,
    parametrisation = 0,
    
    run = 1,
    rep = 1
  )
}

data_stan <- dt %>%
  rename(y = score) %>%
  # mutate(y = replace(y, Week > 12, NA)) %>% # cf. remove test set
  format_data(., "y")

if (run) {
  fit <- stan(file = stan_code,
              data = data_stan,
              iter = n_it,
              chains = n_chains,
              pars = param,
              seed = seed,
              control = list(adapt_delta = case_when(model_name %in% c("SSM", "SSMX") ~ 0.99,
                                                     TRUE ~ 0.9)))
  saveRDS(fit, file = res_file)
  par <- extract_parameters(fit, param, param_ind, param_obs, pt, data_stan)
  saveRDS(par, file = par_file)
} else {
  fit <- readRDS(res_file)
  par <- readRDS(par_file)
}

par0 <- readRDS(par0_file)

# Diagnostics and fit ----------------------------------------------------------------

if (FALSE) {
  
  # shinystan::launch_shinystan(fit)
  check_hmc_diagnostics(fit)
  # max(par[["Rhat"]], na.rm = TRUE)
  
  pairs(fit, pars = setdiff(param_pop, "beta"))
  # pairs(fit, pars = paste0("beta[", 1:5, "]"))
  plot(fit, pars = setdiff(param_pop, "beta"), plotfun = "trace")
  plot(fit, pars = setdiff(param_pop, "beta"), plotfun = "hist")
  
  if (model_name %in% c("MixedAR", "SSM", "SSMX")) {
    # plot(fit, pars = "alpha")
    plot_grid(
      plot_coef(fit, "alpha", pt, limits = c(0, 1), ylab = "Patient"),
      plot_coef(fit, "b", pt, ylab = "Patient"),
      plot_coef(fit, "S_inf", pt, ylab = "Patient"),
      nrow = 1
    )
  }
  if (model_name == "SSMX") {
    plot_grid(
      plot_coef(fit, "beta", colnames(bio), CI = c(.05, 0.95), limits = c(-1, 1)),
      plot_coef(fit, "f", pt, CI = c(0.05, 0.95), limits = c(-1, 1)) +
        labs(x = "Patient", y = "x * beta") +
        theme(axis.text.y = element_blank()),
      labels = "AUTO", nrow = 1, rel_widths = c(.55, .45)
    )
    if (FALSE) {
      ggsave(file.path("Plots", paste0(score, "_covariates.jpg")),
             width = 10, height = 10, units = "cm", dpi = 300, scale = 2)
    }
    
  }
  # print(fit, pars = param_pop)
  
  # Check priors
  param01 <- intersect(param_pop, c("alpha", "mu_alpha", "rho2")) # parameters in 0-1
  HuraultMisc::plot_prior_posterior(par0, par, setdiff(param_pop, param01))
  if (length(param01) > 0) {
    # cf. 0-1 scale
    HuraultMisc::plot_prior_posterior(par0, par, param01) +
      coord_flip(ylim = c(0, 1)) +
      theme(legend.position = "none")
  }
  plot_prior_influence(par0, par, c(param_pop, param_ind))
  # compute_prior_influence(par0, par, param_pop)

  lapply(param_ind, function(x) {PPC_group_distribution(fit, x, 100)}) %>%
    plot_grid(plotlist = .)
  
}

# PPC Trajectories ------------------------------------------------------------

if (FALSE) {
  
  ssi <- full_join(HuraultMisc::extract_distribution(fit, "S_rep", type = "hdi", CI_level = seq(0.1, 0.9, 0.1)),
                   get_index(pt, data_stan),
                   by = "Index")
  pl <- lapply(c(108, 119, 134, 137), # sort(sample(pt, 4, replace = FALSE)), 
               function(pid) {
                 tmp <- dt %>%
                   rename(y = score)
                 if (FALSE) {
                   # Identify training and testing data when the fit is not on the full dataset
                   tmp <- tmp %>%
                     mutate(Validation = case_when(Week <= 12 ~ "Training",
                                                   Week > 12 ~ "Testing"),
                            Validation = fct_relevel(Validation, "Training", "Testing"))
                 }
                 
                 PPC_fanchart(ssi, tmp, pid, score_char$Range) +
                   labs(y = score) # , title = paste("Patient", pid))
               })
  plot_grid(get_legend(pl[[1]] + theme(legend.position = "top")),
            plot_grid(plotlist = lapply(pl,
                                        function(p) {
                                          p + theme(legend.position = "none")
                                        }),
                      nrow = 2, labels = "AUTO"),
            nrow = 2, rel_heights = c(.1, .9))
  # ggsave(file.path("Plots", paste0(score, "_PPC.jpg")), width = 30, height = 20, units = "cm", dpi = 400, bg = "white")
  
}
