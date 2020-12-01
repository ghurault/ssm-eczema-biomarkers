# Notes -------------------------------------------------------------------

# We implement a mixture of K-fold cross-validation (leave N patients out) and forward chaining
# 1/k patients are used for testing and 1-1/k patients for training
# The model is trained on the complete data of training patients and the data up to week w = TrainingWeek (included) for testing patients
# The model is tested on the data after week w for testing patients
# This process is repeated for different w (i.e. we provide more data for the testing patients)
# And this process is repeated for different subsets of training and testing patients

# If k = 1, there is no cross-validation, just forward chaining (predict N step ahead)
# In that case, w cannot be 0!
# If TestingWeek(i) = TrainingWeek(i + 1) then this is prediction one step ahead

# If k > 0 and w = TrainingWeek = 0, this is "predict given initial point"

# Initialisation ----------------------------------------------------------

rm(list = ls()) # Clear Workspace (but better to restart session)

library(HuraultMisc) # Functions shared across projects
library(tidyverse)
library(cowplot)
library(rstan)
rstan_options(auto_write = TRUE) # Save compiled model
options(mc.cores = parallel::detectCores()) # Parallel computing
library(foreach)
library(doParallel)
source("functions.R") # Additional functions

seed <- 462528635 # seed also used for stan
set.seed(seed)

#### OPTIONS
score <- "EASI"
model_name <- "SSM"
k <- 7 # Number of folds for k-fold cross-validation, set to 1 if you don't want k-fold
run <- FALSE
n_it <- 2000
n_chains <- 4
n_cluster <- 7
####

score <- match.arg(score, c("EASI", "SCORAD", "oSCORAD", "POEM"))
model_name <- match.arg(model_name, c("Uniform", "RW", "AR", "MixedAR", "SSM", "SSMX"))

stan_code <- file.path("Models", paste0(model_name, ".stan"))
res_file <- file.path("Results", paste0("val_", score, "_", model_name, ".rds"))
res_dir <- file.path("Results", paste0("val_", score, "_", model_name)) # temporary directory

param <- c("S_pred", "lpd")

if (run & model_name != "Uniform") {
  compiled_model <- rstan::stan_model(stan_code)
}

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

# Validation -------------------------------------------------------------------

stopifnot(k == round(k),
          k > 0 & k < length(pt))
weeks  <- c(0, 2, 4, 8, 12, 24)
if (k > 1) {
  folds <- sample(cut(1:length(pt), breaks = k, labels = FALSE)) # K-fold
  it <- expand.grid(Fold = 1:k, TrainingWeek = weeks[-length(weeks)]) # we can train from week 0 since we test only a subset of patients
} else {
  folds <- rep(1, length(pt))
  it <- expand.grid(Fold = 1:k, TrainingWeek = weeks[c(-1, -length(weeks))])
}

format_data <- function(df, score, idx_test) {
  list(
    N_obs = sum(!is.na(df[, score])),
    N_mis = sum(is.na(df[, score])),
    N_pt = length(unique(df$Patient)),
    max_score = score_char$Range,
    idx_obs = which(!is.na(df[, score])),
    idx_mis = which(is.na(df[, score])),
    S_obs = na.omit(df[, score]),
    
    N_test = length(idx_test),
    idx_test = idx_test,
    S_test = df[idx_test, "S"],
    
    # For horseshoe
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

if (run) {
  
  duration <- Sys.time()
  cl <- makeCluster(n_cluster)
  registerDoParallel(cl)
  
  writeLines(c(""), "log.txt")
  dir.create(res_dir)
  
  out <- foreach(i = 1:nrow(it)) %dopar% {
    w <- it$TrainingWeek[i]
    f <- it$Fold[i]
    
    library(tidyverse)
    library(rstan)
    rstan_options(auto_write = TRUE) # Save compiled model
    options(mc.cores = parallel::detectCores()) # Parallel computing
    source("functions.R")
    
    sink("log.txt", append = TRUE)
    cat(paste0("Starting training at week ", w, ", fold ", f, " \n"))
    
    ###
    
    ## Prepare data
    dt_wf <- dt
    dt_wf$S <- dt_wf[, score]
    dt_wf[, c("SCORAD", "oSCORAD", "EASI", "POEM", "ITCH", "SLEEP")] <- NULL
    
    dt_wf$S_train <- dt_wf$S
    idx_pred <- which((dt_wf$Week %in% weeks) & !is.na(dt_wf$S) & (dt_wf$Week > w) & (dt_wf$Patient %in% pt[folds == f]))
    dt_wf$S_train[idx_pred] <- NA
    
    data_stan <- format_data(dt_wf, "S_train", idx_pred)
    
    perf <- data.frame(Patient = dt_wf$Patient[idx_pred],
                       TrainingWeek = w,
                       TestingWeek = dt_wf$Week[idx_pred],
                       Fold = f,
                       S = dt_wf$S[idx_pred])
    
    if (model_name == "Uniform") {
      
      perf <- perf %>%
        mutate(Mean_pred = score_char$Range / 2,
               lpd = -log(score_char$Range),
               CRPS = scoringRules::crps_unif(perf[["S"]], min = 0, max = score_char$Range),
               Samples = NA)
      
    } else {
      ## Fit
      fit <- sampling(compiled_model,
                      data = data_stan,
                      pars = param,
                      iter = n_it,
                      chains = n_chains,
                      seed = seed,
                      control = list(adapt_delta = case_when(model_name %in% c("SSM", "SSMX") ~ 0.99,
                                                             TRUE ~ 0.9)))
      
      ## Prepare ouput
      lpd <- extract(fit, pars = "lpd")[[1]]
      pred <- extract(fit, pars = "S_pred")[[1]]
      smp <- sapply(1:ncol(pred), function(i) {list(pred[, i])})
      
      perf <- perf %>%
        mutate(Mean_pred = apply(pred, 2, mean), # cf. point prediction (mean)
               lpd = apply(lpd, 2, function(x) {log(mean(exp(x)))}), # marginalise lpd
               CRPS = scoringRules::crps_sample(perf[["S"]], t(pred)),
               Samples = smp)
    }

    ## Save (intermediate results)
    saveRDS(perf, file = file.path(res_dir, paste0("val_", i, ".rds")))
    
    cat(paste0("Ending training at week ", w, ", fold ", f, " \n"))
  }
  stopCluster(cl)
  (duration = Sys.time() - duration)
  
  # Recombine results
  files <- list.files(res_dir)
  if (length(files) < nrow(it)) {
    warning("Number of files (", length(files), ") less than the number of iterations (", max_it + 1, "). Some runs may have failed.")
  }
  res <- do.call(rbind,
                 lapply(files,
                        function(f) {
                          readRDS(file.path(res_dir, f))
                        }))
  saveRDS(res, file = res_file)
} else {
  res <- readRDS(res_file)
}
