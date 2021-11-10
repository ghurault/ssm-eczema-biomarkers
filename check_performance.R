# Notes -------------------------------------------------------------------

# Analyse the predictive performance of the different models
# - Performance one-step-ahead
# - Raw performance estimates
# - Learning curve from meta-model (controlling for prediction horizon)

# Initialisation ----------------------------------------------------------

rm(list = ls()) # Clear Workspace (but better to restart session)

library(HuraultMisc) # Functions shared across projects
library(tidyverse)
library(cowplot)
source("functions.R") # Additional functions

#### OPTIONS
score <- "EASI"
metric <- "lpd"
model_names <- c("Uniform", "RW", "AR", "MixedAR", "SSM")  # "SSMX"
####

score <- match.arg(score, c("EASI", "SCORAD", "oSCORAD", "POEM"))
metric <- match.arg(metric, c("lpd", "CRPS", "ProbAccuracy", "Accuracy", "RMSE"))
stopifnot(all(model_names %in% c("Uniform", "RW", "AR", "MixedAR", "SSM", "SSMX")))
res_files <- file.path("Results", paste0("val_", score, "_", model_names, ".rds"))
stopifnot(all(file.exists(res_files)))

score_char <- data.frame(Score = c("SCORAD", "oSCORAD", "EASI", "POEM"),
                         Range = c(103, 83, 72, 28),
                         MCID = c(8.7, 8.2, 6.6, 3.4)) %>%
  filter(Score == score)

# Process results ---------------------------------------------------------

perf <- do.call(bind_rows,
                lapply(1:length(model_names),
                       function(i) {
                         res <- readRDS(res_files[i])
                         
                         # Probabilistic accuracy
                         if (model_names[i] == "Uniform") {
                           ub <- pmin(score_char$Range, res$S + score_char$MCID)
                           lb <- pmax(0, res$S - score_char$MCID)
                           acc <- (ub - lb) / score_char$Range
                         } else {
                           acc <- sapply(1:nrow(res), function(j) {
                             mean(abs(res$S[j] - res$Samples[j][[1]]) < score_char$MCID)
                           })
                         }

                         res %>%
                           mutate(SquaredError = (S - Mean_pred)^2,
                                  ProbAccuracy = acc,
                                  Accuracy = as.numeric(abs(S - Mean_pred) < score_char$MCID)) %>%
                           mutate(Model = model_names[i]) %>%
                           select(-Samples)
                       })) %>%
  mutate(Model = factor(Model, levels = rev(model_names)))

# One-steap-ahead performance ---------------------------------------------
# Prediction for the next clinical visits
# Prediction horizon differ though

# Select one-step-ahead prediction
cv_osa <- perf %>%
  group_by(Model, TrainingWeek) %>%
  filter(TestingWeek == min(TestingWeek)) %>%
  ungroup()
# Compute performance for each fold
cv_osa <- cv_osa %>%
  group_by(Model, Fold) %>%
  summarise(lpd = mean(lpd),
            CRPS = mean(CRPS),
            ProbAccuracy = mean(ProbAccuracy),
            Accuracy = mean(Accuracy),
            RMSE = sqrt(mean(SquaredError))) %>%
  ungroup()
# Average performance across fold
cv_osa <- cv_osa %>%
  pivot_longer(cols = all_of(c("lpd", "CRPS", "ProbAccuracy", "Accuracy", "RMSE")), names_to = "Metric", values_to = "Value") %>%
  group_by(Model, Metric) %>%
  summarise(Mean = mean(Value), SE = sd(Value) / sqrt(n()))

p1 <- cv_osa %>%
  filter(Metric == metric) %>%
  ggplot(aes(x = Model, y = Mean, ymin = Mean - SE, ymax = Mean + SE)) +
  geom_pointrange() +
  coord_flip() +
  labs(x = "", y = metric) +
  theme_bw(base_size = 15)
if (metric == "Accuracy") {
  p1 <- p1 + scale_y_continuous(limits = c(0, 1))
}
if (metric %in% c("CRPS", "RMSE")) {
  p1 <- p1 + scale_y_continuous(limits = c(0, NA))
}
p1

# Raw performance estimates --------------------------------------------------------

# Compute performance for each fold (and each condition)
cv <- perf %>%
  group_by(Model, TrainingWeek, TestingWeek, Fold) %>%
  summarise(lpd = mean(lpd),
            CRPS = mean(CRPS),
            ProbAccuracy = mean(ProbAccuracy),
            Accuracy = mean(Accuracy),
            RMSE = sqrt(mean(SquaredError))) %>%
  ungroup()
# Average performance across fold
cv <- cv %>%
  pivot_longer(cols = all_of(c("lpd", "CRPS", "ProbAccuracy", "Accuracy", "RMSE")), names_to = "Metric", values_to = "Value") %>%
  group_by(Model, TrainingWeek, TestingWeek, Metric) %>%
  summarise(Mean = mean(Value), SE = sd(Value) / sqrt(n())) %>%
  ungroup()
# Compute prediction horizon
cv <- cv %>%
  mutate(Horizon = TestingWeek - TrainingWeek)

# Performance as a function of prediction horizon, for each training week (Model in colour)
# Alternatively plot as a function of prediction horizon, for each model (training week in colour); harder for model comparison
p2 <- cv %>%
  filter(Metric == metric) %>%
  ggplot(aes(x = Horizon, y = Mean, ymin = Mean - SE, ymax = Mean + SE, colour = factor(Model))) +
  facet_grid(rows = vars(TrainingWeek)) +
  # ggplot(aes(x = Horizon, y = Mean, ymin = Mean - SE, ymax = Mean + SE, colour = factor(TrainingWeek))) +
  # facet_grid(cols = vars(Model)) +
  geom_pointrange() +
  geom_line() +
  scale_color_manual(values = cbbPalette) +
  labs(x = "Prediction Horizon (weeks)",
       y = metric,
       colour = "") +
  scale_x_continuous(breaks = sort(unique(cv[["Horizon"]]))) +
  theme_bw(base_size = 15) +
  theme(panel.grid.minor.x = element_blank())
if (metric == "Accuracy") {
  p2 <- p2 + scale_y_continuous(limits = c(0, 1))
}
if (metric %in% c("CRPS", "RMSE")) {
  p2 <- p2 + scale_y_continuous(limits = c(0, NA))
}
p2

if (FALSE) {
  ggsave(file.path("Plots", paste0(score, "_", metric, "_rawperf.jpg")),
         width = 10, height = 10, units = "cm", dpi = 300, scale = 2)
}

# Learning curves from meta-model -----------------------------------------

cv <- perf %>%
  group_by(Model, TrainingWeek, TestingWeek, Fold) %>%
  summarise(lpd = mean(lpd),
            CRPS = mean(CRPS),
            ProbAccuracy = mean(ProbAccuracy),
            Accuracy = mean(Accuracy),
            RMSE = sqrt(mean(SquaredError))) %>%
  ungroup() %>%
  mutate(Horizon = TestingWeek - TrainingWeek)

estimate_performance <- function(df, metric, adjust_horizon = TRUE) {
  # Estimate learning curves with a meta-model (linear regression)
  #
  # Args:
  # df: Dataframe of performance metric per fold
  # metric: Metric name
  # adjust_horizon:  Whether to adjust for prediciton horizon in the model
  #
  # Returns:
  # Dataframe with columns: TrainingWeek, Horizon, Mean, SE, Variable
  
  stopifnot(is.data.frame(df),
            is.character(metric),
            all(c("TrainingWeek", "Horizon", "Fold", metric) %in% colnames(df)),
            is.logical(adjust_horizon))
  
  f <- paste0(metric, " ~ factor(TrainingWeek) + 0")
  if (adjust_horizon) {
    f <- paste0(f, " + Horizon")
  }
  f <- formula(f)
  
  meta_model <- glm(f,
                    family = "gaussian",
                    data = df)
  lm_fit <- data.frame(TrainingWeek = c(0, 2, 4, 8, 12), Horizon = 2)
  pred <- predict(meta_model, newdata = lm_fit, se.fit = TRUE)
  lm_fit <- lm_fit %>%
    mutate(Mean = pred$fit,
           SE = pred$se.fit,
           Variable = "Fit")
  
  if (adjust_horizon) {
    s <- summary(meta_model)
    lm_horizon <- data.frame(Mean = s$coefficients["Horizon", "Estimate"],
                             SE = s$coefficients["Horizon", "Std. Error"],
                             Variable = "Horizon")
  } else {
    lm_horizon <- data.frame(Mean = 0,
                             SE = 0,
                             Variable = "Horizon")
  }
  
  bind_rows(lm_fit, lm_horizon) %>%
    mutate(Metric = metric)
}

fit_perf <- do.call(rbind,
                    lapply(model_names,
                           function(x) {
                             cv %>%
                               filter(Model == x) %>%
                               estimate_performance(., metric, adjust_horizon = !((x == "Uniform") & (metric == "lpd"))) %>%
                               mutate(Model = x)
                           })) %>%
  mutate(Model = factor(Model, levels = rev(model_names)))
  
p3 <- fit_perf %>%
  filter(Variable == "Fit") %>%
  ggplot(aes(x = TrainingWeek, y = Mean, ymin = Mean - SE, ymax = Mean + SE, colour = Model, fill = Model)) +
  geom_line() +
  # geom_pointrange(position = position_dodge(width = 1)) +
  geom_point() +
  geom_ribbon(alpha = 0.5) +
  scale_colour_manual(values = cbbPalette) +
  scale_fill_manual(values = cbbPalette) +
  scale_x_continuous(breaks = sort(unique(fit_perf[["TrainingWeek"]]))) +
  labs(x = "Training week", y = metric, colour = "", fill = "") +
  theme_bw(base_size = 15) +
  theme(panel.grid.minor.x = element_blank())

p4 <- fit_perf %>%
  filter(Variable == "Horizon") %>%
  ggplot(aes(x = Model, y = Mean, ymin = Mean - SE, ymax = Mean + SE, colour = Model)) +
  geom_pointrange(size = 1.5) +
  scale_colour_manual(values = cbbPalette) +
  labs(x = "", y = paste0(metric, " change with increasing \nprediction horizon of 2 weeks"), colour = "") +
  theme_bw(base_size = 15) +
  theme(legend.position = "bottom",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

plot_grid(p3 + theme(legend.position = "none"),
          p4 + theme(legend.position = "none"),
          get_legend(p3 + theme(legend.position = "right")),
          nrow = 1, rel_widths = c(4, 3, 1), labels = c("A", "B", ""))
if (FALSE) {
  ggsave(file.path("Plots", paste0(score, "_", "metric", "_metaperf.jpg")),
         width = 21, height = 13, units = "cm", dpi = 400, scale = 1.5, bg = "white")
}
