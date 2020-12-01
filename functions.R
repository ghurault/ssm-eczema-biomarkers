cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Data --------------------------------------------------------------------

load_dataset <- function() {
  # Load and prepare systemic therapy dataset
  #
  # - Select week 0 biomarkers
  # - Remove missing biomarkers
  # - Remove missing patients
  # - Log and standardize features
  # - Impute missing biomarkers and demographics
  # - Add row for missing weeks in severity dataframe
  #
  # Returns:
  # List containing the patient and the severity time-series dataframe
  
  library(TanakaData) # # Contains data and data processing functions
  library(dplyr)
  
  # Process biomarkers
  bio <- biomarkers_SystemicTherapy %>%
    # Select biomarkers at week 0
    filter(Week == 0) %>%
    select(-Week) %>%
    # Remove almost missing biomarkers
    select(-IL1a, -GCSF) %>%
    # Log and standardize biomarkers
    mutate(across(-matches("Patient"), ~scale(log10(.x))))
  
  # Process demographics
  demo <- patient_SystemicTherapy %>%
    # Transform to numeric binary variables and standardize Age
    mutate(FLG = as.numeric(FLG),
           Sex = as.numeric(Sex == "M"), # 0 female, 1 male
           Treatment = as.numeric(Treatment == "AZA"), # 0 MTX, 1 AZA
           Age = scale(Age))
  
  # Patient dataframe
  dp <- full_join(demo, bio, by = "Patient")
  # Severity dataframe
  dt <- severity_SystemicTherapy
  
  pt <- intersect(unique(dp[["Patient"]]), unique(dt[["Patient"]]))
  # Exclude patient 140 (mostly missing)
  pt <- pt[pt != 140]
  dp <- dp %>% filter(Patient %in% pt)
  dt <- dt %>% filter(Patient %in% pt)
  
  # Impute missing in dp by 0 (mean/default value for binary data)
  dp <- dp %>%
    mutate(across(-matches("Patient"), ~tidyr::replace_na(., 0)))
  
  # Add rows for missing severity
  dt <- bind_rows(dt,
                  setdiff(expand_grid(Patient = pt,
                                      Week = seq(0, 24, 2)),
                          dt %>% select(Patient, Week)))
  stopifnot(nrow(dt) == length(pt) * 13)
  
  # Reorder dataframes
  dp <- dp %>% arrange(Patient)
  dt <- dt %>% arrange(Patient, Week)
  
  return(list(patient_data = dp, severity_data = dt))
}

# Fitting --------------------------------------------------------------

plot_coef <- function(fit, parName, parLabel = NULL, CI = c(.05, .95), limits = NULL, ylab = "") {
  # Plot patient coefficient estimates from stan model (custom function) 
  #
  # Args:
  # fit: stanfit object
  # parName: name of the patient-dependent parameter in fit
  # parLabel: vector of names for parName
  # CI: vector of length two indicating the credible interval lower and upper bounds
  # limits: vector of length two indicating the range of estimates to plot
  #
  # Returns:
  # Ggplot of patient coefficient estimates
  
  library(ggplot2)
  
  tmp <- rstan::extract(fit, pars = parName)[[1]]
  
  if (is.null(parLabel)) {parLabel <- paste(parName, 1:ncol(tmp) ,sep= "_")}
  
  d <- data.frame(Parameter = factor(parLabel, levels = parLabel, labels = parLabel),
                  Mean = apply(tmp, 2, mean),
                  Lower = apply(tmp, 2, function(x) {quantile(x, probs = min(CI))}),
                  Upper = apply(tmp, 2, function(x) {quantile(x, probs = max(CI))}))
  d$Parameter = factor(d$Parameter, levels = rev(parLabel))
  
  p <- ggplot(data = d, aes(x = Parameter, y = Mean, ymin = Lower, ymax = Upper)) +
    geom_pointrange() +
    labs(y = parName, x = ylab) +
    coord_flip(ylim = limits) +
    theme_bw(base_size = 20) + theme(panel.grid.minor.x = element_blank())
  
  return(p)
}

get_index <- function(pt, data_stan) {
  # Associate (Patient, Week) pairs to corresponding index in the model
  #
  # Args:
  # pt: vector of patients ID (same order as the patient parameters in stanfit)
  # data_stan: data input to the stan function
  #
  # Returns:
  # Dataframe
  
  mult <- with(data_stan, (N_obs +  N_mis)/N_pt)
  max_week <- (mult - 1) * 2
  out <- data.frame(Patient = pt[rep(1:data_stan$N_pt, each = mult)],
                    Week = rep(seq(0, max_week, 2), length(pt)))
  out$Patient <- as.character(out$Patient)
  out$Index <- 1:nrow(out)
  return(out)
}

extract_parameters <- function(fit, param, param_ind, param_obs, pt, data_stan) {
  # Extract parameters' summary
  #
  # Args:
  # fit: stanfit object
  # param: parameters to extract
  # param_ind: individual parameters in param
  # param_obs: observation parameters in param
  # pt: vector of patients ID (same order as the patient parameters in stanfit)
  # data_stan: data input to the stan function
  #
  # Returns:
  # Dataframe containing posterior summary statistics of the parameters 
  
  par <- HuraultMisc::summary_statistics(fit, param)
  par$Patient <- NA
  par$Week <- NA
  
  pt <- as.character(pt)
  
  ## Patient-dependent parameter
  for (i in intersect(param_ind, param)) {
    idx <- which(par$Variable == i)
    par$Patient[idx] <- pt[par$Index[idx]]
  }
  
  ## Patient and time-dependent parameter (observation parameters)
  dict <- get_index(pt, data_stan)
  for (i in intersect(param_obs, param)) {
    idx <- sort(which(par$Variable == i))
    par[idx, c("Patient", "Week")] <- dict[, c("Patient", "Week")]
  }
  
  ## Missing score
  for (i in intersect("S_mis", param)) {
    idx <- which(par$Variable == i)
    id_mis <- data_stan$idx_mis[par$Index[idx]]
    par[idx, c("Patient", "Week")] <- dict[id_mis, c("Patient", "Week")]
  }
  
  # par$Index <- NULL
  return(par)
}

PPC_fanchart <- function(ssi, df = NULL, patient_id, max_score = NULL) {
  # PPC plot with stacked prediction intervals (fan chart) centered around the median
  #
  # Args:
  # ssi: Dataframe summarising predictive distribution as credible intervals (with columns: Lower, Upper, Level, Patient, Day)
  # df: dataframe of observed trajectory (can be NULL, in that case the actual trajectory is not overlapped)
  # patient_id: patient ID
  # max_score: maximum value that the measure can take (for plotting)
  #
  # Returns:
  # Ggplot
  
  library(ggplot2)
  
  stopifnot(is.data.frame(ssi),
            all(c("Lower", "Upper", "Level", "Patient", "Week") %in% colnames(ssi)),
            patient_id %in% unique(ssi[["Patient"]]),
            is.null(df) || is.data.frame(df))
  
  if (is.data.frame(df)) {
    stopifnot(all(c("Patient", "Week", "y") %in% colnames(df)),
              patient_id %in% unique(df[["Patient"]]))
  }
  if (!is.null(max_score)) {
    stopifnot(is.numeric(max_score),
              max_score > 0)
  }
  
  lvl <- sort(unique(ssi[["Level"]]), decreasing = TRUE)
  
  p <- ggplot()
  # Prediction intervals (cf. fill cannot be an aesthetic with a ribbon)
  for (i in 1:length(lvl)) {
    p <- p + geom_ribbon(data = subset(ssi, Patient == patient_id & Level == lvl[i]),
                         aes(x = Week, ymin = Lower, ymax = Upper, fill = Level))
  }
  # Actual trajectory
  if (is.data.frame(df)) {
    sub_df <- subset(df, Patient == patient_id)
    if ("Validation" %in% colnames(df)) {
      p <- p +
        geom_point(data = sub_df, aes(x = Week, y = y, colour = Validation), size = 1)
    } else {
      p <- p +
        geom_point(data = sub_df, aes(x = Week, y = y), size = 1)
    }
  }
  # Formatting
  p <- p +
    scale_x_continuous(expand = expansion(mult = .01)) +
    scale_fill_gradientn(colours = rev(c("#FFFFFF", RColorBrewer::brewer.pal(n = 6, "Blues")))[-1],
                         limits = c(0, 1), breaks = c(.1, .5, .9)) + # seq(0, 1, 0.25)
    scale_colour_manual(values = c("#000000", "#E69F00")) +
    labs(fill = "Confidence level", colour = "") +
    theme_classic(base_size = 15)
  
  if (!is.null(max_score)) {
    p <- p +
      scale_y_continuous(limits = c(0, max_score),
                         breaks = c(round(seq(0, max_score, length.out = 5), -1)[-5], max_score),
                         expand = c(0, 0.01 * max_score))
  }
  
  return(p)
}
