// Autoregressive model of order 1 with patient-dependent parameters

functions {
  real normal_lb_ub_rng(real mu, real sigma, real lb, real ub) {
    // Sample from truncated normal of mean mu, standard deviation sigma, lower bound lb and upper bound ub
    real u;
    real p1;
    real p2;
    if (is_nan(mu) || is_inf(mu)) {
      reject("normal_lb_ub_rng: mu must be finite; ", "found mu = ", mu);
    } else if (is_nan(sigma) || is_inf(sigma) || sigma < 0) {
      reject("normal_lb_ub_rng: sigma must be finite and non-negative; ", "found sigma = ", sigma);
    } else if (lb >= ub) {
      reject("normal_lb_ub_rng: lb must be less than ub; ", "found lb = ", lb, "and ub = ", ub);
    } else {
      p1 = normal_cdf(lb, mu, sigma);
      p2 = normal_cdf(ub, mu, sigma);
      if (p1 >= p2) {
        reject("normal_lb_ub_rng: p1 >= p2. p1 = ", p1, " and p2 = ", p2, ". mu = ", mu, "and sigma = ", sigma);
      } else {
        u = uniform_rng(p1, p2);
      }
    }
    return (sigma * inv_Phi(u)) + mu;  // inverse cdf 
  }
}

data {
  int<lower = 0> N_obs; // Number of non-missing observations
  int<lower = 0> N_mis; // Number of missing observations
  int<lower = 0> N_pt; // Number of patients
  real<lower = 0> max_score; // Maximum value that the score can take
  int<lower = 1, upper = N_obs + N_mis> idx_obs[N_obs]; // Index of non-missing observations
  int<lower = 1, upper = N_obs + N_mis> idx_mis[N_mis]; // Index of missing observations
  real<lower = 0, upper = max_score> S_obs[N_obs]; // Observed score
  
  int<lower = 0, upper = N_mis> N_test; // Number of predictions to evaluate (predictions are set to missing in the data)
  int<lower = 1, upper = N_obs + N_mis> idx_test[N_test]; // Index of test observations
  real<lower = 0, upper = max_score> S_test[N_test]; // Observed score for predictions
  
  int<lower = 0, upper = 1> run; // Switch to evaluate the likelihood
  int<lower = 0, upper = 1> rep; // Switch to generate replications
}

transformed data {
  int N = N_obs + N_mis; // Total number of observations
  int mult = N / N_pt; // Number of timepoints per patient (same by construction)
  int start[N_pt]; // index of first observation for patient each patient
  int end[N_pt]; // index of last observation for patient each patient
  
  for (k in 1:N_pt) {
    start[k] = (k - 1) * mult + 1;
    end[k] = k * mult;
  }
}

parameters {
  real<lower = 0, upper = max_score> S_mis[N_mis]; // Missing S (useful for predictions (cf. bounds))
  real<lower = 0> sigma; // Standard deviation
  
  // Autocorrelation parameter
  real<lower = 0, upper = 1> alpha[N_pt]; // Autocorrelation parameter
  real<lower = 0, upper = 1> mu_alpha; // Population mean of alpha
  real<lower = 0> phi_alpha; // Population pseudo "sample size" of alpha
  
  // Population autoregression mean (for S_lat)
  real mu_inf; // Population mean
  real<lower = 0> sigma_inf; // Population std
  real eta_inf[N_pt]; // Error term
}

transformed parameters {
  real S_meas[N]; // Measured S
  real S_inf[N_pt];
  real b[N_pt];
  
  for (k in 1:N_pt) {
    S_inf[k] = mu_inf + sigma_inf * eta_inf[k];
    b[k] = S_inf[k] * (1 - alpha[k]);
  }
  
  S_meas[idx_obs] = S_obs;
  S_meas[idx_mis] = S_mis;
}

model {
  // implicit (bounded) uniform prior for S_mis
  sigma / max_score ~ lognormal(-log(20), 0.5 * log(5)); // 95% CI is [.01, 0.25] * max_score
  
  eta_inf ~ std_normal();
  mu_inf / max_score ~ normal(0.5, 0.25);
  sigma_inf / max_score ~ normal(0, 0.125);
  
  mu_alpha ~ beta(2, 2);
  phi_alpha ~ lognormal(1 * log(10), 0.5 * log(10)); // mass between 1 and 100
  alpha ~ beta(mu_alpha * phi_alpha, (1 - mu_alpha) * phi_alpha);
  
  for (k in 1:N_pt) {
    to_vector(S_meas[(start[k] + 1):end[k]]) ~ normal(alpha[k] * to_vector(S_meas[start[k]:(end[k] - 1)]) + b[k], sigma); // Vectorise for efficiency
  }
  
}

generated quantities {
  real S_rep[N]; // Replications of S_meas[t + 1] given S_lat[t]
  real S_pred[N_test]; // Predictive sample of S_test
  real lpd[N_test]; // Log predictive density
  
  // Replications
  if (rep == 1) {
    for (k in 1:N_pt) {
      S_rep[start[k]] = S_meas[start[k]];
      for (t in (start[k] + 1):end[k]) {
        S_rep[t] = normal_lb_ub_rng(alpha[k] * S_meas[t - 1] + b[k], sigma, 0, max_score);
      } 
    }
    S_pred = S_rep[idx_test];
  }
  
  // Log predictive density
  {
    int k_test; // Patient index
    real linpred; // Linear predictor
    real Z; // Normalisation constant
    for (i in 1:N_test) {
      k_test =  ((idx_test[i] - 1) % mult) + 1;
      linpred = alpha[k_test] * S_meas[idx_test[i]] + b[k_test];
      Z = normal_cdf(max_score, linpred, sigma) - normal_cdf(0, linpred, sigma);
      lpd[i] = normal_lpdf(S_test[i] | linpred, sigma) - log(Z);
    }
  }
  
}
