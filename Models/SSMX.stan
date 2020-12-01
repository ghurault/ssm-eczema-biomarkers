// State space model with eXogeneous variables (covariates) following a regularised horseshoe prior
// 2 parametrisations are implemented for the Horseshoe
//
// The prior for tau (global shrinkage) is a function of the expected number of "important" features...
// ... and a function of the number of observations (here N_pt, cf. scale_global)
// The slab prior (i.e. for non-zero coefficients) is a scaled-inverse chi-squared distribution...
// ... where the tail is similar to a Student t distribution with slab_df degree of freedom

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
  
  real soft_uniform_lpdf(real x, real lb, real ub) {
    return(log(inv_logit(x - lb) - inv_logit(x - ub)) - log(ub - lb));
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
  
  real<lower = 0> p0; // Horseshoe: guess on the number of non-zero parameters
  real<lower = 0> slab_scale; // Horseshoe: slab scale
  real<lower = 0> slab_df; // Horseshoe: slab degrees of freedom (1 for cauchy, more for closer to gaussian)
  int<lower = 0> N_cov; // Number of covariates at t0
  matrix[N_pt, N_cov] X_cov; // Matrix of covariates at t0
  
  int<lower = 0, upper = 1> parametrisation; // Switch to change parametrisation of the horseshoe
  int<lower = 0, upper = 1> run; // Switch to evaluate the likelihood
  int<lower = 0, upper = 1> rep; // Switch to generate replications
}

transformed data {
  int N = N_obs + N_mis; // Total number of observations
  int mult = N / N_pt; // Number of timepoints per patient (same by construction)
  int start[N_pt]; // index of first observation for patient each patient
  int end[N_pt]; // index of last observation for patient each patient
  
  real scale_global = p0 / (N_cov - p0) / sqrt(N_pt); // Horseshoe: scale for tau
  real nu_local = 1; // Horseshoe: degree of freedom for lambdas prior (for horseshoe it's 1 (cauchy))
  real nu_global = 1; // Horseshoe: degree of freedom for tau prior (1 is cauchy)
  
  for (k in 1:N_pt) {
    start[k] = (k - 1) * mult + 1;
    end[k] = k * mult;
  }
}

parameters {
  real<lower = 0, upper = max_score> S_mis[N_mis]; // Missing S (useful for predictions (cf. bounds))
  real S_lat_eta[N]; // cf. non-centered parametrisation
  
  // Autocorrelation parameter
  real<lower = 0, upper = 1> alpha[N_pt]; // Autocorrelation parameter
  real<lower = 0, upper = 1> mu_alpha; // Population mean of alpha
  real<lower = 0> phi_alpha; // Population pseudo "sample size" of alpha
  
  // Population autoregression mean (for S_lat)
  real mu_inf; // Population mean
  real<lower = 0> sigma_inf; // Population std
  real eta_inf[N_pt]; // Error term
  
  real<lower = 0> sigma_tot; // total noise std
  real<lower = 0, upper = 1> rho2; // proportion of the stochastic noise variance in the total noise variance
  
  // Horseshoe
  vector[N_cov] z; // Horseshoe noise (non-centered)
  real<lower = 0> caux; // Horseshoe: cauchy noise for c
  // Horseshoe parametrisation 0
  real<lower = 0> tau0[1 - parametrisation]; // Horseshoe, parametrisation 0: global shrinkage parameter
  vector<lower = 0>[N_cov] lambda0[1 - parametrisation]; // Horseshoe, parametrisation 0: local shrinkage parameter
  // Horseshoe parametrisation 1
  real<lower = 0> aux1_global[parametrisation]; // Horseshoe: cf. parametrisation of tau
  real<lower = 0> aux2_global[parametrisation]; // Horseshoe: cf. parametrisation of tau
  vector<lower = 0>[N_cov] aux1_local[parametrisation]; // Horseshoe: cf. parametrisation of lambdas
  vector<lower = 0>[N_cov] aux2_local[parametrisation]; // Horseshoe: cf. parametrisation of lambdas
  
}

transformed parameters {
  real S_lat[N]; // Latent S
  real S_inf[N_pt];
  real b[N_pt];
  real sigma_meas = sqrt(rho2) * sigma_tot; // measurement noise std
  real sigma_lat = sqrt(1 - rho2) * sigma_tot; // stochastic noise std
  real MDC = 1.96 * sigma_meas; // minimum detectable change (95% level)
  
  // Horseshoe
  vector<lower = 0>[N_cov] lambda_tilde; // Horseshoe: truncated local shrinkage parameter
  real<lower = 0> c; // Horseshoe: slab scale
  vector[N_cov] beta; // Horseshoe: regularised coefficients
  vector[N_pt] f; // Horseshoe: latent function values (x*beta)
  real<lower = 0> tau; // Horseshoe: global shrinkage parameter
  vector<lower = 0>[N_cov] lambda; // Horseshoe: local shrinkage parameter
  
  for (k in 1:N_pt) {
    S_inf[k] = mu_inf + sigma_inf * eta_inf[k];
    b[k] = S_inf[k] * (1 - alpha[k]);
  }

  // Horseshoe
  if (parametrisation == 0) {
    tau = tau0[1];
    lambda = lambda0[1];
  } else {
    tau = aux1_global[1] * sqrt(aux2_global[1]) * scale_global * sigma_lat;
    lambda = aux1_local[1] .* sqrt(aux2_local[1]);
  }
  c = slab_scale * sqrt(caux);
  lambda_tilde = sqrt(c^2 * square(lambda) ./ (c^2 + tau^2 * square(lambda)));
  beta = z .* lambda_tilde * tau;
  f = X_cov * beta;
  
  for (k in 1:N_pt) {
    S_lat[start[k]] = max_score * (0.5 + 0.25 * S_lat_eta[start[k]]); // prior covering the full range of the score
    for (t in (start[k] + 1):end[k]) {
      S_lat[t] = alpha[k] * S_lat[t - 1] + b[k] + f[k] + sigma_lat * S_lat_eta[t];
    }
  }
}

model {
  S_lat_eta ~ std_normal();
  eta_inf ~ std_normal();
  mu_inf / max_score ~ normal(0.5, 0.25);
  sigma_inf / max_score ~ normal(0, 0.125);
  sigma_tot / max_score ~ lognormal(-log(20), 0.5 * log(5)); // 95% CI is [.01, 0.25] * max_score
  rho2 ~ beta(4, 2); // process noise expected to be small compared to the measurement noise
  
  mu_alpha ~ beta(2, 2);
  phi_alpha ~ lognormal(1 * log(10), 0.5 * log(10)); // mass between 1 and 100
  alpha ~ beta(mu_alpha * phi_alpha, (1 - mu_alpha) * phi_alpha);
  
  for (i in 1:N) {
    // S_lat[i] ~ soft_uniform(-1, max_score + 1);
    S_lat[i] ~ soft_uniform(-.01 * max_score, 1.01 * max_score);
  }

  if (run) {
    for (i in 1:N_obs) {
      S_obs[i] ~ normal(S_lat[idx_obs[i]], sigma_meas) T[0, max_score];
    }
  }
  
  //
  if (parametrisation == 0) {
    lambda0[1] ~ student_t(nu_local, 0, 1);
    tau0[1] ~ student_t(nu_global, 0, scale_global * sigma_lat);
  } else {
    aux1_local[1] ~ std_normal();
    aux2_local[1] ~ inv_gamma (0.5 * nu_local, 0.5 * nu_local);
    aux1_global[1] ~ std_normal();
    aux2_global[1] ~ inv_gamma (0.5 * nu_global , 0.5 * nu_global );
  }
  z ~ std_normal();
  caux ~ inv_gamma (0.5 * slab_df , 0.5 * slab_df);
}

generated quantities {
  real S_rep[N]; // Replications of S_meas[t + 1] given S_lat[t]
  real S_pred[N_test]; // Predictive sample of S_test
  real lpd[N_test]; // Log predictive density
  
  // Replications
  if (rep == 1) {
    for (i in 1:N) {
      S_rep[i] = normal_lb_ub_rng(S_lat[i], sigma_meas, 0, max_score);
    }
    S_pred = S_rep[idx_test];
  }
  
  // Log predictive density
  {
    real Z; // Normalisation constant
    for (i in 1:N_test) {
      Z = normal_cdf(max_score, S_lat[idx_test[i]], sigma_meas) - normal_cdf(0, S_lat[idx_test[i]], sigma_meas);
      lpd[i] = normal_lpdf(S_test[i] | S_lat[idx_test[i]], sigma_meas) - log(Z);
    }
  }
  
  
}
