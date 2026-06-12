// hierarchical_loglogistic4.stan
//
// Hierarchical 4-parameter log-logistic (Hill equation) model.
// curveRcore convention:
//
//   y = a + (d - a) / (1 + (c / x)^b)
//
// Requires x > 0 (raw concentration, NOT log-transformed).
// Symmetric on the log(x) axis. Classic dose-response / Hill model.
//
// This model is primarily useful when is_log_independent = FALSE.
// When the concentration axis IS log10-transformed, loglogistic4 is
// mathematically equivalent to logistic4 and resolve_effective_models()
// prunes it. But it must exist for raw-scale calibration curves
// (flow cytometry, some ELISA configurations, etc.).
//
// Noise model (controlled by use_heteroscedastic_noise):
//   0 = homoscedastic:     sigma_i = sigma_obs   (constant)
//   1 = heteroscedastic:   sigma_i = exp(log_sigma0 + log_sigma_slope * log(|mu_i|))
//                          (power-of-mean; CDAN as in O'Malley 2008)

data {
  int<lower=1> N_obs;
  int<lower=1> N_plates;
  int<lower=1> grainsize;
  array[N_obs] int<lower=1, upper=N_plates> plate_idx;
  vector<lower=0>[N_obs] x;             // raw concentration (must be positive)
  array[N_obs] real y;

  real prior_a_mu;
  real<lower=0> prior_a_sigma;
  real prior_d_mu;
  real<lower=0> prior_d_sigma;
  real prior_log_b_mu;                   // b on log scale for positivity
  real<lower=0> prior_log_b_sigma;
  real prior_log_c_mu;                   // c (EC50) on log scale
  real<lower=0> prior_log_c_sigma;

  int<lower=0> N_blanks;
  array[N_blanks] int<lower=1, upper=N_plates> blank_plate_idx;
  vector[N_blanks] blank_response;

  // Noise model switch: 0 = homoscedastic, 1 = heteroscedastic (CDAN)
  int<lower=0, upper=1> use_heteroscedastic_noise;
  // Heteroscedastic priors (used only when use_heteroscedastic_noise = 1)
  real prior_log_sigma0_mu;
  real<lower=0> prior_log_sigma0_sigma;
  real prior_log_sigma_slope_mu;
  real<lower=0> prior_log_sigma_slope_sigma;
}

parameters {
  real mu_a;
  real<lower=0> sigma_a;
  real mu_d;
  real<lower=0> sigma_d;
  real mu_log_b;
  real<lower=0> sigma_log_b;
  real mu_log_c;                         // c estimated on log scale (EC50 > 0)
  real<lower=0> sigma_log_c;

  vector[N_plates] raw_a;
  vector[N_plates] raw_d;
  vector[N_plates] raw_log_b;
  vector[N_plates] raw_log_c;

  real<lower=0> sigma_obs;
  real<lower=2> nu;
  real<lower=0> sigma_blank;
  // Heteroscedastic noise parameters — always estimated, used in
  // likelihood only when use_heteroscedastic_noise = 1.
  real log_sigma0;
  real log_sigma_slope;
}

transformed parameters {
  vector[N_plates] a = mu_a + sigma_a * raw_a;
  vector[N_plates] d = mu_d + sigma_d * raw_d;
  vector[N_plates] log_b = mu_log_b + sigma_log_b * raw_log_b;
  vector[N_plates] b = exp(log_b);
  vector[N_plates] log_c = mu_log_c + sigma_log_c * raw_log_c;
  vector[N_plates] c_par = exp(log_c);  // EC50 > 0
}

model {
  mu_a ~ normal(prior_a_mu, prior_a_sigma);
  sigma_a ~ normal(0, prior_a_sigma * 0.5);
  mu_d ~ normal(prior_d_mu, prior_d_sigma);
  sigma_d ~ normal(0, prior_d_sigma * 0.5);
  mu_log_b ~ normal(prior_log_b_mu, prior_log_b_sigma);
  sigma_log_b ~ normal(0, 0.5);
  mu_log_c ~ normal(prior_log_c_mu, prior_log_c_sigma);
  sigma_log_c ~ normal(0, 1.0);

  raw_a ~ std_normal();
  raw_d ~ std_normal();
  raw_log_b ~ std_normal();
  raw_log_c ~ std_normal();

  sigma_obs ~ normal(0, prior_a_sigma);
  nu ~ gamma(2, 0.1);
  sigma_blank ~ normal(0, prior_a_sigma);

  // Heteroscedastic noise parameters — always given priors.
  log_sigma0       ~ normal(prior_log_sigma0_mu,     prior_log_sigma0_sigma);
  log_sigma_slope  ~ normal(prior_log_sigma_slope_mu, prior_log_sigma_slope_sigma);

  // curveRcore loglogistic4: y = a + (d - a) / (1 + (c / x)^b)
  for (i in 1:N_obs) {
    int p = plate_idx[i];
    real ratio = pow(c_par[p] / x[i], b[p]);
    real mu_i = a[p] + (d[p] - a[p]) / (1.0 + ratio);
    real sigma_i;
    if (use_heteroscedastic_noise) {
      real log_abs_mu = log(abs(mu_i) + 1e-10);
      sigma_i = exp(log_sigma0 + log_sigma_slope * log_abs_mu);
    } else {
      sigma_i = sigma_obs;
    }
    y[i] ~ student_t(nu, mu_i, sigma_i);
  }

  if (N_blanks > 0)
    blank_response ~ student_t(nu, a[blank_plate_idx], sigma_blank);
}

generated quantities {
  vector[N_obs] y_pred;
  vector[N_obs] log_lik;

  for (i in 1:N_obs) {
    int p = plate_idx[i];
    real ratio = pow(c_par[p] / x[i], b[p]);
    real mu_val = a[p] + (d[p] - a[p]) / (1.0 + ratio);
    real sigma_i;
    if (use_heteroscedastic_noise) {
      real log_abs_mu = log(abs(mu_val) + 1e-10);
      sigma_i = exp(log_sigma0 + log_sigma_slope * log_abs_mu);
    } else {
      sigma_i = sigma_obs;
    }
    y_pred[i]  = student_t_rng(nu, mu_val, sigma_i);
    log_lik[i] = student_t_lpdf(y[i] | nu, mu_val, sigma_i);
  }
}
