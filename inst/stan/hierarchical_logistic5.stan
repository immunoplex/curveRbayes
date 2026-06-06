// hierarchical_logistic5.stan
//
// Hierarchical 5-parameter logistic model for immunoassay standard curves.
// Uses the curveRcore b>0 always-increasing convention:
//
//   y = a + (d - a) / (1 + exp(-(x - c) / b))^g
//
// where x = log10(concentration), y = log10(response) [or raw response].
//
// Hierarchical structure: plate-level parameters are non-centered
// log-normal deviations from population means.
//
// All prior hyperparameters are passed as DATA from R
// (computed by compute_dynamic_priors()), making the model scale-invariant.

data {
  int<lower=1> N_obs;
  int<lower=1> N_plates;
  int<lower=1> grainsize;
  array[N_obs] int<lower=1, upper=N_plates> plate_idx;
  vector[N_obs] x;                      // log10(concentration)
  array[N_obs] real y;                   // response (log10 or raw)

  // Prior hyperparameters (all computed in R)
  real prior_a_mu;
  real<lower=0> prior_a_sigma;
  real prior_d_mu;
  real<lower=0> prior_d_sigma;
  real prior_log_b_mu;                   // b on log scale for positivity
  real<lower=0> prior_log_b_sigma;
  real prior_c_mu;                       // c on fitting scale (log10 conc)
  real<lower=0> prior_c_sigma;
  real<lower=0> prior_log_g_sd;          // g regularisation toward g=1 (log(g)~N(0,sd))
  real<lower=0> prior_log_g_plate_sd;

  // Blank data (optional)
  int<lower=0> N_blanks;
  array[N_blanks] int<lower=1, upper=N_plates> blank_plate_idx;
  vector[N_blanks] blank_response;
}

parameters {
  // Population-level (hyperpriors)
  real mu_a;
  real<lower=0> sigma_a;
  real mu_d;
  real<lower=0> sigma_d;
  real mu_log_b;
  real<lower=0> sigma_log_b;
  real mu_c;
  real<lower=0> sigma_c;
  real mu_log_g;
  real<lower=0> sigma_log_g;

  // Plate-level (non-centered)
  vector[N_plates] raw_a;
  vector[N_plates] raw_d;
  vector[N_plates] raw_log_b;
  vector[N_plates] raw_c;
  vector[N_plates] raw_log_g;

  // Noise model
  real<lower=0> sigma_obs;               // observation noise SD
  real<lower=2> nu;                      // Student-t df
  real<lower=0> sigma_blank;             // blank noise SD
}

transformed parameters {
  // Plate-level parameters (non-centered → actual)
  vector[N_plates] a = mu_a + sigma_a * raw_a;
  vector[N_plates] d = mu_d + sigma_d * raw_d;
  vector[N_plates] log_b = mu_log_b + sigma_log_b * raw_log_b;
  vector[N_plates] b = exp(log_b);      // b > 0 always
  vector[N_plates] c_par = mu_c + sigma_c * raw_c;
  vector[N_plates] log_g = mu_log_g + sigma_log_g * raw_log_g;
  vector[N_plates] g = exp(log_g);      // g > 0 always
}

model {
  // Hyperpriors
  mu_a ~ normal(prior_a_mu, prior_a_sigma);
  sigma_a ~ normal(0, prior_a_sigma * 0.5);

  mu_d ~ normal(prior_d_mu, prior_d_sigma);
  sigma_d ~ normal(0, prior_d_sigma * 0.5);

  mu_log_b ~ normal(prior_log_b_mu, prior_log_b_sigma);
  sigma_log_b ~ normal(0, 0.5);

  mu_c ~ normal(prior_c_mu, prior_c_sigma);
  sigma_c ~ normal(0, prior_c_sigma * 0.5);

  mu_log_g ~ normal(0, prior_log_g_sd);
  sigma_log_g ~ normal(0, prior_log_g_plate_sd);

  // Non-centered random effects
  raw_a ~ std_normal();
  raw_d ~ std_normal();
  raw_log_b ~ std_normal();
  raw_c ~ std_normal();
  raw_log_g ~ std_normal();

  // Noise priors
  sigma_obs ~ normal(0, prior_a_sigma);
  nu ~ gamma(2, 0.1);
  sigma_blank ~ normal(0, prior_a_sigma);

  // Likelihood: curveRcore logistic5 convention
  //   y = a + (d - a) / (1 + exp(-(x - c) / b))^g
  for (i in 1:N_obs) {
    int p = plate_idx[i];
    real z = -(x[i] - c_par[p]) / b[p];
    real mu_i = a[p] + (d[p] - a[p]) / pow(1.0 + exp(z), g[p]);
    y[i] ~ student_t(nu, mu_i, sigma_obs);
  }

  // Blank likelihood (anchors lower asymptote)
  if (N_blanks > 0)
    blank_response ~ student_t(nu, a[blank_plate_idx], sigma_blank);
}

generated quantities {
  vector[N_obs] y_pred;
  vector[N_obs] log_lik;

  for (i in 1:N_obs) {
    int p = plate_idx[i];
    real z = -(x[i] - c_par[p]) / b[p];
    real mu_val = a[p] + (d[p] - a[p]) / pow(1.0 + exp(z), g[p]);
    y_pred[i] = student_t_rng(nu, mu_val, sigma_obs);
    log_lik[i] = student_t_lpdf(y[i] | nu, mu_val, sigma_obs);
  }
}
