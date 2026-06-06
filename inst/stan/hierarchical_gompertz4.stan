// hierarchical_gompertz4.stan
//
// Hierarchical 4-parameter Gompertz model for immunoassay standard curves.
// Uses the curveRcore convention:
//
//   y = a + (d - a) * exp(-exp(-b * (x - c)))
//
// where x = log10(concentration), y = log10(response) [or raw response].
// Intrinsically asymmetric — no need for a 5th parameter.

data {
  int<lower=1> N_obs;
  int<lower=1> N_plates;
  int<lower=1> grainsize;
  array[N_obs] int<lower=1, upper=N_plates> plate_idx;
  vector[N_obs] x;
  array[N_obs] real y;

  // Prior hyperparameters
  real prior_a_mu;
  real<lower=0> prior_a_sigma;
  real prior_d_mu;
  real<lower=0> prior_d_sigma;
  real prior_log_b_mu;
  real<lower=0> prior_log_b_sigma;
  real prior_c_mu;
  real<lower=0> prior_c_sigma;

  // Blank data
  int<lower=0> N_blanks;
  array[N_blanks] int<lower=1, upper=N_plates> blank_plate_idx;
  vector[N_blanks] blank_response;
}

parameters {
  real mu_a;
  real<lower=0> sigma_a;
  real mu_d;
  real<lower=0> sigma_d;
  real mu_log_b;
  real<lower=0> sigma_log_b;
  real mu_c;
  real<lower=0> sigma_c;

  vector[N_plates] raw_a;
  vector[N_plates] raw_d;
  vector[N_plates] raw_log_b;
  vector[N_plates] raw_c;

  real<lower=0> sigma_obs;
  real<lower=2> nu;
  real<lower=0> sigma_blank;
}

transformed parameters {
  vector[N_plates] a = mu_a + sigma_a * raw_a;
  vector[N_plates] d = mu_d + sigma_d * raw_d;
  vector[N_plates] log_b = mu_log_b + sigma_log_b * raw_log_b;
  vector[N_plates] b = exp(log_b);
  vector[N_plates] c_par = mu_c + sigma_c * raw_c;
}

model {
  mu_a ~ normal(prior_a_mu, prior_a_sigma);
  sigma_a ~ normal(0, prior_a_sigma * 0.5);
  mu_d ~ normal(prior_d_mu, prior_d_sigma);
  sigma_d ~ normal(0, prior_d_sigma * 0.5);
  mu_log_b ~ normal(prior_log_b_mu, prior_log_b_sigma);
  sigma_log_b ~ normal(0, 0.5);
  mu_c ~ normal(prior_c_mu, prior_c_sigma);
  sigma_c ~ normal(0, prior_c_sigma * 0.5);

  raw_a ~ std_normal();
  raw_d ~ std_normal();
  raw_log_b ~ std_normal();
  raw_c ~ std_normal();

  sigma_obs ~ normal(0, prior_a_sigma);
  nu ~ gamma(2, 0.1);
  sigma_blank ~ normal(0, prior_a_sigma);

  // Likelihood: curveRcore gompertz4 convention
  for (i in 1:N_obs) {
    int p = plate_idx[i];
    real u = exp(-b[p] * (x[i] - c_par[p]));
    real mu_i = a[p] + (d[p] - a[p]) * exp(-u);
    y[i] ~ student_t(nu, mu_i, sigma_obs);
  }

  if (N_blanks > 0)
    blank_response ~ student_t(nu, a[blank_plate_idx], sigma_blank);
}

generated quantities {
  vector[N_obs] y_pred;
  vector[N_obs] log_lik;

  for (i in 1:N_obs) {
    int p = plate_idx[i];
    real u = exp(-b[p] * (x[i] - c_par[p]));
    real mu_val = a[p] + (d[p] - a[p]) * exp(-u);
    y_pred[i] = student_t_rng(nu, mu_val, sigma_obs);
    log_lik[i] = student_t_lpdf(y[i] | nu, mu_val, sigma_obs);
  }
}
