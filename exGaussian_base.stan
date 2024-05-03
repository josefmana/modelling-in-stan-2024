
data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  int prior_only;  // should the likelihood be ignored?
}

parameters {
  real Intercept;  // temporary intercept for centered predictors
  real<lower=0> sigma;  // dispersion parameter
  real<lower=0> lambda;  // scale parameter
}

transformed parameters {
  real lprior = 0;  // prior contributions to the log posterior
  lprior += student_t_lpdf(Intercept | 3, 742, 296.5);
  lprior += student_t_lpdf(sigma | 3, 0, 296.5)
    - 1 * student_t_lccdf(0 | 3, 0, 296.5);
  lprior += gamma_lpdf(lambda | 1, 0.1);
}

model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] mu = rep_vector(0.0, N);
    mu += Intercept;
    target += exp_mod_normal_lpdf( Y | mu - beta, sigma, inv(lambda) );
  }
  // priors including constants
  target += lprior;
}

generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept;
}