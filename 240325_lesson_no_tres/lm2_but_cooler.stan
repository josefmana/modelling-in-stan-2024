//
// This Stan program defines a Bayesian regression via Matrix multiplication
//

data {
  int<lower=0> N; // number of data items
  int<lower=0> K; // number of predictors
  matrix[N, K] x; // predictor matrix
  vector[N] y;    // outcome vector
  vector[K] beta_mean; // means for slopes' priors
  vector[K] beta_sd; // standard deviations for slopes' priors
}

parameters {
  real mu; // intercept
  vector[K] beta; // vector of predictors
  real<lower=0> sigma; // standard deviation
}

model {
  // likelihoods
  target += normal_lpdf( y | mu + x * beta, sigma );
  
  // priors
  target += normal_lpdf( mu | 0, 1 );
  target += normal_lpdf( beta | beta_mean, beta_sd );
  target += normal_lpdf( sigma | 0, 2 );
  
}


