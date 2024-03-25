//
// This Stan program defines a 'Bayesian t-test' model
//

data {
   int<lower=0> N;
   vector[N] y;
   vector[N] isB;
}

parameters {
  real mu;
  real delta;
  real<lower=0> sigma;
  //real<lower=0> sigma_B;
}

model {
  // likelihoods
  target += normal_lpdf( y | mu + delta * isB, sigma );
  
  // priors
  target += normal_lpdf( mu | 0, 1 );
  target += normal_lpdf( delta | 0, 1.0/2.0 ); // cannot be 1/2!
  target += normal_lpdf( sigma | 0, 1 );
  
}


