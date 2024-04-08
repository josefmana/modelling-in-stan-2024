//
// This Stan program defines a Bayesian regression
//

data {
   int<lower=0> N;
   vector[N] y;
   vector[N] x;
   vector[N] isB;
}

parameters {
  real mu;
  real delta;
  real beta;
  real<lower=0> sigma;
}

model {
  // likelihoods
  target += normal_lpdf( y | mu + beta * x + delta * isB, sigma );
  
  // priors
  target += normal_lpdf( mu | 0, 1 );
  target += normal_lpdf( beta | 0, .5 );
  target += normal_lpdf( delta | 0, .5 );
  target += normal_lpdf( sigma | 0, 2 );
  
}


