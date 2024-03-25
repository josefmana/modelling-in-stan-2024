//
// This Stan program defines a 'Bayesian t-test' model
//

data {
  
  // observations
  int<lower=0> N1;
  int<lower=0> N2;
  vector[N1] yA;
  vector[N2] yB;

}

parameters {
  real mu;
  real delta;
  real<lower=0> sigma_A;
  real<lower=0> sigma_B;
}

model {
  // likelihoods
  target += normal_lpdf( yA | mu, sigma_A );
  target += normal_lpdf( yB | ( mu + delta ), sigma_B );
  
  // priors
  target += normal_lpdf( mu | 0, 1 );
  target += normal_lpdf( delta | 0, 1.0/2.0 );
  target += normal_lpdf( sigma_A | 0, 1 );
  target += normal_lpdf( sigma_B | 0, 1 );
  
}


