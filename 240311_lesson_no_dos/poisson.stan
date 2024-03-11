//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as Poisson distributed
// with lambda rate parameter which has prior for log(lambda)
// parametrised by normal distribution with mean 'mu' and standard deviation 'sigma'.
//

data {
  int<lower=0> N;
  array[N] int y;
  real mu;
  real<lower=0> sigma;
}

parameters {
  real log_lambda;
  //real<lower=0> lambda;
}

model {
  // likelihood
  target += poisson_lpmf( y | exp(log_lambda) );
  // target += poisson_lpmf( y | lambda );
  
  // prior
  target += normal_lpdf( log_lambda | mu, sigma );
  // target += lognormal_lpdf( lambda | mu, sigma );
  
}
