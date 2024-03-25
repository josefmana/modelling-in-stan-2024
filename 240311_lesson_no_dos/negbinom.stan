//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as Negatively Binomial distributed
// with lambda rate parameter which has prior for log(lambda)
// parametrised by normal distribution with mean 'mu' and standard deviation 'sigma'.
//

data {
  int<lower=0> N;
  array[N] int y;
  real mu;
  real<lower=0> sigma;
  real<lower=0> beta;
}

parameters {
  real log_lambda;
  real<lower=0> phi;
}

model {
  // likelihood
  target += neg_binomial_2_lpmf( y | exp(log_lambda), phi );
  
  // prior
  target += normal_lpdf( log_lambda | mu, sigma );
  target += exponential_lpdf( phi | beta );
  
}
