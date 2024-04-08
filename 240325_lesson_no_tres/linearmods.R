
# clear environment
rm( list = ls() )

# prepare environment
library(cmdstanr)
library(posterior)
library(bayesplot)
color_scheme_set("brightblue")


###  Task 1: Bayesian “t-test”

# function to generate data
#gendat0 <- function( mean_mu = 0, sd_mu = 1, mean_sigma = 0, sd_sigma = 2, grps = 2, n = c(30,30) ) {
#  
#  mu <- rnorm( grps, mean_mu, sd_mu )
#  sigma <- abs( rnorm( grps, mean_sigma, sd_sigma ) )
#  
#  output <- list( mu = mu, sigma = sigma, delta = mu[2]-mu[1] )
#  for ( i in 1:grps ) output[[ paste0("y",i)]] <- rnorm( n[i], mu, sigma )
#  
#  return(output)
#  
#}

gendat1 <-
  
  function( mean_mu = 0, sd_mu = 1, mean_delta = 0, sd_delta = .5, mean_sigma = 0, sd_sigma = 2, n = c(33,66) ) {
    
    mu <- rnorm( 1, mean_mu, sd_mu )
    delta <- rnorm( 1, mean_delta, sd_delta )
    sigmaA <- abs( rnorm( 1, mean_sigma, sd_sigma ) )
    sigmaB <- abs( rnorm( 1, mean_sigma, sd_sigma ) )
    
    output <- list( mu = mu,
                    delta = delta,
                    sigmaA = sigmaA,
                    sigmaB = sigmaB,
                    yA = rnorm( n[1], mu, sigmaA ),
                    yB = rnorm( n[2], mu + delta, sigmaB )
                    )
    
    
    return(output)
    
  }

# prepare data set
d0 <- gendat1( n = c(100,100) )
dlist <- with( d0, list( N1 = length(yA), N2 = length(yB), yA = yA, yB = yB ) )

# prepare the model
file <- file.path("240325_lesson_no_tres","ttest.stan")
mod <- cmdstan_model(file)
fit <- mod$sample( data = dlist, chains = 4 ) # fit it

fit$summary()


###  Task 2: Convert to “long format”

# convert the data
dlist2 <- list( N = dlist$N1 + dlist$N2, y = c( dlist$yA, dlist$yB ), isB = c( rep(0,dlist$N1), rep(1,dlist$N2) ) )

# prepare the model
file2 <- file.path("240325_lesson_no_tres","ttest_lm1.stan")
mod2 <- cmdstan_model(file2)
fit2 <- mod2$sample( data = dlist2, chains = 4 ) # fit it

fit2$summary()

###  Task 3: Add a continuous predictor
###  Add a new continuous predictor (a vector of real numbers) for each data point and add an extra coefficient (a new variable in the parameters block) to model its influence.
###  Simulate data and check parameter recovery.

gendat3 <-
  
  function( mean_mu = 0, sd_mu = 1,
            mean_delta = 0, sd_delta = .5,
            mean_sigma = 0, sd_sigma = 2,
            mean_slope = 0, sd_slope = .5,
            n = c(100,100) ) {
    
    mu <- rnorm( 1, mean_mu, sd_mu )
    delta <- rnorm( 1, mean_delta, sd_delta )
    slope <- rnorm( 1, mean_slope, sd_slope )
    sigma <- abs( rnorm( 1, mean_sigma, sd_sigma ) )
    x <- rnorm( n[1]+n[2], 0, 1 )
    
    output <- list( mu = mu,
                    beta = slope,
                    delta = delta,
                    sigmaA = sigma,
                    sigmaB = sigma,
                    x = x,
                    yA = rnorm( n[1], mu + slope * x, sigma ),
                    yB = rnorm( n[2], mu + slope * x + delta, sigma )
                    )
    
    return(output)
    
  }

# prepare data list
d3 <- gendat3( n = c(89,95) )
dlist3 <- with( d3, list( N = length(x), y = c(yA,yB), x = x, isB = c( rep(0,length(yA)), rep(1,length(yB)) ) ) )

# prepare the model
file3 <- file.path("240325_lesson_no_tres","lm2.stan")
mod3 <- cmdstan_model(file3)
fit3 <- mod3$sample( data = dlist3, chains = 4 ) # fit it

# check it
fit3$summary()
unlist( d3[ c("mu","beta","delta","sigmaA","sigmaB") ] )


###  Task 4: Matrix multiplication
###  Make a copy of the model and convert it to matrix multiplication format.

# prepare data
dlist4 <- with( dlist3, list( y = y, N = N, K = 2, x = cbind(isB,x), beta_mean = rep(0,2), beta_sd = rep(.5,2) ) )

# prepare the model
file4 <- file.path("240325_lesson_no_tres","lm2_but_cooler.stan")
mod4 <- cmdstan_model(file4)
fit4 <- mod4$sample( data = dlist4, chains = 4 ) # fit it

# check it
fit3$summary()
fit4$summary()
unlist( d3[ c("mu","delta","beta","sigmaA","sigmaB") ] )


### Task 5: Compare to user’s guide
### Compare what you’ve built with the example linear regression models in Stan’s User’s guide: https://mc-stan.org/docs/stan-users-guide/regression.html#linear-regression
### Generally, don’t be afraid to use the User’s guide as a starting point for your models, it is a great resource!

TRUE

### Task 6: Dummy coding
### Use dummy coding to extend the model to have 3 groups instead of 2 in the categorical predictor. Your Stan code should not need to change.
### Simulate data, fit, check parameter recovery.

gendat6 <-
  
  function( mean_mu = 0, sd_mu = 1,
            mean_delta = 0, sd_delta = .5,
            mean_slope = 0, sd_slope = .5,
            mean_sigma = 0, sd_sigma = 2,
            n = c(100,100,101) ) {
    
    mu <- rnorm( 1, mean_mu, sd_mu )
    slope <- rnorm( 1, mean_slope, sd_slope )
    sigma <- abs( rnorm( 1, mean_sigma, sd_sigma ) )
    k <- length(n)
    
    x <- lapply( 1:k, function(i) rnorm( n[i], 0, 1 ) )
    delta <- sapply( 1:k, function(i) ifelse( i == 1, 0, rnorm( 1, mean_delta, sd_delta ) ) )
    y <- lapply( 1:k, function(i) rnorm( n[i], mu + x[[i]] * slope + delta[i], sigma ) )
    
    output <- list( mu = mu,
                    beta = slope,
                    delta = delta,
                    k = k,
                    n = n,
                    sigma = sigma,
                    x = unlist(x),
                    y = unlist(y)
    )
    
    return(output)
    
  }

# generate data
d6 <- gendat6()

# prepare dummy coded prediction matrix
d6$xmat <- matrix( nrow = sum(d6$n), ncol = d6$k-1 )

for (i in 2:d6$k) {
  if ( i == d6$k ) d6$xmat[ ,i-1] <-  with( d6, c( rep( 0, sum( n[1:(i-1)] ) ), rep( 1, n[i] ) ) )
  else d6$xmat[ ,i-1] <-  with( d6, c( rep( 0, sum( n[1:(i-1)] ) ), rep( 1, n[i] ), rep( 0, sum( n[(i+1):k] ) ) ) )
}

# prepare the list
dlist6 <- with( d6, list( N = length(y), K = k, y = y, x = cbind( xmat, x ), beta_mean = rep(0,k), beta_sd = rep(.5,k) ) )

# prepare the model
file6 <- file.path("240325_lesson_no_tres","lm2_but_cooler.stan")
mod6 <- cmdstan_model(file6)
fit6 <- mod6$sample( data = dlist6, chains = 4, cores = 8 ) # fit it

# check it
fit6$summary()
unlist( d6[ c("mu","delta","beta","sigma") ] )


### Task 7: Negative binomial regression
### Starting with the previous model, create a negative binomial regression model with a log link.
### Simulate data, fit, check parameter recovery (you can use a simpler set of predictors if you prefer).


