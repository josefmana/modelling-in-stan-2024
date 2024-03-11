# This scripts taks a Poisson distribution for a run through statistical workflow with prior and posterior predictive checks.

# clear environment
rm( list = ls() )


# Task 1: prior predictive checks ----

# prepare a function for generating prior data sets from Poisson distribution with normal priors for log(lambda)
genpoiss_norm <-
  
  function( N = 100, mu = 0, sigma = 10, print = F ) {
    
    log_lambda <- rnorm( 1, mu, sigma )
    lambda <- exp(log_lambda)
    
    if(print) return( list( lambda = lambda, y = rpois( N, lambda ) ) )
    else return( rpois( N, lambda ) )
    
  }

# Let us take a vague N(0,10) prior on the logarithm of the mean of a Poisson distribution.
## Run a prior predictive check, is there anything suspicious?

# plotting histograms
plothist <-
  
  function( grid = 5, n = 89, m = 0, s = 10 ) {
    
    par( mfrow = c(grid,grid) )
    for( i in 1:grid^2) hist( genpoiss_norm( N = n, mu = m, sigma = s ), xlab = NULL, ylab = NULL, main = NULL )
    par( mfrow = c(1,1) )
    
  }

# plot some data sets' histograms for fun
plothist( grid = 5, n = 89, m = 0, s = 10 )

# summaries means and SDs from 1000 data sets
summary( replicate( 1e3, mean( genpoiss_norm( 89, 0, 10 ) ) ) )
summary( replicate( 1e3, sd( genpoiss_norm( 89, 0, 10 ) ) ) )

# not nice!

# Find a prior on the mean of Poisson distribution (or its logarithm - you pick) that matches this rough prior knowledge:
## Zeroes are quite unlikely

zeroperc <- function(y) sum( y == 0 ) / length(y) # function that calculates percentage of zeros
hist( replicate( 1e3, zeroperc( genpoiss_norm( 89, 0, 10 ) ) ) ) # get idea how zeros proportion for the initial vague prior

## Most values should be between 2 and 15
## Values above 18 are quite unlikely

betweenperc <- function( y, min=2, max=15 ) sum( y > min & y < max ) / length(y) # percentage of cases 2 < y < 15

# get idea for the initial vague prior
hist( replicate( 1e3, betweenperc( y = genpoiss_norm( 89, 0, 10 ), min = 2, max = 15 ) ) )
hist( replicate( 1e3, betweenperc( y = genpoiss_norm( 89, 0, 10 ), min = 18, max = Inf ) ) )

# prepare a function summarising all the points above
sanitycheck <-
  
  function( rep = 1e3, n = 89, m = 0, s = 10 ) {
    
    yrep <- sapply( 1:rep, function(i) genpoiss_norm( n, m , s ) )
    
    out <- apply( yrep, 2, function(x) c( zeros = zeroperc(x), `2-15` = betweenperc(x,2,15), `> 18` = betweenperc(x,18,Inf) ) )
    return( apply(out,1,mean) )
    
  }

# try a different prior
prior_mu = 1.2
prior_sigma = .8

plothist( grid = 5, n =100, m = prior_mu, s = prior_sigma )
sanitycheck( rep = 1e3, n = 100, m = prior_mu, s = prior_sigma )


# Task 2: Fitting a Poisson ----

# Build a model that takes an array of integers
# Accept arguments for prior for the logarithm of the mean as data
# Simulate 20 values from a Poisson distribution with known mean. Do you recover the values?

# generate data
data <- genpoiss_norm( N = 20, mu = 1.2, sigma = .8, print = T )

# read the model
file <- file.path("240311_lesson_no_dos","poisson.stan")
mod <- cmdstan_model(file)

# prepare a data set with prior
dlist <- list( N = length(data$y), y = data$y, mu = 1.2, sigma = .8 )

# fit it
fit <- mod$sample( data = dlist, chains = 4 )

# summarise posteriors
fit$summary()

# extract posterior draws
draws <- fit$draws( format = "draws_matrix")

# check the mean
summary( exp( draws[,"log_lambda"] ) )
#hist( exp( draws[,"log_lambda"] ), col="lightblue", main = "Posterior distribution of lambda", xlab = NULL )
#abline( v = data$lambda, col="red", lwd = 3, lty = 2 )


# Task 3: Detecting overdispersion with a posterior predictive check ----

# Simulate 30 values with neg. binomial with ϕ=1/2 (size in R’s rnbinom), fit with the model from Task 2.
# Do you recover the simulated mean of the neg. binomial?

# prepare a function for generating prior data sets from Poisson distribution with normal priors for log(lambda)
gennbinom_norm <-
  
  function( N = 30, phi = .5, mu = 1.2, sigma = .8, print = F ) {
    
    log_lambda <- rnorm( 1, mu, sigma )
    lambda <- exp(log_lambda)
    
    if(print) return( list( lambda = lambda, phi = phi, y = rnbinom( N, size = phi, mu = lambda ) ) )
    else return( rnbinom( N, size = phi, mu = lambda ) )
    
  }

# function to generate data, fit model, and check posterior means and SDs
overdispcheck <-
  
  function( obs = 30, showtrue = T ) {
    
    # generate data
    data <- gennbinom_norm( N = obs, phi = .5, mu = 1.2, sigma = .8, print = T )
    dlist <- list( N = length(data$y), y = data$y, mu = 1.2, sigma = .8 )
    
    fit <- mod$sample( data = dlist, chains = 4 ) # fit it
    draws <- fit$draws( format = "draws_matrix") # extract posterior draws
    
    # Try to detect the model-data mismatch with a posterior predictive check.
    
    # prepare a function for calculating  posterior predictions of Poisson model
    poiss_ppc <-
      
      function( post = draws[ ,"log_lambda"], df = dlist$y, fun = "mean", true = data$lambda, show = showtrue ) {
        
        n <-  length(df) # number of observations
        stat_hat <- do.call( fun, list(df) ) # sample estimate of mean
        
        # extract posterior predictions
        ppred <- sapply( post, function(i) rpois( n = n, lambda = exp(i) ) )
        ppred <- apply( ppred, 2, fun ) # extract the function from posteriors
        
        # prepare limits for the histogram
        if (show) lim <- c( min( round( c( stat_hat, true, min(ppred) ), 1 ) ), max( round( c( stat_hat, true, max(ppred) ), 1 ) ) )
        else lim <- c( min( round( c( stat_hat, min(ppred) ), 1 ) ), max( round( c( stat_hat, max(ppred) ), 1 ) ) )
        
        # plot it
        hist( ppred, col="lightblue", main = paste0("Posterior predictive check for ",fun), xlab = NULL, xlim = lim )
        abline( v = stat_hat, col = "blue", lwd = 8, lty = 1 )
        if (show) abline( v = true, col="red", lwd = 4, lty = 2 )
        if (show) legend( "topleft" , c("Sample estimate","Data-generating value"), fill = c("blue","red"), border = rep(NA,4) )
        
      }
    
    # plot the checks
    par( mfrow = c(2,1) )
    poiss_ppc( post = draws[ ,"log_lambda"], df = dlist$y, fun = "mean", true = data$lambda, show = showtrue )
    poiss_ppc( post = draws[ ,"log_lambda"], df = dlist$y, fun = "sd", true = sqrt( with( data, lambda + ( lambda^2 / phi ) ) ), show = showtrue )
    par( mfrow = c(1,1) )

  }

# run it few times
for ( i in 1:10 ) overdispcheck( obs = 30, showtrue = F )

# What is the smallest number of observations we need to reliably detect the problem?
  

# Implement a negative binomial model and see how the check behaves now.