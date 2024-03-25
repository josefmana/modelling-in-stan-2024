# This scripts taks a Poisson distribution for a run through statistical workflow with prior and posterior predictive checks.

# clear environment
rm( list = ls() )

# prepare environment
library(cmdstanr)
library(posterior)
library(bayesplot)
color_scheme_set("brightblue")

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
  
  function( m = mod, obs = 30, showtrue = T ) {
    
    # generate data
    data <- gennbinom_norm( N = obs, phi = .5, mu = 1.2, sigma = .8, print = T )
    dlist <- list( N = length(data$y), y = data$y, mu = 1.2, sigma = .8 )
    
    fit <- m$sample( data = dlist, chains = 4 ) # fit it
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

# What is the smallest number of observations we need to reliably detect the problem?
# run the check few times
for ( i in 1:10 ) overdispcheck( obs = 30, showtrue = F )


# Implement a negative binomial model and see how the check behaves now.
# Use neg_binomial_2_lpmf)
# Put an exponential(1) prior on the overdispersion parameter.

# read the model
file3 <- file.path("240311_lesson_no_dos","negbinom.stan")
mod3 <- cmdstan_model(file3)

# generate data
data3 <- gennbinom_norm( N = 30, phi = .5, mu = 1.2, sigma = .8, print = T )
dlist3 <- list( N = length(data3$y), y = data3$y, mu = 1.2, sigma = .8, beta = 1 )

# fit it
fit3 <- mod3$sample( data = dlist3, chains = 4 )

# summarise posteriors
fit3$summary()

# extract draws
drws3 <- fit3$draws( format = "matrix" )

# prepare posterior predictions
yrep <- t( sapply( 1:nrow(drws3), function(i) rnbinom( dlist3$N, size = drws3[i,3], mu = exp( drws3[i,2] ) ) ) )

ppc_stat( dlist3$y, yrep )
ppc_stat( dlist3$y, yrep, stat = "sd" )


# Task 4: A bit more open-ended exploration ----

# the data
d4 <-
  
  cbind.data.frame(
    y = c( 2, 0, 11, 25, 9, 4, 17, 11, 8, 4, 5, 2, 6, 4, 8, 24, 0, 3, 6, 4, 12, 9, 5, 6, 2, 10, 4, 15, 0, 1, 87, 19, 2, 1, 38, 16, 5, 7, 18, 11, 1, 0, 7, 15, 5, 0, 6, 1, 0, 6, 34, 29, 7, 11, 5, 10, 8, 3, 12, 18, 7, 1, 18, 18, 6, 3, 37, 5, 1, 0, 22, 13, 1, 0, 26, 19, 6, 7, 21, 45 ),
    group  =  c( "A", "B", "A", "B", "A", "B", "A", "B", "A", "B", "A", "B", "A", "B", "A", "B", "A", "B", "A", "B", "A", "B", "A", "B", "A", "B", "A", "B", "A", "B", "A", "B", "A", "B", "A", "B", "A", "B", "A", "B", "A", "B", "A", "B", "A", "B", "A", "B", "A", "B", "A", "B", "A", "B", "A", "B", "A", "B", "A", "B", "A", "B", "A", "B", "A", "B", "A", "B", "A", "B", "A", "B", "A", "B", "A", "B", "A", "B", "A", "B" ),
    type =  c( "X", "X", "Y", "Y", "X", "X", "Y", "Y", "X", "X", "Y", "Y", "X", "X", "Y", "Y", "X", "X", "Y", "Y", "X", "X", "Y", "Y", "X", "X", "Y", "Y", "X", "X", "Y", "Y", "X", "X", "Y", "Y", "X", "X", "Y", "Y", "X", "X", "Y", "Y", "X", "X", "Y", "Y", "X", "X", "Y", "Y", "X", "X", "Y", "Y", "X", "X", "Y", "Y", "X", "X", "Y", "Y", "X", "X", "Y", "Y", "X", "X", "Y", "Y", "X", "X", "Y", "Y", "X", "X", "Y", "Y" )
  )

# prepare a list
dlist4 <- list( N = nrow(d4), y = d4[ ,"y"], mu = 1.2, sigma = 0.8, beta = 1 )

# fit it
fit4 <- mod3$sample( data = dlist4, chains = 4 )

# extract draws
drws4 <- fit4$draws( format = "matrix" )

# prepare posterior predictions
yrep <- t( sapply( 1:nrow(drws4), function(i) rnbinom( dlist4$N, size = drws4[i,3], mu = exp( drws4[i,2] ) ) ) )

# check means and variance
ppc_stat( dlist4$y, yrep )
ppc_stat( dlist4$y, yrep, stat = "sd" )

# group it
ppc_stat_grouped( dlist4$y, yrep, group = d4$group )
ppc_stat_grouped( dlist4$y, yrep, group = d4$group, stat = "sd" )

ppc_stat_grouped( dlist4$y, yrep, group = d4$type )
ppc_stat_grouped( dlist4$y, yrep, group = d4$type, stat = "sd" )


# Task 5: Open-ended exploration 2

# data
d5 <-
  
  cbind.data.frame(
    y = c( 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1 ),
    group = c( "B", "C", "D", "B", "C", "B", "D", "C", "A", "D", "C", "C", "C", "D", "C", "C", "B", "B", "D", "C", "D", "C", "C", "A", "B", "D", "C", "D", "D", "D", "C", "D", "A", "A", "C", "D", "D", "C", "C", "B", "B", "C", "D", "D", "A", "D", "A", "B", "C", "A" )
  )

# prepare a list
dlist5 <- list( N = nrow(d5), y = d5[ ,"y"], mu = 1.2, sigma = 0.8, beta = 1 )

# fit it
fit5 <- mod3$sample( data = dlist5, chains = 4 )

# extract draws
drws5 <- fit5$draws( format = "matrix" )

# prepare posterior predictions
yrep <- t( sapply( 1:nrow(drws4), function(i) rnbinom( dlist5$N, size = drws5[i,3], mu = exp( drws5[i,2] ) ) ) )

# check means and variance
ppc_stat( dlist5$y, yrep )
ppc_stat( dlist5$y, yrep, stat = "sd" )

# group it
ppc_stat_grouped( dlist5$y, yrep, group = d5$group )
ppc_stat_grouped( dlist5$y, yrep, group = d5$group, stat = "sd" )

# finished
