
# clear environment
rm( list = ls() )

# prepare environment
library(cmdstanr)
library(posterior)
library(bayesplot)
color_scheme_set("brightblue")


# Task 1: Bayesian “t-test”

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


# Task 2: Convert to “long format”

# convert the data
dlist2 <- list( N = dlist$N1 + dlist$N2, y = c( dlist$yA, dlist$yB ), isB = c( rep(0,dlist$N1), rep(1,dlist$N2) ) )

# prepare the model
file2 <- file.path("240325_lesson_no_tres","ttest_lm1.stan")
mod2 <- cmdstan_model(file2)
fit2 <- mod2$sample( data = dlist2, chains = 4 ) # fit it

fit2$summary()

