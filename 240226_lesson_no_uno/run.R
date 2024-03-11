# This is a script running Stan models built in the first lesson.
# The steps are based on https://mc-stan.org/cmdstanr/articles/cmdstanr.html#compiling-a-model

# clear environment
rm( list = ls() )

# TASK 1 ----

# prepare environment
library(cmdstanr)
library(posterior)
library(bayesplot)
color_scheme_set("brightblue")


# TASK 2 ----

# read the model
file <- file.path("240226_lesson_no_uno","initial.stan")
mod <- cmdstan_model(file)

# set-up data (N = 10 data points from Normal(mu,1) )
N <- 10
mu <- rnorm(1,0,1) # draw mu from its prior
dlist <- list( N = N, y = rnorm(N,mu,1) ) # generate data

# fit it
fit <-
  
  mod$sample(
    data = dlist,
    #seed = 87542,
    chains = 4
  )

# summarise posteriors
fit$summary() # default
fit$summary( "mu", pr_lt_zero = ~ mean(. < 0) ) # use a formula to summarize arbitrary functions, e.g. Pr(mu < 0)

# summarise with some more stuff
fit$summary(
  variables = NULL,
  posterior::default_summary_measures(),
  extra_quantiles = ~posterior::quantile2( ., probs = c(.0275, .975) ),
  pr_lt_zero = ~ mean(. < 0)
)

# extract draws
draws_ar <- fit$draws() # as array
draws_df <- fit$draws(format = "df") # as data.frame
mcmc_hist( fit$draws("mu") ) # plot draws


# TASK 3 ----

# read the model
file3 <- file.path("240226_lesson_no_uno","sigma.stan")
mod3 <- cmdstan_model(file3)

# generate data from the model (N = 10 data points from Normal(mu,sigma) )
N2 <- 55
mu <- rnorm( 1, 0, 1 )
sigma <- abs( rnorm( 1, 0,2 ) ) # absolute value of rnorm to get half normal
dlist2 <- list( N = N2, y = rnorm( N2, mu, sigma ) )

# fit it
fit3 <-
  
  mod3$sample(
    data = dlist2,
    #seed = 87542,
    chains = 4
  )

# summarise posteriors
fit3$summary()


# TASK 4 ----

dlist2 <- list( N = 1, y = rnorm( 1, mu, sigma ) )

# fit it
fit3 <-
  
  mod3$sample(
    data = dlist2,
    #seed = 87542,
    chains = 4
  )

# summarise posteriors
fit3$summary()

