
rm( list = ls() )
options( mc.cores = parallel::detectCores() )

library(brms)
library(tidyverse)
library(bayesplot)
library(cmdstanr)

# read data
d0 <- read.csv( "/Users/josefmana/Desktop/crusades/dbs/combSTIM/_data/ssrt_lab.csv", sep = "," )

# pivot it wider
d1 <- d0 %>% filter( block != 0 ) %>% select( id, block, cond, signal, rt1, trueSOA ) %>% rename( "rt" = "rt1", "ssd" = "trueSOA" )

# separate files for "go" and "stop" trials
Dgo <- d1 %>% filter( signal == "nosignal" )
Dsr <- d1 %>% filter( signal == "signal" & !is.na(rt) )
Dna <- d1 %>% filter( signal == "signal" & is.na(rt) )

# print brms exGaussian model
#stancode( bf( rt ~ 1, sigma ~ 1, beta ~ 1 ), data = Dgo, family = exgaussian() )
stancode( bf( rt ~ 1 ), data = Dgo, family = exgaussian() )

# https://github.com/twiecki/stopsignal.git/src/stop_likelihoods.pyx
# http://dx.doi.org/10.1037/a0030543

# found out the parametrisation is different in Stan
# (https://mc-stan.org/docs/2_22/functions-reference/exponentially-modified-normal-distribution.html) compared to the
# source article (http://dx.doi.org/10.1037/a0030543)

# also found out that default brms priors are way too restrictive:

## tried generate prior predictive distributions vai rexgaussian() and with brms defaults got distributions shifted to
## a lot of different places between 300-900 ms mean but always very close to the mean

## tried to fit some models with much wider priors (inspired by the original article) using brms only (so far) to seem
## how well/bad will the model behave

# BRMS MODELS ----

df <- Dgo[ complete.cases(Dgo$rt) , ]
df$idcond <- paste( df$id, df$cond, sep = "_" )

## common mean ----

p0 <- c( prior( normal(0,1e6), class = Intercept, lb = 1, ub = 1000 ),
         prior( normal(0,1e4), class = Intercept, dpar = beta, lb = 1, ub = 300 ),
         prior( normal(0,1e4), class = Intercept, dpar = sigma, lb = 1, ub = 300 )
         )

m0 <- brm( bf( rt ~ 1, sigma ~ 1, beta ~ 1 ),
           family = exgaussian(link_sigma = "identity", link_beta = "identity"),
           prior = p0,
           data = df
           )

ppc_stat_grouped( y = df$rt, yrep = posterior_predict( m0, newdata = df ), group = df$idcond, stat = "median" )


## stratified by participant ----

# not working very consistently - quite often on chains ends up stuck

p1 <- c( prior( normal(780,100), class = Intercept ),
         prior( normal(0,1e4), class = b, lb = -300, ub = 300 ),
         prior( normal(0,1e4), class = Intercept, dpar = sigma, lb = 1, ub = 300 ),
         prior( normal(0,1e4), class = Intercept, dpar = beta, lb = 1, ub = 300 )
         )

m1 <- brm( bf( rt ~ 1 + id, sigma ~ 1, beta ~ 1 ),
           family = exgaussian(link_sigma = "identity", link_beta = "identity"),
           prior = p1,
           data = df
           )

ppc_stat_grouped( y = df$rt, yrep = posterior_predict( m1, newdata = df ), group = df$idcond, stat = "median" )


## participant/condition interaction ----

p2 <- c( prior( normal(780,100), class = Intercept ),
         prior( normal(0,1e4), class = b, lb = -300, ub = 300 ),
         prior( normal(0,1e4), class = Intercept, dpar = sigma, lb = 1, ub = 300 ),
         prior( normal(0,1e4), class = Intercept, dpar = beta, lb = 1, ub = 300 )
         )

m2 <- brm( bf( rt ~ 1 + id * cond, sigma ~ 1, beta ~ 1 ),
           family = exgaussian(link_sigma = "identity", link_beta = "identity"),
           prior = p2,
           data = df
           )

ppc_stat_grouped( y = df$rt, yrep = posterior_predict( m2, newdata = df ), group = df$idcond, stat = "median" )


## fixed condition/varying participant ----

p3 <- c( prior( normal(780,100), class = Intercept ),
         prior( normal(0,1e4), class = b, lb = -300, ub = 300 ),
         prior( normal(0,1e4), class = Intercept, dpar = sigma, lb = 1, ub = 300 ),
         prior( normal(0,1e4), class = Intercept, dpar = beta, lb = 1, ub = 300 ),
         prior( normal(0,1e4), class = sd, lb = 0, ub = 300 )
         )

m3 <- brm( bf( rt ~ 1 + cond + (1 | id ), sigma ~ 1, beta ~ 1 ),
           family = exgaussian(link_sigma = "identity", link_beta = "identity"),
           prior = p3,
           data = df
           )

ppc_stat_grouped( y = df$rt, yrep = posterior_predict( m3, newdata = df ), group = df$idcond, stat = "median" )

