# This script is used to run ExGaussian model of SSRT data implemented in Stan.

rm( list = ls() )
options( mc.cores = parallel::detectCores() )

# recommended running this is a fresh R session or restarting current session
# install.packages( "cmdstanr", repos = c( "https://mc-stan.org/r-packages/", getOption("repos") ) )

library(here)
library(tidyverse)
library(brms) # for rexgaussian()
library(bayesplot)
library(cmdstanr)

color_scheme_set("viridisA")
theme_set( theme_bw(base_size = 14) )


# DATA PRE-PROCESSING ----

d0 <- read.csv( here("_raw","ssrt_lab.csv"), sep = "," ) # read data

# pivot it wider
d1 <- d0 %>%
  
  filter( block != 0 ) %>%
  select( id, block, cond, correct, signal, rt1, trueSOA ) %>%
  rename( "rt" = "rt1", "ssd" = "trueSOA" ) %>%
  mutate( across( c("rt","ssd"), ~ .x/1e3 ) ) # re-scale from ms to s

# GO response times should be generally slower than STOP-RESPOND response times
d1 %>%
  
  # some formatting shinanigans
  filter( complete.cases(rt) ) %>%
  mutate(
    `Response type: ` = ifelse( signal == "nosignal", "GO", "SIGNAL-RESPOND" ),
    `Stimulation type: ` = ifelse( cond == "ctrl", "High frequency", "Combined frequency" )
  ) %>%
  
  # plot it
  ggplot() +
  aes(x = rt, colour = `Response type: `, linetype = `Stimulation type: `  ) +
  geom_density(linewidth = 1.25) +
  labs(x = "Response time (s)", y = "Density") +
  scale_colour_manual( values = c("grey","black") ) +
  facet_wrap(~id, ncol = 2) +
  theme( legend.position = "bottom", panel.grid = element_blank() )

# GO response times should be generally slower than STOP-RESPOND response times, no.2
d1 %>%
  
  # some formatting shinanigans
  filter( complete.cases(rt) ) %>%
  mutate(
    `Response type: ` = ifelse( signal == "nosignal", "GO", "SIGNAL-RESPOND" ),
    `Stimulation type: ` = ifelse( cond == "ctrl", "High frequency", "Combined frequency" )
  ) %>%
  
  # plot it
  ggplot() +
  aes(x = rt, fill = `Stimulation type: `  ) +
  geom_density(alpha = .3) +
  labs(x = "Response time (s)", y = "Density") +
  facet_grid( id ~ `Response type: ` ) +
  theme( legend.position = "bottom", panel.grid = element_blank() )


# INDIVIDUAL MODEL ----

# prepare the model
mod_indi <- cmdstan_model( here("ExGaussian_individual.stan") )

# function for manual initial values setting
ifun <- function() list(
  
  Int_mu_go_0 = runif(1,-2,0),
  Int_sigma_go_0 = runif(1,-2,0),
  Int_lambda_go_0 = runif(1,-2,0),
  Int_mu_stop_0 = runif(1,-2,0),
  Int_sigma_stop_0 = runif(1,-2,0),
  Int_lambda_stop_0 = runif(1,-2,0)
  
)

# extract number of participants
k <- unique(d1$id)

# prepare a list for the fits
fit_many_indi <- lapply( setNames(k,k), function(i) list() )

# fit it
for ( i in k ) for ( x in c("ctrl","exp") ) {
  
  # data
  dGO <- subset(d1, id == i & cond == x & complete.cases(rt)) %>% filter( signal == "nosignal" )
  dSR <- subset(d1, id == i & cond == x & complete.cases(rt) & (rt>ssd) ) %>% filter( signal == "signal" & !is.na(rt) )
  dNA <- subset(d1, id == i & cond == x) %>% filter( signal == "signal" & is.na(rt) )
  
  # input file
  dlist <- list(
    
    # data
    Y_0_go = dGO$rt, N_0_go = nrow(dGO),
    Y_0_sr = dSR$rt, N_0_sr = nrow(dSR), SSD_0_sr = dSR$ssd,
    N_0_na = nrow(dNA), SSD_0_na = dNA$ssd,
    
    # priors
    mu_go_0_p = c(-0.4,0.2), sigma_go_0_p = c(-2.0,0.2), lambda_go_0_p = c(-2.0,0.2),
    mu_stop_0_p = c(-1.0,0.2), sigma_stop_0_p = c(-2.0,0.2), lambda_stop_0_p = c(-2.0,0.2)
    
  )
  
  # fitting proper
  fit_many_indi[[i]][[x]] <- mod_indi$sample(data = dlist, chains = 4, save_warmup = T, init = ifun)

}

# With current priors, model, and data, fitting IPN390 control condition block (k = "IPN390", x = "ctrl")
# leads to relatively larger percentage of divergent transitions compared to other data
# and occasional non-converging chain.
