# This script is used to run ExGaussian model of SSRT data implemented in Stan.

rm( list = ls() )
options( mc.cores = parallel::detectCores() )

# recommended running this is a fresh R session or restarting current session
# install.packages( "cmdstanr", repos = c( "https://mc-stan.org/r-packages/", getOption("repos") ) )
# install_cmdstan(version = "2.34.1")
# set_cmdstan_path()
#
# Important note: I was able to run the models successfully using cmdstanr version 2.34.1 on three different machines
# (one MacStudio and two MacBooks Pro, all with M1 or M2 processors), however, newer versions of cmdstanr stopped the chains
# due to the numerical integral not converging (using the same Stan code) and chains finishing unexpectedly!

library(here)
library(tidyverse)
library(bayesplot)
library(cmdstanr)
library(ggh4x)

color_scheme_set("viridisA")
theme_set( theme_bw(base_size = 12) )

source("ExGaussian_fake_data_simulation.R") # read data generating function


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
  aes(x = rt, fill = `Response type: `  ) +
  geom_density(alpha = .3) +
  labs(x = "Response time (s)", y = "Density") +
  facet_grid2(id ~ `Stimulation type: `, scales = "free_y", independent = F) +
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


## POSTERIOR PREDICTIVE CHECK ----

### ---- posterior prediction via densities ----
ppc_density <- function(data, type, preds, cols, tit) lapply(
  
  names(preds),
  function(i)
    
    cbind.data.frame(
      subset(data, id == i & signal == type),
      preds[[i]][ sample(1:4e3, 1e2), which( subset(data, id == i)$signal == type) ] %>% t()
    )
  
) %>%
  
  do.call( rbind.data.frame, . ) %>%
  pivot_longer( cols = c( "rt",as.character(1:100) ), names_to = "sample", values_to = "Response time (s)" ) %>%
  mutate( source = if_else(sample == "rt", "observed", "predicted") ) %>%
  
  ggplot() +
  aes(x = `Response time (s)`, colour = source, size = source, group = sample) +
  geom_density() +
  scale_size_manual( values = c(1.15,.15) ) +
  scale_colour_manual(values = cols) +
  facet_wrap( ~ id, ncol = 2, scales = "free" ) +
  labs(
    title = tit,
    subtitle = "Thick lines represent observed data, thin lines represent model posterior predictions"
  ) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    plot.title = element_text(hjust = .5, face = "bold"),
    plot.subtitle = element_text(hjust = .5)
  )


### ---- observed data ----

d2 <-
  subset(d1, correct != "missed") %>%
  mutate( signal = if_else(signal == "signal", 1, 0) ) %>%
  filter( !( rt < ssd & signal == 1 & !is.na(rt) ) )


### ---- parameters for control condition ----

# prepare model fits
fit_con <- lapply(
  
  setNames( names(fit_many_indi), names(fit_many_indi) ),
  function(x)
    fit_many_indi[[x]]$ctrl
  
)

# extract estimates
estCON <- lapply(
  
  names(fit_con),
  function(i)
    
    fit_con[[i]]$draws(format = "data.frame") %>%
    select( starts_with("Int"), .chain, .iteration ) %>%
    mutate(ID = i, .before = 1)
  
) %>%
  
  do.call( rbind.data.frame, . ) %>%
  pivot_longer( starts_with("Int"), values_to = "value", names_to = "name" ) %>%
  mutate( name = sub("Int_", "", ( sub("_0","",name) ) ) ) %>%
  mutate( type = sub(".*_", "", name), parameter = sub("_.*", "", name), ID = as.character(ID) )

# extract posterior predictions
ppredCON <- lapply(
  
  setNames( names(fit_con), names(fit_con) ), # one for each participant
  function(i) {
    
    # extract posterior draws for subject i
    df <-
      subset(estCON, ID == i) %>%
      select(.chain, .iteration, name, value) %>%
      pivot_wider(names_from = name, values_from = value)
    
    sapply(
      
      1:nrow(df),
      function(j) {
        
        print( paste0("Participant ",i,", sample #",j) )
        
        return( with(
          
          df,
          ssrt_data_sim(
            alpha_go = c(mu_go[j], 0),
            alpha_stop = c(mu_stop[j], 0),
            beta_go = c(sigma_go[j], 0),
            beta_stop = c(sigma_stop[j], 0),
            gamma_go = c(lambda_go[j], 0),
            gamma_stop = c(lambda_stop[j], 0),
            tau_go = c(0, 0),
            tau_stop = c(0, 0),
            zeta_go = c(0, 0),
            zeta_stop = c(0, 0),
            epsilon_go = c(0, 0),
            epsilon_stop = c(0, 0),
            N = 1,
            df = subset(d2, id == i & cond == "ctrl") %>% mutate(id = 1, trial = 1:nrow(.), rt = NA)
          )$data$rt
          
        ) )
        
      }
    ) %>% t()

  }
)


### ---- parameters for experimental condition ----

# prepare model fits
fit_exp <- lapply(
  
  setNames( names(fit_many_indi), names(fit_many_indi) ),
  function(x)
    fit_many_indi[[x]]$exp
  
)

# extract estimates
estEXP <- lapply(
  
  names(fit_exp),
  function(i)
    
    fit_exp[[i]]$draws(format = "data.frame") %>%
    select( starts_with("Int"), .chain, .iteration ) %>%
    mutate(ID = i, .before = 1)
  
) %>%
  
  do.call( rbind.data.frame, . ) %>%
  pivot_longer( starts_with("Int"), values_to = "value", names_to = "name" ) %>%
  mutate( name = sub("Int_", "", ( sub("_0","",name) ) ) ) %>%
  mutate( type = sub(".*_", "", name), parameter = sub("_.*", "", name), ID = as.character(ID) )

# extract posterior predictions
ppredEXP <- lapply(
  
  setNames( names(fit_exp), names(fit_exp) ), # one for each participant
  function(i) {
    
    # extract posterior draws for subject i
    df <-
      subset(estEXP, ID == i) %>%
      select(.chain, .iteration, name, value) %>%
      pivot_wider(names_from = name, values_from = value)
    
    sapply(
      
      1:nrow(df),
      function(j) {
        
        print( paste0("Participant ",i,", sample #",j) )
        
        return( with(
          
          df,
          ssrt_data_sim(
            alpha_go = c(mu_go[j], 0),
            alpha_stop = c(mu_stop[j], 0),
            beta_go = c(sigma_go[j], 0),
            beta_stop = c(sigma_stop[j], 0),
            gamma_go = c(lambda_go[j], 0),
            gamma_stop = c(lambda_stop[j], 0),
            tau_go = c(0, 0),
            tau_stop = c(0, 0),
            zeta_go = c(0, 0),
            zeta_stop = c(0, 0),
            epsilon_go = c(0, 0),
            epsilon_stop = c(0, 0),
            N = 1,
            df = subset(d2, id == i & cond == "exp") %>% mutate(id = 1, trial = 1:nrow(.), rt = NA)
          )$data$rt

        ) )
 
      }

    ) %>% t()
  }
)


### ---- show densities ----

ppc_density(subset(d2, cond == "ctrl"), type = 0, preds = ppredCON, cols = c("red4","lightpink"), tit = "GO TRIALS")
ppc_density(subset(d2, cond == "ctrl"), type = 1, preds = ppredCON, cols = c("blue4","lightblue"), tit = "SIGNAL-RESPOND TRIALS")
ppc_density(subset(d2, cond == "exp"), type = 0, preds = ppredEXP, cols = c("red4","lightpink"), tit = "GO TRIALS")
ppc_density(subset(d2, cond == "exp"), type = 1, preds = ppredEXP, cols = c("blue4","lightblue"), tit = "SIGNAL-RESPOND TRIALS")
