# This script is used to simulate SSRT data from the assumed data-generating process
# and fit Stan models with the same structure on the data

rm( list = ls() ) # clear environment
options( mc.cores = parallel::detectCores() ) # set-up multiple cores

# recommended running this is a fresh R session or restarting current session
# install.packages( "cmdstanr", repos = c( "https://mc-stan.org/r-packages/", getOption("repos") ) )

library(here)
library(tidyverse)
library(bayesplot)
library(cmdstanr)
library(ggh4x)

color_scheme_set("viridisA")
theme_set( theme_bw() )

source("ExGaussian_fake_data_simulation.R") # read data generating function


# UTILS ----

## plotting recovery results (posteriors vs true data-generating value) ----
reco_plot <- function(data, col = "red", tit = "Go Runner") data %>%    

  ggplot() +
  aes(x = value) +
  geom_histogram(fill = "grey50") +
  geom_vline( aes( xintercept = data$true_value ), colour = col, linewidth = 1.2 ) +
  facet_grid2(ID ~ parameter, scales = "free", independent = T) +
  labs(title = tit) +
  theme( plot.title = element_text(face = "bold", size = 14, hjust = .5) )


## posterior prediction for densities ----
ppc_density <- function(data, preds, cols, tit) lapply(
  
  1:length(preds),
  function(i)
    
    cbind.data.frame(
      subset(data, id == i),
      preds[[i]][ sample(1:4e3, 1e2), subset(data, id == i)$trial ] %>% t()
    )
  
) %>%
  
  do.call( rbind.data.frame, . ) %>%
  pivot_longer( cols = c( "rt",as.character(1:100) ), names_to = "sample", values_to = "Response time (s)" ) %>%
  mutate( source = if_else(sample == "rt", "observed", "predicted"), id = paste0("participant #",id) ) %>%
  
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
  
  
# GENERATE DATA ----

# use default priors and default experiment length for 8 subjects
d0 <- ssrt_data_sim(N = 8, S = 1:18)

## SOME SANITY CHECKS ----

# check for negative response times
summary(d0$data$rt)
sum(d0$data$rt < 0, na.rm = T) / nrow(d0$data)

# check success rates in STOP-SIGNAL trials
# should be approximately 50/50 conditional on subject
t( sapply( unique(d0$data$id), function(i) table( subset(d0$data, signal == 1 & id == i)$response ) ) )

# GO response times should be generally slower than STOP-RESPOND response times
# under the assumption of the horse race model
d0$data %>%
  
  # some formatting shinanigans
  filter( complete.cases(rt) ) %>%
  mutate( `Response type: ` = ifelse( signal == "0", "GO", "SIGNAL-RESPOND" ) ) %>%
  mutate( subject = paste0("Subject #",id) ) %>%
  
  # plot it
  ggplot() +
  aes(x = rt, fill = `Response type: `) +
  geom_density(alpha = .5) +
  labs(x = "Response time (s)", y = "Density") +
  facet_wrap(~subject, ncol = 2) +
  theme(legend.position = "bottom")


# FITTING ----

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


## single model, many participants ----

# loop through d0 data and fit the same individual model to each participant
# extract number of participants
k <- length( unique(d0$data$id) )

# prepare a list for the fits
fit_many_indi <- list()

# fit it
for ( i in 1:k ) {
  
  # data
  dGO <- subset(d0$data, id == i) %>% filter( signal == 0 )
  dSR <- subset(d0$data, id == i) %>% filter( signal == 1 & !is.na(rt) )
  dNA <- subset(d0$data, id == i) %>% filter( signal == 1 & is.na(rt) )
  
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
  fit_many_indi[[i]] <- mod_indi$sample(data = dlist, chains = 4, save_warmup = T, init = ifun, seed = 87542)
  
}

# plot trace plots, each for ca 3 s
for (i in 1:k) {

  print( mcmc_trace( fit_many_indi[[i]]$draws(inc_warmup = T), n_warmup = 1e3 ) )
  Sys.sleep(3)

}

# save the results for immediate use if needed/wanted
saveRDS( object = fit_many_indi, file = here("fake_data_fits.rds") )


# RECOVERY CHECKS ----

# extract data-generating values
real <- with(
  
  d0$parameters,
  lapply(
    
    setNames( c("mu","sigma","lambda"), c("mu","sigma","lambda") ),
    function(i)
      
      sapply(
        1:k,
        function(j)
          c( go = get(i)[type == "global intercept" & racer == "go"] + get(i)[type == "subject-level variability" & racer == "go"] * get(i)[type == "varying effect" & racer == "go" & id == j],
             stop = get(i)[type == "global intercept" & racer == "stop"] + get(i)[type == "subject-level variability" & racer == "stop"] * get(i)[type == "varying effect" & racer == "stop" & id == j]
          )
        
      ) %>%
      
      as.data.frame() %>%
      rownames_to_column("type") %>%
      mutate( parameter = i, .before = 1 )
  )
) %>%
  
  do.call( rbind.data.frame, . ) %>%
  pivot_longer( cols = starts_with("V"), names_to = "ID", values_to = "true_value" ) %>%
  mutate( ID = sub("V","",ID) )

# extract models' posteriors
estimated <- lapply(
  
  1:k,
  function(i)
    
    fit_many_indi[[i]]$draws(format = "data.frame") %>%
    select( starts_with("Int"), .chain, .iteration ) %>%
    mutate(ID = i, .before = 1)
  
) %>%

  do.call( rbind.data.frame, . ) %>%
  pivot_longer( starts_with("Int"), values_to = "value", names_to = "name" ) %>%
  mutate( name = sub("Int_", "", ( sub("_0","",name) ) ) ) %>%
  mutate( type = sub(".*_", "", name), parameter = sub("_.*", "", name), ID = as.character(ID) )

# prepare a single data set with true data generating values and posteriors alike
pars <-
  left_join( estimated, real, by = c("ID","type","parameter") ) %>%
  mutate( parameter = factor(parameter, levels = c("mu","sigma","lambda"), ordered = T) )

# check recovery of the parameters across participants
reco_plot(subset(pars, type == "go"), col = "red", tit = "Go Runner") # GO runner
reco_plot(subset(pars, type == "stop"), col = "blue", tit = "Stop Runner") #STOP runner


# POSTERIOR PREDICTION ----

# extract a set of posterior predictions
ppred <- lapply(
  
  1:k, # one for each participant
  function(i) {
    
    # extract posterior draws for subject i
    df <-
      subset(estimated, ID == i) %>%
      select(.chain, .iteration, name, value) %>%
      pivot_wider(names_from = name, values_from = value)
    
    sapply(
      
      1:nrow(df),
      function(j) {
        
        print( paste0("Participant #",i,", sample #",j) )
        
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
              df = subset(d0$data, id == i) %>% mutate(id = 1, rt = NA)
            )$data$rt
          
          ) )

      }
      
    ) %>% t()
    
    
  }
)

# print posterior predictions (density overlaps)
# not using bayesplot because its function cannot use vectors with NAs
# in the signal-respond trials, it makes little sense to "align" NAs between observed data and posterior predictions
ppc_density(subset(d0$data, signal == 0), preds = ppred, cols = c("red4","lightpink"), tit = "GO TRIALS")
ppc_density(subset(d0$data, signal == 1), preds = ppred, cols = c("blue4","lightblue"), tit = "SIGNAL-RESPOND TRIALS")

