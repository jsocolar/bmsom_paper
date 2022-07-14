# Script to fit bMSOM and data-augmented versions of occupancy model to West
# Andes dataset

# packages
library(cmdstanr); library(flocker)

# fit the bMSOM
## data & model
data_package <- readRDS("/Users/jacobsocolar/Dropbox/Work/Occupancy/biogeographicMSOM/wandes_data/wandes_bsd9_package.RDS")
wandes_mod <- cmdstan_model("/Users/jacobsocolar/Dropbox/Work/Code/Occupancy/biogeographicMSOM/stan_files/occupancy_v9_wandes.stan", 
                            force_recompile = T)

## sample
wandes_samples <- wandes_mod$sample(data = data_package$data, 
                               chains = 4,
                               parallel_chains = 4,
                               #threads_per_chain = 1,
                               refresh = 1,
                               iter_sampling = 1000,
                               iter_warmup = 1000,
                               save_warmup = T,
                               step_size = .0015,
                               max_treedepth = 9,
                               output_dir = "/Users/jacobsocolar/Dropbox/Work/Occupancy/biogeographicMSOM/wandes_stan_outputs",
                               adapt_engaged = T)


# attempt to fit the data-augmented MSOM
# First attempt
fd <- readRDS("/Users/Jacob/Desktop/fd.RDS")

user_prior <- c(brms::set_prior("normal(-3, 1)", class = "Intercept"), 
                brms::set_prior("logistic(0,1)", class = "Intercept", dpar = "Omega"),
                brms::set_prior("normal(-7, 2.5)", class = "Intercept", dpar = "occ"), 
                brms::set_prior("normal(0, 5)", coef = "elev", dpar = "occ"),
                brms::set_prior("normal(0, 5)", coef = "elev2", dpar = "occ"),
                brms::set_prior("normal(0, 1)", coef = "pasture", dpar = "occ"),
                brms::set_prior("normal(0, .75)", coef = "pasture"),
                brms::set_prior("normal(0, .5)", coef = "time"),
                brms::set_prior("normal(0, .25)", coef = "obsDE"),
                brms::set_prior("normal(0, 3)", class = "sd", group = "cluster:species", dpar = "occ"),
                brms::set_prior("normal(0, 2)", class = "sd", group = "species", coef = "Intercept", dpar = "occ"),
                brms::set_prior("normal(0, 3)", class = "sd", group = "species", coef = "elev", dpar = "occ"),
                brms::set_prior("normal(0, 3)", class = "sd", group = "species", coef = "elev2", dpar = "occ"),
                brms::set_prior("normal(0, 1)", class = "sd", group = "species", coef = "pasture", dpar = "occ"),
                brms::set_prior("normal(0, 2)", class = "sd", group = "species", coef = "Intercept"),
                brms::set_prior("normal(0, 1)", class = "sd", group = "species", coef = "pasture"),
                brms::set_prior("normal(0, 1)", class = "sd", group = "species", coef = "time")
)

da_fit <- 
  flock(f_occ = ~ elev + elev2 + pasture + 
          (1 |int| species) + 
          (0 + elev + elev2 + pasture || species) +
          (1 | cluster:species),
        f_det = ~ pasture + time + obsDE + 
          (1 |int| species) +
          (0 + pasture + time || species),
        flocker_data = fd, 
        prior = user_prior,
        backend = "cmdstanr",
        cores = 1, chains = 1, refresh = 1, iter = 2)
saveRDS(da_fit, "/Users/Jacob/Desktop/da_fit.RDS")

# An abortive run of the following model, which generated the 

