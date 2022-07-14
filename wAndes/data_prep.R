# Script to prepare data for model fitting

wandes_birds <- readRDS("/Users/jacobsocolar/Dropbox/Work/Occupancy/biogeographicMSOM/wandes_birds.RDS")

# clip to elevations 
wandes_birds_2 <- wandes_birds[wandes_birds$elev_sp_standard > 0 & 
                               wandes_birds$elev_sp_standard < 1, ]

# Explore good functional form for distance covariate:
hist(wandes_birds_2$distance_from_range)
hist(wandes_birds_2$distance_from_range[wandes_birds_2$Q==1])
quantile(wandes_birds_2$distance_from_range[wandes_birds_2$Q==1], .05)
quantile(wandes_birds_2$distance_from_range[wandes_birds_2$Q==1], .95)
min(wandes_birds_2$distance_from_range)
max(wandes_birds_2$distance_from_range)
breaks <- seq(from = -600000, to = 160000, length.out = 600)
midpoints <- breaks[-length(breaks)] + (breaks[2]-breaks[1])/2

n0 <- n1 <- vector()
for(i in 1:length(midpoints)){
  n0[i] <- sum(wandes_birds_2$distance_from_range > breaks[i] & wandes_birds_2$distance_from_range < breaks[i+1] & wandes_birds_2$Q == 0)
  n1[i] <- sum(wandes_birds_2$distance_from_range > breaks[i] & wandes_birds_2$distance_from_range < breaks[i+1] & wandes_birds_2$Q == 1)
}
range_data <- data.frame(prop_det = n1/(n1+n0), distance = midpoints)
range_data$logit_prop <- boot::logit(range_data$prop_det)
plot(logit_prop ~ distance, data = range_data)
plot(logit_prop ~ distance, data = range_data[n1>1,])

# High noise at very small distances, probably related to smallish sample sizes 
# and the idiosyncracies of just a few species that are this distant anywhere in 
# the data.
# Then a nice downwards trend near zero, and pair of obnoxious outliers well 
# out-of-range
wandes_birds_2$species[wandes_birds_2$distance_from_range > 60000 & wandes_birds_2$Q == 1]
# These outliers are entirely confined to a single species, which simply will 
# have a very atypical slope for the range covariate in the model. We'll focus 
# on getting something that fits nicely between -60000 and 60000

breaks <- seq(from = -60000, to = 60000, length.out = 100)
midpoints <- breaks[-length(breaks)] + (breaks[2]-breaks[1])/2

n0 <- n1 <- vector()
for(i in 1:length(midpoints)){
  n0[i] <- sum(wandes_birds_2$distance_from_range > breaks[i] & wandes_birds_2$distance_from_range < breaks[i+1] & wandes_birds_2$Q == 0)
  n1[i] <- sum(wandes_birds_2$distance_from_range > breaks[i] & wandes_birds_2$distance_from_range < breaks[i+1] & wandes_birds_2$Q == 1)
}
range_data <- data.frame(prop_det = n1/(n1+n0), distance = midpoints)
range_data$logit_prop <- boot::logit(range_data$prop_det)
plot(logit_prop ~ distance, data = range_data)
# The dip down on the left side involves just a handful of species in bins with 
# only a single detection each (high noise.  We can ignore these for now)

plot(logit_prop ~ distance, data = range_data[n1>1,])
wandes_birds_2$species[wandes_birds_2$distance_from_range > -60000 & wandes_birds_2$distance_from_range < -40000 & wandes_birds_2$Q == 1]

funnyroot <- function(x, n){
  x_pos <- 2*as.numeric(x > 0) - 1
  out <- x_pos * abs(x)^n
  return(out)
}

plot(logit_prop ~ distance, data = range_data)

# Other than this dip, we already look linear-ish, but perhaps we're a bit concave down.
# We have a strong a priori expectation for concave-downness and an asymptote at low values
# so we'll start by exponentiating in search of a formulation that places hierarchical priors 
# on a common "core-of-range" scenario

range_data$dist_trans <- exp(range_data$distance)
plot(logit_prop ~ dist_trans, data = range_data)

# Rescale the units
range_data$dist_trans <- exp(range_data$distance/10000)
plot(logit_prop ~ dist_trans, data = range_data)

range_data$dist_trans <- exp(range_data$distance/20000)
plot(logit_prop ~ dist_trans, data = range_data)

range_data$dist_trans <- exp(range_data$distance/30000)
plot(logit_prop ~ dist_trans, data = range_data)

range_data$dist_trans <- exp(range_data$distance/50000)
plot(logit_prop ~ dist_trans, data = range_data)

# By the time we start to look reasonable, we've pulled the minimum values way too far above zero.

range_data$dist_trans <- exp(funnyroot(range_data$distance/10000, .5))
plot(logit_prop ~ dist_trans, data = range_data)

# This is looking much better. We can't completely fix those two points on the far right,
# but this seems like an ok compromise.

range_data$dist_trans <- boot::inv.logit(range_data$distance/14900)
plot(logit_prop ~ dist_trans, data = range_data)
plot(logit_prop ~ distance, data = range_data)

wandes_birds$dist_trans <- boot::inv.logit(wandes_birds$distance_from_range/10000)

# function to translate effects-coded predictor to 1/2 indices
grt_b <- function(predictor){
  if(sum(predictor != 1 & predictor != -1) != 0){stop('predictor must be effects coded')}
  predictor[predictor == 1] <- 2
  predictor[predictor == -1] <- 1
  return(predictor)
}
ec <- c(-1,1)

# create data for stan model
wandes_bsd9 <- list(
  # Grainsize for reduce_sum
  grainsize = 1,
  # Dimensions
  # Random effects
  n_spCl = length(unique(wandes_birds$sp_cl)),
  n_spSr = length(unique(wandes_birds$sp_sr)),
  n_sp = length(unique(wandes_birds$species)),
  n_fam = length(unique(wandes_birds$Family)),
  n_spObs = length(unique(as.vector(wandes_birds$sp_obs_matrix[!is.na(as.vector(wandes_birds$sp_obs_matrix))]))),
  # Dataset size
  n_tot = nrow(wandes_birds),
  n_visit_max = max(wandes_birds$nv),
  # Integer data matrix
  integer_data = data.frame(
    # Detection data
    det_data = wandes_birds$det_data,  # 1-4
    # Q and nv
    Q = wandes_birds$Q,               # 5
    nv = wandes_birds$nv,             # 6
    # Random effect IDs
    id_spCl = as.numeric(as.factor(wandes_birds$sp_cl)),     # 7
    id_spSr = as.numeric(as.factor(wandes_birds$sp_sr)),     # 8
    id_sp = as.numeric(as.factor(wandes_birds$species)),     # 9
    id_fam = as.numeric(as.factor(wandes_birds$Family)),     # 10
    id_spObs = wandes_birds$sp_obs_matrix,                   # 11-14
    # Covariate indices
    lowland = grt_b(ec[wandes_birds$lowland+1]),             # 15
    pasture = grt_b(ec[wandes_birds$pasture+1]),                          # 16
    mountain_barrier = grt_b(ec[wandes_birds$mountain_limited+1]),        # 17
    valley_barrier = grt_b(ec[wandes_birds$valley_limited+1]),            # 18
    forestPresent = grt_b(ec[wandes_birds$forest_present+1]),             # 19
    forestSpecialist = grt_b(ec[wandes_birds$forest_specialist+1]),       # 20
    tfSpecialist = grt_b(ec[wandes_birds$tf_specialist+1]),               # 21
    dryForestPresent = grt_b(ec[wandes_birds$dry_forest_present+1]),      # 22
    floodDrySpecialist = grt_b(ec[wandes_birds$flood_dry_specialist+1]),  # 23
    aridPresent = grt_b(ec[wandes_birds$arid_present+1]),                 # 24
    migratory = grt_b(ec[as.numeric(!is.na(wandes_birds$start1))+1]),     # 25
    dietInvert = grt_b(ec[as.numeric(wandes_birds$Diet.5Cat == "Invertebrate")+1]),    # 26
    dietCarn = grt_b(ec[as.numeric(wandes_birds$Diet.5Cat == "VertFishScav")+1]),      # 27
    dietFruitNect = grt_b(ec[as.numeric(wandes_birds$Diet.5Cat == "FruiNect")+1]),     # 28
    dietGran = grt_b(ec[as.numeric(wandes_birds$Diet.5Cat == "PlantSeed")+1]),         # 29
    mountainBarrier_x_pasture = grt_b(ec[wandes_birds$mountain_limited + 1] * ec[wandes_birds$pasture + 1]),        # 30
    valleyBarrier_x_pasture = grt_b(ec[wandes_birds$valley_limited + 1] * ec[wandes_birds$pasture + 1]),            # 31
    forestPresent_x_pasture = grt_b(ec[wandes_birds$forest_present + 1] * ec[wandes_birds$pasture + 1]),            # 32
    forestSpecialist_x_pasture = grt_b(ec[wandes_birds$forest_specialist + 1] * ec[wandes_birds$pasture + 1]),      # 33
    tfSpecialist_x_pasture = grt_b(ec[wandes_birds$tf_specialist + 1] * ec[wandes_birds$pasture + 1]),              # 34
    dryForestPresent_x_pasture = grt_b(ec[wandes_birds$dry_forest_present + 1] * ec[wandes_birds$pasture + 1]),     # 35
    floodDrySpecialist_x_pasture = grt_b(ec[wandes_birds$flood_dry_specialist + 1] * ec[wandes_birds$pasture + 1]), # 36
    aridPresent_x_pasture = grt_b(ec[wandes_birds$arid_present + 1] * ec[wandes_birds$pasture + 1]),                # 37
    migratory_x_pasture = grt_b(ec[as.numeric(!is.na(wandes_birds$start1))+1] * ec[wandes_birds$pasture + 1]),      # 38
    dietInvert_x_pasture = grt_b(ec[as.numeric(wandes_birds$Diet.5Cat == "Invertebrate")+1] * ec[wandes_birds$pasture + 1]),   # 39
    dietCarn_x_pasture = grt_b(ec[as.numeric(wandes_birds$Diet.5Cat == "VertFishScav")+1] * ec[wandes_birds$pasture + 1]),     # 40
    dietFruitNect_x_pasture = grt_b(ec[as.numeric(wandes_birds$Diet.5Cat == "FruiNect")+1] * ec[wandes_birds$pasture + 1]),    # 41
    dietGran_x_pasture = grt_b(ec[as.numeric(wandes_birds$Diet.5Cat == "PlantSeed")+1] * ec[wandes_birds$pasture + 1]),        # 42
    obsSM = matrix(grt_b(ec[wandes_birds$obsSM+1]), ncol=4),   # 43-46
    obsJG = matrix(grt_b(ec[wandes_birds$obsJG+1]), ncol=4),   # 47-50
    obsDE = matrix(grt_b(ec[wandes_birds$obsDE+1]), ncol=4)    # 51-54
  ),
  # Continuous distance-to-range (does not compress well)
  distance_to_range = as.vector(wandes_birds$dist_trans),
  # continuous covariates
  relev = wandes_birds$relev,
  relev2 = wandes_birds$relev2,
  lowland_x_relev = ec[wandes_birds$lowland+1] * wandes_birds$relev,
  lowland_x_relev2 = ec[wandes_birds$lowland+1] * wandes_birds$relev2,
  elevMedian = wandes_birds$elev_median_scaled,
  elevBreadth = wandes_birds$elev_breadth_scaled,
  mass = wandes_birds$log_mass_scaled,
  elevMedian_x_forestPresent = wandes_birds$elev_median_scaled * ec[wandes_birds$forest_present + 1],
  elevMedian_x_forestSpecialist = wandes_birds$elev_median_scaled * ec[wandes_birds$forest_specialist + 1],
  elevMedian_x_pasture = wandes_birds$elev_median_scaled * ec[wandes_birds$pasture + 1],
  elevBreadth_x_pasture = wandes_birds$elev_breadth_scaled * ec[wandes_birds$pasture + 1],
  mass_x_pasture = wandes_birds$log_mass_scaled * ec[wandes_birds$pasture + 1],
  elevMedian_x_forestPresent_x_pasture = wandes_birds$elev_median_scaled * ec[wandes_birds$forest_present + 1] * ec[wandes_birds$pasture + 1],
  elevMedian_x_forestSpecialist_x_pasture = wandes_birds$elev_median_scaled * ec[wandes_birds$forest_specialist + 1] * ec[wandes_birds$pasture + 1],
  time = wandes_birds$time,
  time_x_elev = sweep(wandes_birds$time, MARGIN = 1, wandes_birds$elev_median_scaled, FUN = `*`)
)

wandes_bsd9_means_and_sds <- list(time_mean = mean(c(wandes_birds$hps1, wandes_birds$hps2, wandes_birds$hps3, wandes_birds$hps4), na.rm = T), 
                                  time_sd =  sd(c(wandes_birds$hps1, wandes_birds$hps2, wandes_birds$hps3, wandes_birds$hps4), na.rm = T),
                                  relev_offset = .5, relev_sd = sd(wandes_birds$elev_sp_standard),
                                  elev_median_mean =  mean(wandes_birds$elev_median), elev_median_sd = sd(wandes_birds$elev_median),
                                  elev_breadth_mean = mean(wandes_birds$elev_breadth), elev_breadth_sd = sd(wandes_birds$elev_breadth),
                                  log_mass_mean = mean(log(wandes_birds$BodyMass.Value)), log_mass_sd = sd(log(wandes_birds$BodyMass.Value)), 
                                  distance_to_range_denom = 10000)
wandes_bsd9_package <- list(data = wandes_bsd9,
                            means_and_sds = wandes_bsd9_means_and_sds)
saveRDS(wandes_bsd9_package, "/Users/jacobsocolar/Dropbox/Work/Occupancy/biogeographicMSOM/wandes_data/wandes_bsd9_package.RDS")

# create data that plays nicely with the posterior predictions script
wandes_bsd6 <- list(
  # Grainsize for reduce_sum
  grainsize = 1,
  # Dimensions
  n_spCl = length(unique(wandes_birds$sp_cl)),
  n_spSr = length(unique(wandes_birds$sp_sr)),
  n_sp = length(unique(wandes_birds$species)),
  n_fam = length(unique(wandes_birds$Family)),
  n_spObs = length(unique(as.vector(wandes_birds$sp_obs_matrix[!is.na(as.vector(wandes_birds$sp_obs_matrix))]))),
  n_tot = nrow(wandes_birds),
  n_visit_max = max(wandes_birds$nv),
  # Detection matrix
  det_data = wandes_birds$det_data,
  # Q and nv
  Q = wandes_birds$Q,
  nv = wandes_birds$nv,
  # Random effect IDs
  id_spCl = as.numeric(as.factor(wandes_birds$sp_cl)),
  id_spSr = as.numeric(as.factor(wandes_birds$sp_sr)),
  id_sp = as.numeric(as.factor(wandes_birds$species)),
  id_fam = as.numeric(as.factor(wandes_birds$Family)),
  id_spObs = wandes_birds$sp_obs_matrix,
  # Covariates
  relev = wandes_birds$relev,
  relev2 = wandes_birds$relev2,
  lowland = ec[wandes_birds$lowland+1],
  pasture = ec[wandes_birds$pasture+1],
  mountain_barrier = ec[wandes_birds$mountain_limited+1],
  valley_barrier = ec[wandes_birds$valley_limited+1],
  elevMedian = wandes_birds$elev_median_scaled,
  elevBreadth = wandes_birds$elev_breadth_scaled,
  forestPresent = ec[wandes_birds$forest_present+1],
  forestSpecialist = ec[wandes_birds$forest_specialist+1],
  tfSpecialist = ec[wandes_birds$tf_specialist+1],
  dryForestPresent = ec[wandes_birds$dry_forest_present+1],
  floodDrySpecialist = ec[wandes_birds$flood_dry_specialist+1],
  aridPresent = ec[wandes_birds$arid_present+1],
  migratory = ec[as.numeric(!is.na(wandes_birds$start1))+1],
  mass = wandes_birds$log_mass_scaled,
  dietInvert = ec[as.numeric(wandes_birds$Diet.5Cat == "Invertebrate")+1],
  dietCarn = ec[as.numeric(wandes_birds$Diet.5Cat == "VertFishScav")+1],
  dietFruitNect = ec[as.numeric(wandes_birds$Diet.5Cat == "FruiNect")+1],
  dietGran = ec[as.numeric(wandes_birds$Diet.5Cat == "PlantSeed")+1],
  distance_to_range = as.vector(wandes_birds$dist_trans),
  time = wandes_birds$time,
  obsDE = matrix(ec[wandes_birds$obsDE+1], ncol=4)
)

wandes_bsd6_means_and_sds <- list(time_mean = mean(c(wandes_birds$hps1, wandes_birds$hps2, wandes_birds$hps3, wandes_birds$hps4), na.rm = T), 
                                     time_sd =  sd(c(wandes_birds$hps1, wandes_birds$hps2, wandes_birds$hps3, wandes_birds$hps4), na.rm = T),
                                     relev_offset = .5, relev_sd = sd(wandes_birds$elev_sp_standard),
                                     elev_median_mean =  mean(wandes_birds$elev_median), elev_median_sd = sd(wandes_birds$elev_median),
                                     elev_breadth_mean = mean(wandes_birds$elev_breadth), elev_breadth_sd = sd(wandes_birds$elev_breadth),
                                     log_mass_mean = mean(log(wandes_birds$BodyMass.Value)), log_mass_sd = sd(log(wandes_birds$BodyMass.Value)), 
                                     distance_to_range_denom = 10000)
wandes_bsd6_package <- list(data = wandes_bsd6,
                                means_and_sds = wandes_bsd6_means_and_sds)
saveRDS(wandes_bsd6_package, "/Users/jacobsocolar/Dropbox/Work/Occupancy/biogeographicMSOM/wandes_data/wandes_bsd6_package.RDS")

# Create 
wandes_species <- unique(wandes_birds$species[wandes_birds$Q == 1])
wandes_points <- unique(wandes_birds$point)
length(wandes_points) == length(unique(wandes_birds$point[wandes_birds$species %in% wandes_species]))

obs <- array(dim = c(length(wandes_points), 4, length(wandes_species)))
for (pt in 1:dim(obs)[1]) {
  print(pt)
  for (sp in 1:dim(obs)[3]) {
    obs[pt, ,sp] <- as.numeric(wandes_birds[wandes_birds$point == wandes_points[pt] &
                                 wandes_birds$species == wandes_species[sp],
                                c("v1", "v2", "v3", "v4")])
  }
}
obs[is.na(obs)] <- 0 #zero-fill species that were previously biogeographically clipped


get_unique_vals <- function(col_name) {
  out <- rep(NA, length(wandes_points))
  for (i in 1:length(wandes_points)) {
    out[i] <- unique(wandes_birds[wandes_birds$point == wandes_points[i],col_name])
  }
  return(out)
}

site_covs <- data.frame(elev = scale(get_unique_vals("elev_ALOS")))
site_covs$elev2 <- site_covs$elev^2
site_covs$cluster <- get_unique_vals("cluster")
site_covs$subregion <- get_unique_vals("subregion")
site_covs$pasture <- get_unique_vals("pasture")
site_covs$pasture[site_covs$pasture == 0] <- -1

wb2 <- wandes_birds[!duplicated(wandes_birds$point), ]
event_covs <- list()
event_covs$time <- event_covs$obsDE <- matrix(nrow = length(wandes_points), ncol = 4)
for (i in 1:length(wandes_points)) {
  event_covs$time[i, ] <- wb2$time[wb2$point == wandes_points[i]]
  event_covs$obsDE[i, ] <- wb2$obsDE[wb2$point == wandes_points[i]]
}
event_covs$obsDE[event_covs$obsDE==0] <- -1

fd <- make_flocker_data_augmented(obs = obs, n_aug = 1000, 
                                  site_covs = site_covs, event_covs = event_covs)

saveRDS(fd, "/Users/Jacob/Desktop/fd.RDS")

obs2 <- as.matrix(obs[,,1])
for(i in 2:dim(obs)[3]){
  obs2 <- rbind(obs2, as.matrix(obs[,,i]))
}
obs2 <- rbind(obs2, matrix(0, nrow = 1000*dim(obs)[1], ncol = dim(obs)[2]))

site_covs2 <- data.frame(elev=site_covs$elev, elev2 = site_covs$elev2, 
                         cluster = site_covs$cluster, 
                         subregion = site_covs$subregion,
                         pasture = site_covs$pasture,
                         species = rep(c(1:(1000+dim(obs)[3])), each = dim(obs)[1]))
event_covs2 <- event_covs
for(i in 2:(1000+dim(obs)[3])){
  event_covs2$obsDE <- rbind(event_covs2$obsDE, event_covs$obsDE)
  event_covs2$time <- rbind(event_covs2$time, event_covs$time)
}


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
da_fit <- readRDS("/Users/Jacob/Desktop/da_fit.RDS")

use_data <- brms::standata(da_fit)

use_data$N_unit <- da_fit$data$n_unit[1]

identical(unname(use_data$J_3[1:use_data$N_unit]), unname(use_data$J_3[use_data$N_unit + (1:use_data$N_unit)]))
use_data$X_occ <- use_data$X_occ[1:use_data$N_unit, ]
use_data$Z_1_occ_2 <- use_data$Z_1_occ_2[1:use_data$N_unit]
use_data$Z_3_occ_1 <- use_data$Z_3_occ_1[1:use_data$N_unit]
use_data$J_3 <- use_data$J_3[1:use_data$N_unit]
use_data$J_4 <- use_data$J_4[1:use_data$N_unit]
use_data$Z_4_occ_1 <- use_data$Z_4_occ_1[1:use_data$N_unit]
use_data$Z_4_occ_2 <- use_data$Z_4_occ_2[1:use_data$N_unit]
use_data$Z_4_occ_3 <- use_data$Z_4_occ_3[1:use_data$N_unit]

mymod <- cmdstanr::read_cmdstan_csv("/Users/jacob/Desktop/model_a1fa1c9353a014dbff6d8ad83ea10a95-202109190701-2-904c78.csv")
init1 <- posterior::as_draws_rvars(mymod$post_warmup_draws[1,1,])
init1_1 <- lapply(init1, posterior::draws_of)
names(init1_1)[1] <- "Intercept_Omega"
init1_1$Intercept_Omega <- boot::logit(.99)

inv_1 <- c(1, mymod$inv_metric[[1]])


mymod2 <- cmdstanr::read_cmdstan_csv("/Users/jacob/Desktop/model_a1fa1c9353a014dbff6d8ad83ea10a95-202109190701-1-904c78.csv")
init2 <- posterior::as_draws_rvars(mymod2$warmup_draws[900,1,])
init2_1 <- lapply(init2, posterior::draws_of)
names(init2_1)[1] <- "Intercept_Omega"
init2_1$Intercept_Omega <- boot::logit(.99)


mymod3 <- cmdstanr::read_cmdstan_csv("/Users/jacob/Desktop/model_a1fa1c9353a014dbff6d8ad83ea10a95-202109190701-3-904c78.csv")
init3 <- posterior::as_draws_rvars(mymod3$warmup_draws[380,1,])
init3_1 <- lapply(init3, posterior::draws_of)
names(init3_1)[1] <- "Intercept_Omega"
init3_1$Intercept_Omega <- boot::logit(.99)

mymod4 <- cmdstanr::read_cmdstan_csv("/Users/jacob/Desktop/model_a1fa1c9353a014dbff6d8ad83ea10a95-202109190701-4-904c78.csv")
init4 <- posterior::as_draws_rvars(mymod4$warmup_draws[900,1,])
init4_1 <- lapply(init4, posterior::draws_of)
names(init4_1)[1] <- "Intercept_Omega"
init4_1$Intercept_Omega <- boot::logit(.99)

rm(list = c("mymod4", "mymod3", "mymod2", "mymod", 
        "init4", "init3", "init2", "init1"))
gc()
data_aug <- cmdstanr::cmdstan_model("/Users/jacob/Desktop/da_aug_mod2.stan")

da_fit2 <- data_aug$sample(data = use_data, refresh = 1, save_warmup = T,
                           output_dir = "/Users/Jacob/Desktop", parallel_chains = 4,
                           inv_metric = inv_1, step_size = .02,
                           init = list(init1_1, init2_1, init3_1, init4_1), window = 50, term_buffer = 100)


fd2 <- flocker::make_flocker_data(obs = obs2, unit_covs = site_covs2,
                                  event_covs = event_covs2)
saveRDS(fd2, "/Users/jacob/Dropbox/Work/Occupancy/biogeographicMSOM/fd2.RDS")


user_prior2 <- c(brms::set_prior("logistic(0,1)", class = "Intercept"), 
                brms::set_prior("logistic(0,1)", class = "Intercept", dpar = "occ"), 
                brms::set_prior("normal(0, 5)", coef = "elev", dpar = "occ"),
                brms::set_prior("normal(0, 5)", coef = "elev2", dpar = "occ"),
                brms::set_prior("normal(0, 1)", coef = "pasture", dpar = "occ"),
                brms::set_prior("normal(0, .75)", coef = "pasture"),
                brms::set_prior("normal(0, .5)", coef = "time"),
                brms::set_prior("normal(0, .25)", coef = "obsDE"),
                #brms::set_prior("normal(0, 3)", class = "sd", group = "subregion:species", dpar = "occ"),
                brms::set_prior("normal(0, 3)", class = "sd", group = "cluster:species", dpar = "occ"),
                brms::set_prior("normal(0, 2)", class = "sd", group = "species", coef = "Intercept", dpar = "occ"),
                brms::set_prior("normal(0, 3)", class = "sd", group = "species", coef = "elev", dpar = "occ"),
                brms::set_prior("normal(0, 3)", class = "sd", group = "species", coef = "elev2", dpar = "occ"),
                brms::set_prior("normal(0, 1)", class = "sd", group = "species", coef = "pasture", dpar = "occ"),
                brms::set_prior("normal(0, 2)", class = "sd", group = "species", coef = "Intercept"),
                brms::set_prior("normal(0, 1)", class = "sd", group = "species", coef = "pasture"),
                brms::set_prior("normal(0, 1)", class = "sd", group = "species", coef = "time")
)

nda_fit <- 
  flocker::flock(f_occ = ~ elev + elev2 + pasture + 
          (1 |int| species) + 
          (0 + elev + elev2 + pasture || species) +
          (1 | cluster:species),
        f_det = ~ pasture + time + obsDE + 
          (1 |int| species) +
          (0 + pasture + time || species),
        flocker_data = fd2, 
        prior = user_prior2,
        backend = "cmdstanr",
        cores = 4, refresh = 1,
        save_warmup = F,
        output_dir = "/Users/Jacob/Desktop")



mymod <- cmdstanr::read_cmdstan_csv("/Users/jacob/Dropbox/Work/Occupancy/biogeographicMSOM/wandes_stan_outputs/penultimate_na_model/model_a1fa1c9353a014dbff6d8ad83ea10a95-202109190701-2-904c78.csv")
init1 <- posterior::as_draws_rvars(mymod$post_warmup_draws[1,1,-1])
init1_1 <- lapply(init1, posterior::draws_of)
inv_1 <- mymod$inv_metric[[1]]
rm(mymod)
gc()

mymod2 <- cmdstanr::read_cmdstan_csv("/Users/jacob/Dropbox/Work/Occupancy/biogeographicMSOM/wandes_stan_outputs/penultimate_na_model/model_a1fa1c9353a014dbff6d8ad83ea10a95-202109190701-1-904c78.csv")
init2 <- posterior::as_draws_rvars(mymod2$warmup_draws[900,1,-1])
init2_1 <- lapply(init2, posterior::draws_of)
rm(mymod2)
gc()

mymod3 <- cmdstanr::read_cmdstan_csv("/Users/jacob/Dropbox/Work/Occupancy/biogeographicMSOM/wandes_stan_outputs/penultimate_na_model/model_a1fa1c9353a014dbff6d8ad83ea10a95-202109190701-3-904c78.csv")
init3 <- posterior::as_draws_rvars(mymod3$warmup_draws[380,1,-1])
init3_1 <- lapply(init3, posterior::draws_of)
rm(mymod3)
gc()

mymod4 <- cmdstanr::read_cmdstan_csv("/Users/jacob/Dropbox/Work/Occupancy/biogeographicMSOM/wandes_stan_outputs/penultimate_na_model/model_a1fa1c9353a014dbff6d8ad83ea10a95-202109190701-4-904c78.csv")
init4 <- posterior::as_draws_rvars(mymod4$warmup_draws[900,1,-1])
init4_1 <- lapply(init4, posterior::draws_of)
rm(mymod4)
gc()


fd2 <- readRDS("/Users/jacob/Dropbox/Work/Occupancy/biogeographicMSOM/fd2.RDS")

nda_fit2 <- 
  flocker::flock(f_occ = ~ elev + elev2 + pasture + 
                   (1 |int| species) + 
                   (0 + elev + elev2 + pasture || species) +
                   (1 | cluster:species),
                 f_det = ~ pasture + time + obsDE + 
                   (1 |int| species) +
                   (0 + pasture + time || species),
                 flocker_data = fd2, 
                 prior = user_prior2,
                 backend = "cmdstanr",
                 cores = 4, refresh = 1,
                 save_warmup = F,
                 inv_metric = inv_1, step_size = .03,
                 init = list(init1_1, init2_1, init3_1, init4_1),  
                 adapt_engaged = F, iter = 1000, warmup = 0,
                 output_dir = "/Users/Jacob/Desktop")

saveRDS(nda_fit2, "/Users/Jacob/Dropbox/Work/Occupancy/biogeographicMSOM/nda_fit2.RDS")

nda_summary <- summary(nda_fit2)
saveRDS(nda_summary, "/Users/Jacob/Dropbox/Work/Occupancy/biogeographicMSOM/nda_summary.RDS")

nda_summary2 <- posterior::summarise_draws(nda_fit2$fit)

post_Z <- flocker::get_Z(nda_fit2)
saveRDS(post_Z, "/Users/Jacob/Dropbox/Work/Occupancy/biogeographicMSOM/post_Z.RDS")

elevs <- get_unique_vals("elev_ALOS")
pasture <- get_unique_vals("pasture")

pt_richness <- pt_richness_no <- data.frame(elevs = get_unique_vals("elev_ALOS"), 
                          pasture = get_unique_vals("pasture"))

get_richness_rep <- function (index) {
  pz_rep <- post_Z[index,]
  richness_rep <- rep(NA, 146)
  for(i in 1:146) {
    richness_rep[i] <- sum(pz_rep[i+(c(1:1313)*146)])
  }
  return(richness_rep)
}


get_richness_rep <- function (index) {
  pz_rep <- post_Z[index,]
  richness_rep <- rep(NA, 146)
  for(i in 1:146) {
    richness_rep[i] <- sum(pz_rep[i+(c(1:1313)*146)])
  }
  return(richness_rep)
}

get_richness_rep_unobserved_only <- function (index) {
  pz_rep <- post_Z[index,]
  richness_rep <- rep(NA, 146)
  for(i in 1:146) {
    richness_rep[i] <- sum(pz_rep[i+(c(314:1313)*146)])
  }
  return(richness_rep)
}



for(i in 1:100){
  print(i)
  pt_richness[,i+2] <- get_richness_rep(i*4)
}


for(i in 1:100){
  print(i)
  pt_richness_no[,i+2] <- get_richness_rep_unobserved_only(i*4)
}

pt_richness$post_mean <- NA
for(i in 1:146){pt_richness$post_mean[i] <- mean(as.numeric(pt_richness[i,3:102]))}

saveRDS(pt_richness, "/Users/Jacob/Dropbox/Work/Occupancy/biogeographicMSOM/pt_richness.RDS")

pt_richness_no$post_mean <- NA
for(i in 1:146){pt_richness_no$post_mean[i] <- mean(as.numeric(pt_richness_no[i,3:102]))}
saveRDS(pt_richness_no, "/Users/Jacob/Dropbox/Work/Occupancy/biogeographicMSOM/pt_richness_no.RDS")


forest_rich <- pt_richness[pt_richness$pasture == 0, ]
plot(post_mean ~ elevs, data = forest_rich, pch = 16)
lw1 <- loess(post_mean ~ elevs, data=forest_rich)
j <- order(forest_rich$elevs)
lines(forest_rich$elevs[j],lw1$fitted[j],col="red",lwd=3)


pasture_rich <- pt_richness[pt_richness$pasture == 1, ]
plot(post_mean ~ elevs, data = pasture_rich, pch = 16)
lw3 <- loess(post_mean ~ elevs, data=pasture_rich)
j <- order(pasture_rich$elevs)
lines(pasture_rich$elevs[j],lw3$fitted[j],col="red",lwd=3)

plot(seq(1300, 2500, 10), predict(lw1, seq(1300, 2500, 10)) - predict(lw3, seq(1300, 2500, 10)), 
     type = "l", lwd = 3, col = "red", ylim = c(0, 80),
     xlab = "elevation", ylab = "richness difference")
