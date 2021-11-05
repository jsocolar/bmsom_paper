# Run full BBS occupancy analysis
setwd("/Users/jacob/Dropbox/Work/Occupancy/biogeographicMSOM")
library("brms")

# read data
flattened_data <- readRDS('flattened_data_2018.RDS')
flattened_data$elev_scaled <- as.numeric(scale(flattened_data$elev))
flattened_data$distance_transformed <- 
  boot::inv.logit((flattened_data$distance)/200000)

# demonstrate distance transform
dist_bins <- seq(min(flattened_data$distance), max(flattened_data$distance), length.out = 200)
dist_bin_midpoints <- dist_bins[-length(dist_bins)] + diff(dist_bins)/2
p <- vector()
for(i in seq_along(dist_bin_midpoints)){
  p[i] <- mean(flattened_data$Q[flattened_data$distance > dist_bins[i] & 
                                  flattened_data$distance < dist_bins[i+1]])
}
plot(boot::logit(p) ~ dist_bin_midpoints, ylab = "logit proportion", xlab = "distance")
funnyroot <- function(x, n){
  x_pos <- 2*as.numeric(x > 0) - 1
  out <- x_pos * abs(x)^n
  return(out)
}

#dbm3 <- exp(funnyroot(dist_bin_midpoints/500000, .7))
dbm3 <- boot::inv.logit((dist_bin_midpoints)/200000)
dbm3 <- dist_bin_midpoints
dbm3[dbm3 > 0] <- 6*dbm3[dbm3>0]
plot(boot::logit(p) ~ dbm3, ylab = "logit proportion", xlab = "transformed distance")
plot(boot::logit(p) ~ dist_bin_midpoints, ylab = "logit proportion", xlab = "transformed distance")
plot(p ~ dist_bin_midpoints, ylab = "logit proportion", xlab = "transformed distance")



plot(boot::logit(p)[dbm3 < .1] ~ dbm3[dbm3 < .1], ylab = "logit proportion", xlab = "transformed distance")

fd3 <- flattened_data[,c("N", "species", "site", "elev_scaled", "distance", "distance_transformed")]
fd3$elev_scaled2 <- fd3$elev_scaled^2
fd3$trials <- 5

prior1 <- c(set_prior("logistic(0, 1)", class = "Intercept"), 
            set_prior("normal(0, 10)", class = "sd"), 
            set_prior("normal(0, 10)", class = "sd", dpar = "zi"), 
            set_prior("normal(0, 10)", class = "b", dpar = "zi"))

# traditional MSOM
bbs_naive <- brm(bf(N | trials(trials) ~ (1 |b| species), 
                 zi ~ elev_scaled +
                   (1 |b| species) + 
                   (0 + elev_scaled || species)), 
                 data = fd3, family = zero_inflated_binomial(), 
                 prior = prior1, cores = 4, backend = 'cmdstanr')
saveRDS(bbs_naive, "/Users/jacob/Dropbox/Work/Occupancy/biogeographicMSOM/stan_outputs/brms_bbs/naive_final.RDS")

# bMSOM
bbs_dist <- brm(bf(N | trials(trials) ~ (1 |b| species), 
                   zi ~ elev_scaled + distance_transformed + 
                     (1 |b| species) + 
                     (0 + elev_scaled || species) +
                     (0 + distance_transformed || species)), 
                data = fd3, family = zero_inflated_binomial(), 
                prior = prior1, backend = 'cmdstanr', cores = 4,
                term_buffer = 100)
saveRDS(bbs_dist, "/Users/jacob/Dropbox/Work/Occupancy/biogeographicMSOM/stan_outputs/brms_bbs/dist_final.RDS")

# Clipped bMSOM
fd4 <- fd3[fd3$distance < 400000, ]
bbs_dist_clip <- brm(bf(N | trials(trials) ~ (1 |b| species), 
                   zi ~ elev_scaled + distance_transformed + 
                     (1 |b| species) + 
                     (0 + elev_scaled || species) +
                     (0 + distance_transformed || species)), 
                data = fd4, family = zero_inflated_binomial(), 
                prior = prior1, backend = 'cmdstanr', cores = 4)
saveRDS(bbs_dist_clip, "/Users/jacob/Dropbox/Work/Occupancy/biogeographicMSOM/stan_outputs/brms_bbs/dist_clip_400K_final.RDS")
