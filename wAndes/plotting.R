library(cmdstanr)
library(posterior)

wandes <- read_cmdstan_csv(list.files("/Users/jacob/Dropbox/Work/Occupancy/biogeographicMSOM/wandes_stan_outputs", full.names = T))
wandes_draws <- as_draws_df(wandes$post_warmup_draws)

source("/Users/jacob/Dropbox/Work/Code/Occupancy/biogeographicMSOM/final_analysis/wandes/get_posterior_z_wandes.R")
source("/Users/jacob/Dropbox/Work/Code/colombiaBeta/bird_analysis_plotting/posterior_predictive_checks/discrepancy_functions.R")
bird_data <- readRDS("/Users/jacob/Dropbox/Work/Occupancy/biogeographicMSOM/wandes_data/wandes_bsd6_package.RDS")


# create z_info object for computing posterior Z (see get_posterior_z.R)
z_info <- data.frame(bird_data$data[8:39])
z_info$point <- wandes_birds$point
z_info$species <- wandes_birds$species
z_info$cl_q_real <- cluster_q(z_info, z_info$Q)[z_info$id_spCl]


pz <- matrix(nrow = nrow(z_info), ncol = 100)

for(i in 1:100){
  print(i)
  pz[,i] <- get_Z_probs(wandes_draws, 60*i, z_info)$Z_prob
}  

colnames(pz) <- paste0("iter", 1:100)

elev_df <- cbind(wandes_birds[, c("point", "elev_ALOS", "Q")], pz)

species_observed <- unique(wandes_birds$species[wandes_birds$Q == 1])
species <- unique(wandes_birds$species)
species_no <- species[!(species %in% species_observed)]

elev_df_no <- elev_df[wandes_birds$species %in% species_no, ]

saveRDS(elev_df, "/Users/Jacob/Desktop/elev_df.RDS")
saveRDS(elev_df_no, "/Users/Jacob/Desktop/elev_df_no.RDS")


df_forest <- elev_df[wandes_birds$pasture == 0, ]

obs_rich <- aggregate(df_forest[,3], list(df_forest$point), sum)
mod_rich <- rowMeans(aggregate(df_forest[,4:103], list(df_forest$point), sum)[,2:101])

elevs_forest <- aggregate(df_forest$elev_ALOS, list(df_forest$point), unique)

forest_df <- data.frame(elevation = elevs_forest$x, mod_rich = mod_rich, obs_rich = obs_rich$x)

plot(mod_rich ~ elevation, data = forest_df, pch = 16)
lw1 <- loess(mod_rich ~ elevation, data=forest_df)
j <- order(forest_df$elevation)
lines(forest_df$elevation[j],lw1$fitted[j],col="red",lwd=3)

plot(obs_rich ~ elevation, data = forest_df, pch = 16)
lw2 <- loess(obs_rich ~ elevation, data=forest_df)
lines(forest_df$elevation[j],lw2$fitted[j],col="red",lwd=3)


df_pasture <- elev_df[wandes_birds$pasture == 1, ]

obs_rich <- aggregate(df_pasture[,3], list(df_pasture$point), sum)
mod_rich <- rowMeans(aggregate(df_pasture[,4:103], list(df_pasture$point), sum)[,2:101])

elevs_pasture <- aggregate(df_pasture$elev_ALOS, list(df_pasture$point), unique)

pasture_df <- data.frame(elevation = elevs_pasture$x, mod_rich = mod_rich, obs_rich = obs_rich$x)

plot(mod_rich ~ elevation, data = pasture_df, pch = 16)
lw3 <- loess(mod_rich ~ elevation, data=pasture_df)
j <- order(pasture_df$elevation)
lines(pasture_df$elevation[j],lw3$fitted[j],col="red",lwd=3)


plot(obs_rich ~ elevation, data = pasture_df, pch = 16)
lw4 <- loess(obs_rich ~ elevation, data=pasture_df)
lines(pasture_df$elevation[j],lw4$fitted[j],col="red",lwd=3)


pt_richness <- readRDS("/Users/Jacob/Dropbox/Work/Occupancy/biogeographicMSOM/pt_richness.RDS")
forest_rich <- pt_richness[pt_richness$pasture == 0, ]
pasture_rich <- pt_richness[pt_richness$pasture == 1, ]



dev.off()
par(mfrow = c(4,2))
par(mar = c(4,5,0.5,2))

# Forest observed richness
plot(obs_rich ~ elevation, data = forest_df, pch = 16, 
     ylab = "observed richness", xaxt = 'n', xlab = "",
     xlim = c(1300, 2700))
lw_f_obs <- loess(obs_rich ~ elevation, data=forest_df)
j <- order(forest_df$elevation)
lines(forest_df$elevation[j],lw_f_obs$fitted[j],col="red",lwd=3)

# Pasture observed richness
plot(obs_rich ~ elevation, data = pasture_df, pch = 16, ylab = "",
     xaxt = 'n', xlab = "",
     xlim = c(1250, 2550))
lw_p_obs <- loess(obs_rich ~ elevation, data=pasture_df)
j <- order(pasture_df$elevation)
lines(pasture_df$elevation[j],lw_p_obs$fitted[j],col="red",lwd=3)

# forest modeled richness bMSOM
plot(mod_rich ~ elevation, data = forest_df, pch = 16, xaxt = "n", ylab = "bMSOM richness",
     xlim = c(1300, 2700), xaxt = 'n', xlab = "")
lw_f_bmsom <- loess(mod_rich ~ elevation, data=forest_df)
j <- order(forest_df$elevation)
lines(forest_df$elevation[j],lw_f_bmsom$fitted[j],col="red",lwd=3)

# pasture modeled richness bMSOM
plot(mod_rich ~ elevation, data = pasture_df, pch = 16, xaxt = "n", xlab = "", ylab = "",
     xlim = c(1250, 2550))
lw_p_bmsom <- loess(mod_rich ~ elevation, data=pasture_df)
j <- order(pasture_df$elevation)
lines(pasture_df$elevation[j],lw_p_bmsom$fitted[j],col="red",lwd=3)

# forest modeled richness da
plot(post_mean ~ elevs, data = forest_rich, pch = 16, xlim = c(1300, 2700),
     xaxt = 'n', xlab = "", ylab = "da MSOM richness")
lw_f_da <- loess(post_mean ~ elevs, data=forest_rich)
j <- order(forest_rich$elevs)
lines(forest_rich$elevs[j],lw_f_da$fitted[j],col="red",lwd=3)

# pasture modeled richness da
plot(post_mean ~ elevs, data = pasture_rich, pch = 16, xlim = c(1250, 2550),
     xaxt = 'n', xlab = "", ylab = "")
lw_p_da <- loess(post_mean ~ elevs, data=pasture_rich)
j <- order(pasture_rich$elevs)
lines(pasture_rich$elevs[j],lw_p_da$fitted[j],col="red",lwd=3)


# forest bMSOM minus forest da
plot(seq(1300, 2700, 10), predict(lw_f_bmsom, seq(1300, 2700, 10)) - predict(lw_f_da, seq(1300, 2700, 10)), 
     type = "l", lwd = 3, col = "red", ylim = c(-25, 25), xlim = c(1300, 2700),
     xlab = "elevation", ylab = "richness difference")
for (i in 1:100){
  mod_rich <- aggregate(df_forest[,3+i], list(df_forest$point), sum)[,2]
  forest_df_i <- data.frame(elevation = elevs_forest$x, mod_rich = mod_rich)
  lw_nf <- loess(mod_rich ~ elevation, data=forest_df_i)
  
  mod_rich2 <- pt_richness[pt_richness$pasture == 0, i + 2]
  forest_df_i2 <- data.frame(elevation = pt_richness$elevs[pt_richness$pasture == 0], mod_rich = mod_rich2)
  lw_nf2 <- loess(mod_rich ~ elevation, data=forest_df_i2)
  
  lines(seq(1300, 2700, 10), predict(lw_nf, seq(1300, 2700, 10)) - predict(lw_nf2, seq(1300, 2700, 10)), type = "l")
  
}
lines(seq(1300, 2700, 10), predict(lw_f_bmsom, seq(1300, 2700, 10)) - predict(lw_f_da, seq(1300, 2700, 10)), type = "l", lwd = 3, col = "red")


# pasture bMSOM minus forest da
plot(seq(1250, 2550, 10), predict(lw_p_bmsom, seq(1250, 2550, 10)) - predict(lw_p_da, seq(1250, 2550, 10)), 
     type = "l", lwd = 3, col = "red", ylim = c(-40, 40),
     xlab = "elevation", ylab = "")
for (i in 1:100){
  mod_rich <- aggregate(df_pasture[,3+i], list(df_pasture$point), sum)[,2]
  pasture_df_i <- data.frame(elevation = elevs_pasture$x, mod_rich = mod_rich)
  lw_np <- loess(mod_rich ~ elevation, data=pasture_df_i)
  
  mod_rich2 <- pt_richness[pt_richness$pasture == 1, i + 2]
  pasture_df_i2 <- data.frame(elevation = pt_richness$elevs[pt_richness$pasture == 1], mod_rich = mod_rich2)
  lw_np2 <- loess(mod_rich ~ elevation, data=pasture_df_i2)
  
  lines(seq(1250, 2550, 10), predict(lw_np, seq(1250, 2550, 10)) - predict(lw_np2, seq(1250, 2550, 10)), type = "l")
  
}
lines(seq(1250, 2550, 10), predict(lw_p_bmsom, seq(1250, 2550, 10)) - predict(lw_p_da, seq(1250, 2550, 10)), type = "l", lwd = 3, col = "red")



dev.off()


# forest bMSOM rep
plot(seq(1300, 2650, 10), predict(lw_f_bmsom, seq(1300, 2650, 10)), 
     type = "l", lwd = 3, col = "red", ylim = c(70, 130),
     xlab = "elevation", ylab = "richness difference")
for (i in 1:100){
  mod_rich <- aggregate(df_forest[,3+i], list(df_forest$point), sum)[,2]
  forest_df_i <- data.frame(elevation = elevs_forest$x, mod_rich = mod_rich)
  lw_nf <- loess(mod_rich ~ elevation, data=forest_df_i)
  
  lines(seq(1300, 2650, 10), predict(lw_nf, seq(1300, 2650, 10)), type = "l")
  
}
lines(seq(1300, 2650, 10), predict(lw_f_bmsom, seq(1300, 2650, 10)), type = "l", lwd = 3, col = "red")



# forest da rep
plot(seq(1300, 2650, 10), predict(lw_f_da, seq(1300, 2650, 10)), 
     type = "l", lwd = 3, col = "red", ylim = c(70, 130),
     xlab = "elevation", ylab = "richness difference")
for (i in 1:100){
  mod_rich2 <- pt_richness[pt_richness$pasture == 0, i + 2]
  forest_df_i2 <- data.frame(elevation = pt_richness$elevs[pt_richness$pasture == 0], mod_rich = mod_rich2)
  lw_nf2 <- loess(mod_rich ~ elevation, data=forest_df_i2)
  
  lines(seq(1300, 2650, 10), predict(lw_nf2, seq(1300, 2650, 10)), type = "l")
  
}
lines(seq(1300, 2650, 10), predict(lw_f_da, seq(1300, 2650, 10)), type = "l", lwd = 3, col = "red")









# pasture bMSOM rep
plot(seq(1300, 2650, 10), predict(lw_p_bmsom, seq(1300, 2650, 10)), 
     type = "l", lwd = 3, col = "red", ylim = c(30, 120),
     xlab = "elevation", ylab = "richness difference")
for (i in 1:100){
  mod_rich <- aggregate(df_pasture[,3+i], list(df_pasture$point), sum)[,2]
  pasture_df_i <- data.frame(elevation = elevs_pasture$x, mod_rich = mod_rich)
  lw_np <- loess(mod_rich ~ elevation, data=pasture_df_i)
  
  lines(seq(1300, 2650, 10), predict(lw_np, seq(1300, 2650, 10)), type = "l")
  
}
lines(seq(1300, 2650, 10), predict(lw_p_bmsom, seq(1300, 2650, 10)), type = "l", lwd = 3, col = "red")



# pasture da rep
plot(seq(1300, 2650, 10), predict(lw_p_da, seq(1300, 2650, 10)), 
     type = "l", lwd = 3, col = "red", ylim = c(30, 120),
     xlab = "elevation", ylab = "richness difference")
for (i in 1:100){
  mod_rich2 <- pt_richness[pt_richness$pasture == 1, i + 2]
  pasture_df_i2 <- data.frame(elevation = pt_richness$elevs[pt_richness$pasture == 1], mod_rich = mod_rich2)
  lw_np2 <- loess(mod_rich ~ elevation, data=pasture_df_i2)
  
  lines(seq(1300, 2650, 10), predict(lw_np2, seq(1300, 2650, 10)), type = "l")
  
}
lines(seq(1300, 2650, 10), predict(lw_p_da, seq(1300, 2650, 10)), type = "l", lwd = 3, col = "red")
