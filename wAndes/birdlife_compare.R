library(cmdstanr)
library(posterior)

# A list of traits as given by the BirdLife Datazone for each species of bird.
load("/Users/JacobSocolar/Dropbox/Work/Useful_data/BirdlifeTraits/birdlife_traits.Rdata")

wandes_samples <- read_cmdstan_csv(list.files("/Users/jacobsocolar/Dropbox/Work/Occupancy/biogeographicMSOM/wandes_stan_outputs",
                                              full.names = T))
samples <- as_draws_df(wandes_samples$post_warmup_draws)

wandes_data <- readRDS("/Users/JacobSocolar/Dropbox/Work/Occupancy/biogeographicMSOM/wandes_data/wandes_bsd9_package.RDS")

length(unique(wandes_data$data$integer_data$id_sp[wandes_data$data$integer_data$Q == 1]))
length(unique(wandes_data$data$integer_data$id_sp))

birds <- readRDS("/Users/jacobSocolar/Dropbox/Work/Colombia/Data/Analysis/birds.RDS")
wandes_birds <- birds[grepl("wandes", birds$subregion),]

w_sp_obs_columns <- c(wandes_birds$sp_obs1, wandes_birds$sp_obs2, wandes_birds$sp_obs3, wandes_birds$sp_obs4)
w_sp_obs_columns[grep("__NA", w_sp_obs_columns)] <- NA
w_sp_obs_matrix <- matrix(as.integer(as.factor(w_sp_obs_columns)), ncol=4)
w_sp_obs_matrix[is.na(w_sp_obs_matrix)] <- 0
wandes_birds$sp_obs_matrix <- w_sp_obs_matrix
traitdata <- wandes_birds[!duplicated(wandes_birds$species),]
integerdata <- wandes_data$data$integer_data[!duplicated(wandes_data$data$integer_data$id_sp), ]


forest_int <- pasture_offset <- matrix(data = NA, nrow = 910, ncol = 6000)

ec <- c(-1,1)
for(sp in 1:910){
  print(sp)
  forest_int[sp,] <- t(as.vector(samples$mu_b0 +
                                    samples[,paste0('b0_sp_raw[',sp,']')] * samples$sigma_b0_sp + samples[,paste0('b0_fam_raw[',integerdata$id_fam[sp],']')] * samples$sigma_b0_fam +
                                    
                                    samples[,'b1_lowland']*ec[integerdata$lowland[sp]] +
                                    
                                    samples[,'b3_mountain_barrier']*ec[integerdata$mountain_barrier[sp]] + samples[,'b3_valley_barrier']*ec[integerdata$valley_barrier[sp]] + 
                                    
                                    samples[,'b3_elevMedian']*traitdata$elev_median_scaled[sp] + samples[,'b3_elevBreadth']*traitdata$elev_breadth_scaled[sp] +
                                    
                                    samples[,'b3_forestPresent']*ec[integerdata$forestPresent[sp]] + samples[,'b3_forestSpecialist']*ec[integerdata$forestSpecialist[sp]] + 
                                    samples[,'b3_tfSpecialist']*ec[integerdata$tfSpecialist[sp]] + samples[,'b3_dryForestPresent']*ec[integerdata$dryForestPresent[sp]] +
                                    samples[,'b3_floodDrySpecialist']*ec[integerdata$floodDrySpecialist[sp]] + 
                                    samples[,'b3_aridPresent']*ec[integerdata$aridPresent[sp]] +
                                    
                                    samples[,'b3_migratory']*ec[integerdata$migratory[sp]] + samples[,'b3_mass']*traitdata$log_mass_scaled[sp] +
                                    
                                    samples[,'b3_dietInvert']*ec[integerdata$dietInvert[sp]] + samples[,'b3_dietCarn']*ec[integerdata$dietCarn[sp]] + 
                                    samples[,'b3_dietFruitNect']*ec[integerdata$dietFruitNect[sp]] + samples[,'b3_dietGran']*ec[integerdata$dietGran[sp]] +
                                    
                                    samples[,'b3_x_elevMedian_forestPresent']*traitdata$elev_median_scaled[sp]*ec[integerdata$forestPresent[sp]] +
                                    samples[,'b3_x_elevMedian_forestSpecialist']*traitdata$elev_median_scaled[sp]*ec[integerdata$forestSpecialist[sp]]))
  
  
  pasture_offset[sp,] <- t(as.vector(samples$mu_d0 +
                                        samples[,paste0('b2_pasture_sp_raw[',sp,']')] * samples$sigma_b2_pasture_sp + samples[,paste0('b2_pasture_fam_raw[',integerdata$id_fam[sp],']')] * samples$sigma_b2_pasture_fam +
                                        samples[,'b4_mountain_barrier']*ec[integerdata$mountain_barrier[sp]] + samples[,'b4_valley_barrier']*ec[integerdata$valley_barrier[sp]] + 
                                        
                                        samples[,'b4_elevMedian']*traitdata$elev_median_scaled[sp] + samples[,'b4_elevBreadth']*traitdata$elev_breadth_scaled[sp] +
                                        
                                        samples[,'b4_forestPresent']*ec[integerdata$forestPresent[sp]] + samples[,'b4_forestSpecialist']*ec[integerdata$forestSpecialist[sp]] + 
                                        samples[,'b4_tfSpecialist']*ec[integerdata$tfSpecialist[sp]] + samples[,'b4_dryForestPresent']*ec[integerdata$dryForestPresent[sp]] +
                                        samples[,'b4_floodDrySpecialist']*ec[integerdata$floodDrySpecialist[sp]] +
                                        samples[,'b4_aridPresent']*ec[integerdata$aridPresent[sp]] +
                                        
                                        samples[,'b4_migratory']*ec[integerdata$migratory[sp]] + samples[,'b4_mass']*traitdata$log_mass_scaled[sp] +
                                        
                                        samples[,'b4_dietInvert']*ec[integerdata$dietInvert[sp]] + samples[,'b4_dietCarn']*ec[integerdata$dietCarn[sp]] + 
                                        samples[,'b4_dietFruitNect']*ec[integerdata$dietFruitNect[sp]] + samples[,'b4_dietGran']*ec[integerdata$dietGran[sp]] +
                                        
                                        samples[,'b4_x_elevMedian_forestPresent']*traitdata$elev_median_scaled[sp]*ec[integerdata$forestPresent[sp]]+
                                        samples[,'b4_x_elevMedian_forestSpecialist']*traitdata$elev_median_scaled[sp]*ec[integerdata$forestSpecialist[sp]]))
}


birdlife_compare <- data.frame(pasture = apply(pasture_offset, 1, median), birdlife_dep = NA)
for(i in 1:nrow(birdlife_compare)){
  tryCatch(birdlife_compare$birdlife_dep[i] <- birdlife_traits$forest_dep[birdlife_traits$names == gsub("_", " ", traitdata$species[i])], error = function(x){NA})
}
birdlife_compare$birdlife_high <- birdlife_compare$birdlife_dep == "High"
birdlife_compare$birdlife_low <- birdlife_compare$birdlife_dep %in% c("Low", "Does not normally occur in forest")
birdlife_compare$Q <- integerdata$Q

bco <- birdlife_compare[birdlife_compare$Q==1,]
bcn <- birdlife_compare[birdlife_compare$Q==0,]

boxplot(bco$pasture[bco$birdlife_high], bco$pasture[bco$birdlife_low],
        bcn$pasture[bcn$birdlife_high], bcn$pasture[bcn$birdlife_low],
        ylab = "pasture index",
        names = c("observed high", "observed low", "never-observed high", "never-observed low"))
