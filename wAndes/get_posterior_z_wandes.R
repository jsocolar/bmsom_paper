# Functions to extract a posterior iteration for several quantities of interest

# get_psi_components extracts the components of the linear predictor for occupancy that do not vary by point
# (i.e. one forest value and one pasture offset per species)
# get_theta_component does the same for the linear predictor for detection
# get_z_components calls get_psi_components and get_theta_components and returns all results in a list.
# get_prediction_components calls get_psi_components and get_theta_component, and additionally returns the coefficients for terms
# that vary by point (elevation, distance-to-range) and the hyperparameters for random effects that vary
# by point.
# get_Z_probs gives a posterior iteration for psi and Z at the locations of the actual sampling points
# get_data_rep simulates a posterior dataset

##### get_Z_probs #####
get_Z_probs <- function(draws, iter, z_info, spatial_effect = "include"){
  # Get the posterior occupancy probabilities for each species-point.
  
  # ARGUMENTS
  # draws: a draws_df from occupancy_v6 or higher
  # iter: a posterior iteration
  # z_info: the $data part of a data object used in fitting, with points and 
  #         species names appended from the birds dataframe. This should be 
  #         from bird_stan_data6_package; the data are stable.
  # spatial_effect: "include":  include fitted spatial effects
  #                 "exclude":  exclude fitted spatial effects
  #                 "resample": resample spatial effects from hyperparameters
  
  # OUTPUT
  # output: dataframe with columns for the point, the species, the observed Q, the fitted psi, 
  #           the fitted conditional probability of at least one detection, the posterior
  #           probability that Z == 1, and the conditional detection probabilities of each visit
  
  # note: internally called helper functions defined below
  
  zc <- get_prediction_components(draws, iter, z_info)
  if(!all.equal(zc$id_sp, 1:max(zc$id_sp))){stop("zc is mis-ordered; fix by ordering but make sure you understand why this happened; it's not supposed to happen")}
  
  logit_psi_no_spatial <- 
    # forest values if we're in forest; pasture values if we're in pasture
    zc$logit_psi_0[z_info$id_sp] + zc$logit_psi_pasture_offset[z_info$id_sp]*z_info$pasture +
    # elevation term
    zc$b1_relev_sp[z_info$id_sp] * z_info$relev + zc$b1_relev2_sp[z_info$id_sp] * z_info$relev2 +
    (zc$b1_x_lowland_relev * zc$lowland)[z_info$id_sp] * z_info$relev + 
    (zc$b1_x_lowland_relev2 * zc$lowland)[z_info$id_sp] * z_info$relev2 +
    # distance_to_range term
    zc$b5_distance_to_range_sp[z_info$id_sp] * z_info$distance_to_range
  
  if(spatial_effect == "include"){
    logit_psi <- 
      logit_psi_no_spatial +
      as.numeric(draws[iter, "sigma_b0_spCl"]) * as.numeric(draws[iter, paste0("b0_spCl_raw[", 
                                                                               format(z_info$id_spCl, scientific = F, trim = "true", justify = "none"),
                                                                               "]")]) +  # the format() statement above is necessary because as.character() returns "1e+05" and "2e+05" for those two indices
      as.numeric(draws[iter, "sigma_b0_spSr"]) * as.numeric(draws[iter, paste0("b0_spSr_raw[", 
                                                                               as.character(z_info$id_spSr),"]")])
    
  }else if(spatial_effect == "exclude"){
    logit_psi <- logit_psi_no_spatial
  }else if(spatial_effect == "resample"){
    cluster_effects <- rnorm(max(z_info$id_spCl), 0, zc$sigma_sp_cl)
    subregion_effects <- rnorm(max(z_info$id_spSr), 0, zc$sigma_sp_sr)
    logit_psi <- logit_psi_no_spatial + cluster_effects[z_info$id_spCl] + subregion_effects[z_info$id_spSr]
  }
  
  psi <- boot::inv.logit(logit_psi)
  
  logit_theta_matrix <- replicate(4, zc$logit_theta_0[z_info$id_sp] + zc$logit_theta_pasture_offset[z_info$id_sp]*z_info$pasture) +
    sweep(as.matrix(z_info[,c("time.1", "time.2", "time.3", "time.4")]), MARGIN = 1, zc$d2_time_sp[z_info$id_sp], `*`) +
    zc$d2_obsDE[1] * as.matrix(z_info[,c("obsDE.1", "obsDE.2", "obsDE.3", "obsDE.4")]) +
    sweep(as.matrix(z_info[,c("time.1", "time.2", "time.3", "time.4")]), MARGIN = 1, zc$d3_x_time_elevMedian[1] * z_info$elevMedian, `*`) +
    as.numeric(draws[iter, "sigma_d0_spObs"])*
    cbind(as.numeric(draws[iter, paste0("d0_spObs_raw[", as.character(z_info$id_spObs.1), "]")]),
          as.numeric(draws[iter, paste0("d0_spObs_raw[", as.character(z_info$id_spObs.1), "]")]),
          as.numeric(draws[iter, paste0("d0_spObs_raw[", as.character(z_info$id_spObs.1), "]")]),
          as.numeric(draws[iter, paste0("d0_spObs_raw[", as.character(z_info$id_spObs.1), "]")])
    )
  
  theta_matrix <- boot::inv.logit(logit_theta_matrix)
  
  # get the likelihood of zero detections given Z == 1
  zero_det_lik <- apply(cbind(z_info$nv, theta_matrix), 1, FUN = function(x){prod(1 - x[2:(1+x[1])])})
  
  # get the relative probabilities of Z=1 and Z=0 give Q=0
  z1_lik <- psi*zero_det_lik
  z0_lik <- 1-psi
  # Normalize to get the posterior probability that Z=1
  z1_prob <- z1_lik/(z1_lik + z0_lik)
  
  df_out <- data.frame(point = z_info$point, species = z_info$species, Q = z_info$Q, psi = psi)
  df_out$pdet <- 1 - zero_det_lik
  df_out$Z_prob <- df_out$Q + (1 - df_out$Q)*z1_prob
  df_out$theta_matrix <- theta_matrix
  
  return(df_out)
}

##### get_data_rep #####
get_data_rep <- function(z_prob_df){
  # Function simulate a posterior data replicate
  
  # ARGUMENT: output of get_Z_probs()
  
  # OUTPUT: 3-element list
  #         element "post" is a simulated posterior replicate of whether each point has at least one detection, conditioning on the observed Q
  #         element "mixed" is a simulated replicate of whether each point has at least one detection, without conditioning on the observed Q
  #         element "histories" is a simulated detection history at each point, conditioning on the observed Q
  
  Z_sample <- rbinom(nrow(z_prob_df), 1, z_prob_df$Z_prob)
  D_sample <- rbinom(nrow(z_prob_df), 1, z_prob_df$pdet)
  data_rep_post <- Z_sample*D_sample
  data_rep_mixed <- rbinom(nrow(z_prob_df), 1, z_prob_df$psi*z_prob_df$pdet)
  cond_detmatrix <- matrix(rbinom(length(z_prob_df$theta_matrix), 1, z_prob_df$theta_matrix), ncol = 4)
  detmatrix <- sweep(cond_detmatrix, MARGIN = 1, Z_sample, `*`)
  histories <- apply(detmatrix, 1, function(x){paste(x, collapse = "")})
  return(list(post=data_rep_post, mixed=data_rep_mixed, histories=histories))
}


##### Helper Functions #####
get_psi_components <- function(draws, iter, z_info){
  # Function to get the species-specific components of logit_psi that do not vary across points with a habitat class
  
  # Get the unique species
  b <- z_info[!duplicated(z_info$id_sp),]
  
  # Get the expectation at "pasture = 0" (then we'll add the pasture offset if we're in pasture and subtract the 
  # pasture offset if we're in forest)
  logit_psi_0 <- 
    # Intercept
    as.numeric(draws[iter, "mu_b0"]) +
    as.numeric(draws[iter, paste0("b0_sp_raw[", b$id_sp, "]")])*as.numeric(draws[iter, "sigma_b0_sp"]) +
    as.numeric(draws[iter, paste0("b0_fam_raw[", b$id_fam, "]")])*as.numeric(draws[iter, "sigma_b0_fam"]) +
    # Constant part of elevation relationship
    as.numeric(draws[iter, "b1_lowland"])*b$lowland + 
    # Biogeography
    as.numeric(draws[iter, "b3_mountain_barrier"])*b$mountain_barrier + as.numeric(draws[iter, "b3_valley_barrier"])*b$valley_barrier +
    as.numeric(draws[iter, "b3_elevMedian"])*b$elevMedian + as.numeric(draws[iter, "b3_elevBreadth"])*b$elevBreadth +
    as.numeric(draws[iter, "b3_forestPresent"])*b$forestPresent + as.numeric(draws[iter, "b3_forestSpecialist"])*b$forestSpecialist +
    as.numeric(draws[iter, "b3_tfSpecialist"])*b$tfSpecialist + as.numeric(draws[iter, "b3_dryForestPresent"])*b$dryForestPresent +
    as.numeric(draws[iter, "b3_floodDrySpecialist"])*b$floodDrySpecialist + 
    as.numeric(draws[iter, "b3_aridPresent"])*b$aridPresent + as.numeric(draws[iter, "b3_migratory"])*b$migratory + 
    as.numeric(draws[iter, "b3_x_elevMedian_forestPresent"])*b$forestPresent*b$elevMedian +
    as.numeric(draws[iter, "b3_x_elevMedian_forestSpecialist"])*b$forestSpecialist*b$elevMedian +
    # Functional
    as.numeric(draws[iter, "b3_mass"])*b$mass + as.numeric(draws[iter, "b3_dietInvert"])*b$dietInvert + as.numeric(draws[iter, "b3_dietCarn"])*b$dietCarn +
    as.numeric(draws[iter, "b3_dietFruitNect"])*b$dietFruitNect + as.numeric(draws[iter, "b3_dietGran"])*b$dietGran 
  
  # Get the magnitude of the pasture offset
  logit_psi_pasture_offset <- 
    # Main effect
    as.numeric(draws[iter, "mu_b2_pasture"]) +
    as.numeric(draws[iter, paste0("b2_pasture_sp_raw[", b$id_sp, "]")])*as.numeric(draws[iter, "sigma_b2_pasture_sp"]) + 
    as.numeric(draws[iter, paste0("b2_pasture_fam_raw[", b$id_fam, "]")])*as.numeric(draws[iter, "sigma_b2_pasture_fam"]) +
    # Biogeographic interactions
    as.numeric(draws[iter, "b4_mountain_barrier"])*b$mountain_barrier + as.numeric(draws[iter, "b4_valley_barrier"])*b$valley_barrier +
    as.numeric(draws[iter, "b4_elevMedian"])*b$elevMedian + as.numeric(draws[iter, "b4_elevBreadth"])*b$elevBreadth +
    as.numeric(draws[iter, "b4_forestPresent"])*b$forestPresent + as.numeric(draws[iter, "b4_forestSpecialist"])*b$forestSpecialist +
    as.numeric(draws[iter, "b4_tfSpecialist"])*b$tfSpecialist + as.numeric(draws[iter, "b4_dryForestPresent"])*b$dryForestPresent +
    as.numeric(draws[iter, "b4_floodDrySpecialist"])*b$floodDrySpecialist + 
    as.numeric(draws[iter, "b4_aridPresent"])*b$aridPresent + as.numeric(draws[iter, "b4_migratory"])*b$migratory + 
    as.numeric(draws[iter, "b4_x_elevMedian_forestPresent"])*b$forestPresent*b$elevMedian +
    as.numeric(draws[iter, "b4_x_elevMedian_forestSpecialist"])*b$forestSpecialist*b$elevMedian +
    # Functional interactions
    as.numeric(draws[iter, "b4_mass"])*b$mass + as.numeric(draws[iter, "b4_dietInvert"])*b$dietInvert + as.numeric(draws[iter, "b4_dietCarn"])*b$dietCarn +
    as.numeric(draws[iter, "b4_dietFruitNect"])*b$dietFruitNect + as.numeric(draws[iter, "b4_dietGran"])*b$dietGran
  
  return(data.frame(id_sp = b$id_sp, logit_psi_0 = logit_psi_0, logit_psi_pasture_offset = logit_psi_pasture_offset))
}

get_theta_components <- function(draws, iter, z_info){
  # Function to get the species-specific components of logit_theta that do not vary across points or visits within a habitat class
  b <- z_info[!duplicated(z_info$id_sp),]
  
  logit_theta_0 <- 
    # Intercept
    as.numeric(draws[iter, "mu_d0"]) +
    as.numeric(draws[iter, paste0("d0_sp_raw[", b$id_sp, "]")])*as.numeric(draws[iter, "sigma_d0_sp"]) + 
    as.numeric(draws[iter, paste0("d0_fam_raw[", b$id_fam, "]")])*as.numeric(draws[iter, "sigma_d0_fam"]) +
    as.numeric(draws[iter, "d3_mass"])*b$mass + as.numeric(draws[iter, "d3_elevMedian"])*b$elevMedian +
    as.numeric(draws[iter, "d3_migratory"])*b$migratory + as.numeric(draws[iter, "d3_dietCarn"])*b$dietCarn
  
  logit_theta_pasture_offset <- 
    # Main effect
    as.numeric(draws[iter, "mu_d1_pasture"]) +
    as.numeric(draws[iter, paste0("d1_pasture_sp_raw[", b$id_sp, "]")])*as.numeric(draws[iter, "sigma_d1_pasture_sp"]) + 
    as.numeric(draws[iter, paste0("d1_pasture_fam_raw[", b$id_fam, "]")])*as.numeric(draws[iter, "sigma_d1_pasture_sp"])
  # No interaction terms
  
  return(data.frame(id_sp = b$id_sp, logit_theta_0 = logit_theta_0, logit_theta_pasture_offset = logit_theta_pasture_offset))
}

get_z_components <- function(draws, iter, z_info){
  # Function that calls get_psi_components and get_theta_components and packages up the results in a single data.frame
  
  b <- z_info[!duplicated(z_info$id_sp),]
  psi_c <- get_psi_components(draws, iter, z_info)
  theta_c <- get_theta_components(draws, iter, z_info)
  if(all.equal(psi_c$id_sp, theta_c$id_sp)){
    return(cbind(psi_c, theta_c[,c("logit_theta_0", "logit_theta_pasture_offset")]))
  }else{stop('psi and theta have different species ids!')}
}

get_prediction_components <- function(draws, iter, z_info){
  # Function that calls get_z_components and packages the result together with the coefficient values
  # (fixed effects) or standard deviations (random effects) for covariates that vary by points or visits 
  
  zc <- get_z_components(draws, iter, z_info)
  zc$lowland <- z_info$lowland[!duplicated(z_info$id_sp)]
  zc$b1_relev_sp <- as.numeric(draws[iter, "mu_b1_relev"]) +
    as.numeric(draws[iter, paste0("b1_relev_sp_raw[", zc$id_sp, "]")])*as.numeric(draws[iter, "sigma_b1_relev_sp"])
  zc$b1_relev2_sp <- as.numeric(draws[iter, "mu_b1_relev2"]) + 
    as.numeric(draws[iter, paste0("b1_relev2_sp_raw[", zc$id_sp, "]")])*as.numeric(draws[iter, "sigma_b1_relev2_sp"])
  zc$b1_x_lowland_relev <- as.numeric(draws[iter, "b1_x_lowland_relev"])
  zc$b1_x_lowland_relev2 <- as.numeric(draws[iter, "b1_x_lowland_relev2"])
  zc$b5_distance_to_range_sp <- 
    as.numeric(draws[iter, "mu_b5_distance_to_range"]) +
    as.numeric(draws[iter, paste0("b5_distance_to_range_sp_raw[", zc$id_sp, "]")])*as.numeric(draws[iter, "sigma_b5_distance_to_range_sp"])
  zc$sigma_sp_cl <- as.numeric(draws[iter, "sigma_b0_spCl"])
  zc$sigma_sp_sr <- as.numeric(draws[iter, "sigma_b0_spSr"])
  zc$d2_time_sp <-  
    as.numeric(draws[iter, "mu_d2_time"]) + 
    as.numeric(draws[iter, paste0("d2_time_sp_raw[", zc$id_sp, "]")])*as.numeric(draws[iter, "sigma_d2_time_sp"])
  zc$d2_obsDE <- as.numeric(draws[iter, "d2_obsDE"])
  zc$d3_x_time_elevMedian <- as.numeric(draws[iter, "d3_x_time_elevMedian"])
  zc$sigma_sp_obs <- as.numeric(draws[iter, "sigma_d0_spObs"])
  return(zc)
}
