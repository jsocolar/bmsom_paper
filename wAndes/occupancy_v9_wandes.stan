// This is a Stan model for the full Colombia bird dataset, version 9.0, which is built on 7.0
// Changes:   switching to slice the occupancy intercept
//            better naming of data containers

// function to form a matrix with the same dimensions as ind, whose elements i,j are given by cov_u[ind[i,j]].
// This strategy for vectorizing the operation is due to Juho Timonen in Stan Slack post on 3 March 2021.
functions{
  matrix rt_mat(
    int r, // number of rows
    int c, // number of columns
    int[,] ind, // indices
    vector cov_u // unique covariate values
  ){
    int a_flat[r*c] = to_array_1d(ind);
    matrix[r,c] out = to_matrix(cov_u[a_flat], r, c, 0);
    return(out);
  }
  
  
  real partial_sum(
    // Function arguments:
      // Data slicing and indexing
        // Data slice    
          real[] intercept_slice, // slice of longest parameter vector
        // cutpoints for slicing                       
          int start,             // the starting row of the data slice
          int end,              // the ending row of the data slice  
      // main data  
        int[,] integer_data,      // data array
      // numbers of random effect levels
        int n_spCl,             // number of species-clusters
        int n_spSr,             // number of species-subregions
        int n_sp,               // number of species
        int n_fam,              // number of families
        int n_spObs,            // number of species-observers
      // max number of visits
        int n_visit_max,
      // Additional data
        // Continuous covariates
          // Elevations
            vector relev,    
            vector relev2, 
            vector lowland_x_relev,
            vector lowland_x_relev2,
          // Biogeography and traits
            vector elevMedian,
            vector elevBreadth,
            vector elevMedian_x_forestPresent,
            vector elevMedian_x_forestSpecialist,
            vector mass,
            vector elevMedian_x_pasture,
            vector elevBreadth_x_pasture,
            vector elevMedian_x_forestPresent_x_pasture,
            vector elevMedian_x_forestSpecialist_x_pasture,
            vector mass_x_pasture,
          // Nuisance
            matrix time,
            matrix time_x_elev,
          // Distance to range
            vector distance_to_range,
      // Parameters             
        // Occupancy
          // Intercept is sliced 
          // Slopes
            // Elevation effects
              real mu_b1_relev, 
              vector b1_relev_sp,    // slopes by species (zero-centered)
              real mu_b1_relev2,
              vector b1_relev2_sp,   // slopes by species (zero-centered)
              real b1_lowland,          // modified logit-quadratic lowland spp: intercept
              real b1_x_lowland_relev,  //    linear term
              real b1_x_lowland_relev2, //    quadratic term
          // Pasture effects
            real mu_b2_pasture,
            // Taxonomic effects (zero-centered)
            vector b2_pasture_sp,  // slopes by species
            vector b2_pasture_fam, // slopes by family 
          // Trait effects; covariates vary by species and all slopes are universal
            // Biogeography
              // Barriers
                real b3_mountain_barrier,
                real b3_valley_barrier,
              // Elevations
                real b3_elevMedian,
                real b3_elevBreadth,
              // Habitats
                real b3_forestPresent,
                real b3_forestSpecialist,
                real b3_tfSpecialist,
                real b3_dryForestPresent,
                real b3_floodDrySpecialist,
                real b3_aridPresent,
              // Migratory
                real b3_migratory,
              // Interactions
                real b3_x_elevMedian_forestPresent,
                real b3_x_elevMedian_forestSpecialist,
            // Functional traits
              real b3_mass,
              // diet (effects-coded; reference category = omnivore)
                real b3_dietInvert,
                real b3_dietCarn,
                real b3_dietFruitNect,
                real b3_dietGran,
          // pasture-x-trait interactions (same traits as b3).
            // Biogeography
              // Barriers
                real b4_mountain_barrier, 
                real b4_valley_barrier,
              // Elevations
                real b4_elevMedian,
                real b4_elevBreadth,
              // Habitats
                real b4_forestPresent,
                real b4_forestSpecialist,
                real b4_tfSpecialist,
                real b4_dryForestPresent,
                real b4_floodDrySpecialist,
                real b4_aridPresent,
              // Migratory
                real b4_migratory,
              // Interactions
                real b4_x_elevMedian_forestPresent,
                real b4_x_elevMedian_forestSpecialist,
            // Functional traits
              real b4_mass,
              // diet (effects coded; reference category = omnivore)
                real b4_dietInvert,
                real b4_dietCarn,
                real b4_dietFruitNect,
                real b4_dietGran,
          // Range maps
            real mu_b5_distance_to_range,
            vector b5_distance_to_range_sp, // species effects (zero-centered)
        // Detection
          // Intercept
            real mu_d0,
            // Taxonomic effects (zero-centered)
              vector d0_sp,          // species
              vector d0_fam,         // family
            // Species-observer effect (zero-centered)
              row_vector d0_spObs,
          // Slopes
            // Pasture effects
              real mu_d1_pasture,
              // Taxonomic effects (zero-centered)
                vector d1_pasture_sp,  // species
                vector d1_pasture_fam, // family
            // Nuisance effects
              // Time
                real mu_d2_time,
                vector d2_time_sp,     // species effects (zero-centered)
              // Observer
                real d2_obsDE,         // fixed slope for observer Edwards
              // Trait effects (universal slopes)
                real d3_mass,
                real d3_elevMedian,
                real d3_migratory,
                real d3_dietCarn,
                real d3_x_time_elevMedian
  ){   // End function arguments, begin computation
    // indexing variables
      int len = 1 + end - start;
    // data slice
      int data_slice[len,54] = integer_data[start:end,];
    // generate binary contrasts vector
      vector[2] binary_contrasts = [-1,1]';

    // variables for computation
      vector[len] lp;             // container for row-wise log-probability increments
      vector[len] logit_psi;      // rowise logit_psi
      row_vector[n_visit_max] logit_theta;  // single-row logit_theta
            // For points with < n_visit_max (4) visits, trailing elements of 
              // logit_theta are nonsense, computed from trailing 0s on the visit
              //  covariates. det_data has trailing -1s, ensuring that nonexistent 
              //  visits cannot accidentally slip into analysis.
      vector[len] logit_theta_vector; // vectorizable part of logit_theta
      matrix[len, n_visit_max] logit_theta_matrix;  // visit-specific logit_theta terms
      matrix[len, n_visit_max] obsDE = rt_mat(len, n_visit_max, data_slice[,51:54], d2_obsDE * binary_contrasts);

    // Computation:
      logit_psi =
        // Intercept
          to_vector(intercept_slice) +
        // Slopes
          // Elevation
            // Linear
              (mu_b1_relev + b1_relev_sp[data_slice[,9]]) .* relev[start:end] +   
            // Quadratic
              (mu_b1_relev2 + b1_relev2_sp[data_slice[,9]]) .* relev2[start:end] +
            // Lowland
              (b1_lowland*binary_contrasts)[data_slice[,15]] + b1_x_lowland_relev*lowland_x_relev[start:end] +
              b1_x_lowland_relev2*lowland_x_relev2[start:end] +
          // Pasture
            (mu_b2_pasture + b2_pasture_sp[data_slice[,9]] + b2_pasture_fam[data_slice[,10]]) .* binary_contrasts[data_slice[,16]] + // The random coefficient terms are slower than they probably could be.  Right now we multiply n_tot coefficient values by n_tot covariate values, but there's a LOT of redundancy here.
          // Biogeography
            // Barriers
              (b3_mountain_barrier * binary_contrasts)[data_slice[,17]] + (b3_valley_barrier * binary_contrasts)[data_slice[,18]] +
            // Elevations
              b3_elevMedian*elevMedian[start:end] + b3_elevBreadth*elevBreadth[start:end] +
            // Habitats
              (b3_forestPresent * binary_contrasts)[data_slice[,19]] + (b3_forestSpecialist * binary_contrasts)[data_slice[,20]] +
              (b3_tfSpecialist * binary_contrasts)[data_slice[,21]] + (b3_dryForestPresent * binary_contrasts)[data_slice[,22]] +
              (b3_floodDrySpecialist * binary_contrasts)[data_slice[,23]] + (b3_aridPresent * binary_contrasts)[data_slice[,24]] +
            // Migratory
              (b3_migratory * binary_contrasts)[data_slice[,25]] +
            // Interactions
              b3_x_elevMedian_forestPresent*elevMedian_x_forestPresent[start:end] +
              b3_x_elevMedian_forestSpecialist*elevMedian_x_forestSpecialist[start:end] +
          // Traits
            b3_mass*mass[start:end] +
            // Diet
              (b3_dietInvert * binary_contrasts)[data_slice[,26]] + (b3_dietCarn * binary_contrasts)[data_slice[,27]] +
              (b3_dietFruitNect * binary_contrasts)[data_slice[,28]] + (b3_dietGran * binary_contrasts)[data_slice[,29]] +
          // Pasture interactions
            // Biogeography
              // Barriers
                (b4_mountain_barrier * binary_contrasts)[data_slice[,30]] + (b4_valley_barrier * binary_contrasts)[data_slice[,31]] +
              // Elevations
                b4_elevMedian*elevMedian_x_pasture[start:end] + b4_elevBreadth*elevBreadth_x_pasture[start:end] +
              // Habitats
                (b4_forestPresent * binary_contrasts)[data_slice[,32]] + (b4_forestSpecialist * binary_contrasts)[data_slice[,33]] +
                (b4_tfSpecialist * binary_contrasts)[data_slice[,34]] + (b4_dryForestPresent * binary_contrasts)[data_slice[,35]] +
                (b4_floodDrySpecialist * binary_contrasts)[data_slice[,36]] + (b4_aridPresent * binary_contrasts)[data_slice[,37]] +
          // Migratory
            (b4_migratory * binary_contrasts)[data_slice[,38]] +
              // Interactions
                b4_x_elevMedian_forestPresent*elevMedian_x_forestPresent_x_pasture[start:end] +
                b4_x_elevMedian_forestSpecialist*elevMedian_x_forestSpecialist_x_pasture[start:end] +
              // Traits
                b4_mass*mass_x_pasture[start:end] +
                // Diet
                  (b4_dietInvert * binary_contrasts)[data_slice[,39]] + (b4_dietCarn * binary_contrasts)[data_slice[,40]] +
                  (b4_dietFruitNect * binary_contrasts)[data_slice[,41]] + (b4_dietGran * binary_contrasts)[data_slice[,42]] +
          // Range maps
            (mu_b5_distance_to_range + b5_distance_to_range_sp[data_slice[,9]]) .* distance_to_range[start:end];

      logit_theta_vector = 
        // Intercepts
          mu_d0 +
          // Taxonomic
            d0_sp[data_slice[,9]] + d0_fam[data_slice[,10]] + 
          // Species-observer: see loop
        // Slopes
          // Pasture
            (mu_d1_pasture + d1_pasture_sp[data_slice[,9]] + d1_pasture_fam[data_slice[,10]]) .* binary_contrasts[data_slice[,16]] +
          // Nuisance
            // Time: see matrix part
            // Observer: see matrix part
            // Trait
              d3_mass*mass[start:end] + d3_elevMedian*elevMedian[start:end] + 
              (d3_migratory * binary_contrasts)[data_slice[,25]] + 
              (d3_dietCarn * binary_contrasts)[data_slice[,27]];
        
      logit_theta_matrix =
        // Time
          diag_pre_multiply(mu_d2_time + d2_time_sp[data_slice[,9]], time[start:end]) + 
        // Observer
           obsDE + 
        // Interaction
          d3_x_time_elevMedian*time_x_elev[start:end];
        
      for (r in 1:len) {  // loop over species-points
        logit_theta = logit_theta_vector[r] + logit_theta_matrix[r,1:data_slice[r,6]] + d0_spObs[data_slice[r,11:(10 + data_slice[r,6])]];
        // likelihood
        if (data_slice[r,5] == 1) {
          lp[r] = log_inv_logit(logit_psi[r]) +  // likelihood of occupancy
          bernoulli_logit_lpmf(data_slice[r, 1:data_slice[r,6]] | logit_theta);   // conditional likelihood of observed history
        } else {
          lp[r] = log_sum_exp(
            log_inv_logit(logit_psi[r]) + // likelihood of occupancy
            sum(log1m_inv_logit(logit_theta)), // conditional likelihood of all-zero detection history
            log1m_inv_logit(logit_psi[r])); // likelihood of non-occupancy
        }
      } 
    return sum(lp);
  } // end of function
} // Close the functions block

data {
  // grainsize for reduce_sum 
    int<lower=1> grainsize;  
  // dimensions
    // random effects
      int<lower=1> n_spCl; // number of unique species:cluster 
      int<lower=1> n_spSr; // number of unique species:subregion
      int<lower=1> n_sp; // number of species
      int<lower=1> n_fam; // number of families
      int<lower=1> n_spObs; // number of unique species:observer
    // dataset size
      int<lower=1> n_tot; // number of total rows (number of unique species:point)
      int<lower=1> n_visit_max; // maximum number of visits to a point
  // main integer data
    int integer_data[n_tot, 54];
  // Continuous covariates
    vector[n_tot] distance_to_range;
    vector[n_tot] relev;
    vector[n_tot] relev2;
    vector[n_tot] lowland_x_relev;
    vector[n_tot] lowland_x_relev2;
    vector[n_tot] elevMedian;
    vector[n_tot] elevBreadth;
    vector[n_tot] mass;
    vector[n_tot] elevMedian_x_forestPresent;
    vector[n_tot] elevMedian_x_forestSpecialist;
    vector[n_tot] elevMedian_x_pasture;
    vector[n_tot] elevBreadth_x_pasture;
    vector[n_tot] mass_x_pasture;
    vector[n_tot] elevMedian_x_forestPresent_x_pasture;
    vector[n_tot] elevMedian_x_forestSpecialist_x_pasture;
    matrix[n_tot, n_visit_max] time;
    matrix[n_tot, n_visit_max] time_x_elev;
} // Close the data block

parameters {
  // Occupancy
    // Intercepts
      // Overall
        real mu_b0;
      // zero-centered taxonomic
        real<lower=0> sigma_b0_taxonomic;
        real<lower=0,upper=1> p_b0_taxonomic_sp;
        vector[n_sp] b0_sp_raw;
        vector[n_fam] b0_fam_raw;
      // zero-centered spatial
        real<lower=0> sigma_b0_spatial;
        real<lower=0,upper=1> p_b0_spatial_cl;
        vector[n_spCl] b0_spCl_raw;
        vector[n_spSr] b0_spSr_raw;
    // Slopes
      // zero-centered relev by species
        real mu_b1_relev;
        real<lower=0> sigma_b1_relev_sp;
        vector[n_sp] b1_relev_sp_raw;
      // zero-centered relev2 by species
        real mu_b1_relev2;
        real<lower=0> sigma_b1_relev2_sp;
        vector[n_sp] b1_relev2_sp_raw;
      // lowland interactions (universal)
        real b1_lowland;
        real b1_x_lowland_relev;
        real b1_x_lowland_relev2;
      // Pasture
        real mu_b2_pasture;
        // zero-centered taxonomic
          real<lower=0> sigma_b2_taxonomic;
          real<lower=0,upper=1> p_b2_taxonomic_sp;
          vector[n_sp] b2_pasture_sp_raw;
          vector[n_fam] b2_pasture_fam_raw;
      // Trait effects
        // Biogeography
          // Barriers
            real b3_mountain_barrier;
            real b3_valley_barrier;
          // Elevations
            real b3_elevMedian;
            real b3_elevBreadth;
          // Habitat
            real b3_forestPresent;
            real b3_forestSpecialist;
            real b3_tfSpecialist;
            real b3_dryForestPresent;
            real b3_floodDrySpecialist;
            real b3_aridPresent;
          // Migration
            real b3_migratory;
          // Interactions
            real b3_x_elevMedian_forestPresent;
            real b3_x_elevMedian_forestSpecialist;
        // Functional traits
          real b3_mass;
          real b3_dietInvert;
          real b3_dietCarn;
          real b3_dietFruitNect;
          real b3_dietGran;
      // Pasture-x-trait interactions
        // Biogeography
          // Barriers
            real b4_mountain_barrier; 
            real b4_valley_barrier;
          // Elevations
            real b4_elevMedian;
            real b4_elevBreadth;
          // Habitat
            real b4_forestPresent;
            real b4_forestSpecialist;
            real b4_tfSpecialist;
            real b4_dryForestPresent;
            real b4_floodDrySpecialist;
            real b4_aridPresent;
          // Migration
            real b4_migratory;
          // Interactions
            real b4_x_elevMedian_forestPresent;
            real b4_x_elevMedian_forestSpecialist;
        // Functional traits
          real b4_mass;
          real b4_dietInvert;
          real b4_dietCarn;
          real b4_dietFruitNect;
          real b4_dietGran;
      // Range maps
        real mu_b5_distance_to_range;
        real<lower=0> sigma_b5_distance_to_range_sp;
        vector[n_sp] b5_distance_to_range_sp_raw;
  // Detection
    // Intercepts
      real mu_d0;
      // zero-centered taxonomic
        real<lower=0> sigma_d0_taxonomic;
        real<lower=0,upper=1> p_d0_taxonomic_sp;
        vector[n_sp] d0_sp_raw;
        vector[n_fam] d0_fam_raw;
      // zero-centered species/observer
        real<lower=0> sigma_d0_spObs;
        row_vector[n_spObs] d0_spObs_raw;
    // Slopes
      // Pasture
        real mu_d1_pasture;
        // zero-centered taxonomic
          real<lower=0> sigma_d1_taxonomic;
          real<lower=0,upper=1> p_d1_taxonomic_sp;
          vector[n_sp] d1_pasture_sp_raw;
          vector[n_sp] d1_pasture_fam_raw;
      // zero-centered time by species
        real mu_d2_time;
        real<lower=0> sigma_d2_time_sp;
        vector[n_sp] d2_time_sp_raw;
      // observer
        real d2_obsDE;       
      // Trait effects
        real d3_mass;
        real d3_elevMedian;
        real d3_migratory;
        real d3_dietCarn;
        real d3_x_time_elevMedian;
} // End parameters block

transformed parameters {
  // Standard deviations for effects with combined variance terms
    real<lower=0> sigma_b0_spCl = sqrt(sigma_b0_spatial^2 * p_b0_spatial_cl);
    real<lower=0> sigma_b0_spSr = sqrt(sigma_b0_spatial^2 * (1-p_b0_spatial_cl));
    real<lower=0> sigma_b0_sp = sqrt(sigma_b0_taxonomic^2 * p_b0_taxonomic_sp);
    real<lower=0> sigma_b0_fam = sqrt(sigma_b0_taxonomic^2 * (1-p_b0_taxonomic_sp));
    real<lower=0> sigma_b2_pasture_sp = sqrt(sigma_b2_taxonomic^2 * p_b2_taxonomic_sp);
    real<lower=0> sigma_b2_pasture_fam = sqrt(sigma_b2_taxonomic^2 * (1-p_b2_taxonomic_sp));
    real<lower=0> sigma_d0_sp = sqrt(sigma_d0_taxonomic^2 * p_d0_taxonomic_sp);
    real<lower=0> sigma_d0_fam = sqrt(sigma_d0_taxonomic^2 * (1-p_d0_taxonomic_sp));
    real<lower=0> sigma_d1_pasture_sp = sqrt(sigma_d1_taxonomic^2 * p_d1_taxonomic_sp);
    real<lower=0> sigma_d1_pasture_fam = sqrt(sigma_d1_taxonomic^2 * (1-p_d1_taxonomic_sp));
}

model {
  // Manual non-centering to avoid sticky boundaries
    vector[n_sp] b0_sp = b0_sp_raw * sigma_b0_sp;
    vector[n_fam] b0_fam = b0_fam_raw * sigma_b0_fam;
    vector[n_spCl] b0_spCl = b0_spCl_raw * sigma_b0_spCl;
    vector[n_spSr] b0_spSr = b0_spSr_raw * sigma_b0_spSr;
    vector[n_sp] b1_relev_sp = b1_relev_sp_raw * sigma_b1_relev_sp;
    vector[n_sp] b1_relev2_sp = b1_relev2_sp_raw * sigma_b1_relev2_sp;
    vector[n_sp] b2_pasture_sp = b2_pasture_sp_raw * sigma_b2_pasture_sp;
    vector[n_fam] b2_pasture_fam = b2_pasture_fam_raw * sigma_b2_pasture_fam;
    vector[n_sp] b5_distance_to_range_sp = b5_distance_to_range_sp_raw * sigma_b5_distance_to_range_sp;
    vector[n_sp] d0_sp = d0_sp_raw * sigma_d0_sp;
    vector[n_fam] d0_fam = d0_fam_raw * sigma_d0_fam;
    row_vector[n_spObs] d0_spObs = d0_spObs_raw * sigma_d0_spObs;
    vector[n_sp] d1_pasture_sp = d1_pasture_sp_raw * sigma_d1_pasture_sp;
    vector[n_fam] d1_pasture_fam = d1_pasture_fam_raw * sigma_d1_pasture_fam;
    vector[n_sp] d2_time_sp = d2_time_sp_raw * sigma_d2_time_sp;
    
    real intercept_expanded[n_tot] = to_array_1d(mu_b0 + 
                                              b0_spCl[integer_data[,7]] + b0_spSr[integer_data[,8]] +
                                              b0_sp[integer_data[,9]] + b0_fam[integer_data[,10]]);

  // Priors and Jacobian adjustments
    // Occupancy
      // Intercepts
        mu_b0 ~ normal(-7, 2.5);
        // Spatial
          sigma_b0_spatial ~ normal(0, 3);
          // implicit uniform prior on p_b0_spatial_cl
          b0_spCl_raw ~ std_normal();
          b0_spSr_raw ~ std_normal();
        // Taxonomic
          sigma_b0_taxonomic ~ normal(0, 2);
          // implicit uniform prior on p_b0_taxonomic_sp
          b0_sp_raw ~ std_normal();
          b0_fam_raw ~ std_normal();
      // Slopes
        // Elevations
          // Linear
            mu_b1_relev ~ normal(0, 5);
            sigma_b1_relev_sp ~ normal(0, 2);
            b1_relev_sp_raw ~ std_normal();
          // Quadratic
            mu_b1_relev2 ~ normal(0, 5);
            sigma_b1_relev2_sp ~ normal(0, 2);
            b1_relev2_sp_raw ~ std_normal();
          // Lowland
            b1_lowland ~ normal(0, 1);
            b1_x_lowland_relev ~ normal(0, 5);
            b1_x_lowland_relev2 ~ normal(0, 5);
        // Pasture
          mu_b2_pasture ~ normal(0, 1);
          // Taxonomic
            sigma_b2_taxonomic ~ normal(0, 1);
            // implicit uniform prior on p_b2_taxonomic_sp
            b2_pasture_sp_raw ~ std_normal();
            b2_pasture_fam_raw ~ std_normal();
        // Biogeographic
          b3_mountain_barrier ~ normal(0, 1);
          b3_valley_barrier ~ normal(0, 1);
          b3_elevMedian ~ normal(0, 1);
          b3_elevBreadth ~ normal(0, 1);
          b3_forestPresent ~ normal(0, 1);
          b3_forestSpecialist ~ normal(0, 1);
          b3_tfSpecialist ~ normal(0, 1);
          b3_dryForestPresent ~ normal(0, 1);
          b3_floodDrySpecialist ~ normal(0, 1);
          b3_aridPresent ~ normal(0, 1);
          b3_migratory ~ normal(0, 1);
          b3_x_elevMedian_forestPresent ~ normal(0, 1);
          b3_x_elevMedian_forestSpecialist ~ normal(0, 1);
        // Functional
          b3_mass ~ normal(0, .5);
          b3_dietInvert ~ normal(0, .5);
          b3_dietCarn ~ normal(0, .5);
          b3_dietFruitNect ~ normal(0, .5);
          b3_dietGran ~ normal(0, .5);
        // Pasture interactions
          // Biogeographic
            b4_mountain_barrier ~ normal(0, 1);
            b4_valley_barrier ~ normal(0, 1);
            b4_elevMedian ~ normal(0, 1);
            b4_elevBreadth ~ normal(0, 1);
            b4_forestPresent ~ normal(0, 1);
            b4_forestSpecialist ~ normal(0, 1);
            b4_tfSpecialist ~ normal(0, 1);
            b4_dryForestPresent ~ normal(0, 1);
            b4_floodDrySpecialist ~ normal(0, 1);
            b4_aridPresent ~ normal(0, 1);
            b4_migratory ~ normal(0, 1);
            b4_x_elevMedian_forestPresent ~ normal(0, 1);
            b4_x_elevMedian_forestSpecialist ~ normal(0, 1);
          // Functional
            b4_mass ~ normal(0, .5);
            b4_dietInvert ~ normal(0, .5);
            b4_dietCarn ~ normal(0, .5);
            b4_dietFruitNect ~ normal(0, .5);
            b4_dietGran ~ normal(0, .5);
        // Distance-to-range
          mu_b5_distance_to_range ~ normal(0, 5);
          sigma_b5_distance_to_range_sp ~ normal(0, 2);
          b5_distance_to_range_sp_raw ~ std_normal();
    // Detection
      // Intercepts
        mu_d0 ~ normal(-3,1);
        // Taxonomic
          sigma_d0_taxonomic ~ normal(0, 2);
          // Implicit uniform prior on p_d0_taxonomic
          d0_sp_raw ~ std_normal();
          d0_fam_raw ~ std_normal();
        // Observer
          sigma_d0_spObs ~ normal(0, 1);
          d0_spObs_raw ~ std_normal();
      // Slopes
        // Pasture
          mu_d1_pasture ~ normal(0,.75);
          // Taxonomic
            sigma_d1_taxonomic ~ normal(0,1);
            // Implicit uniform prior on p_d1_taxonomic
            d1_pasture_sp_raw ~ std_normal();
            d1_pasture_fam_raw ~ std_normal();
        // Time
          mu_d2_time ~ normal(0, .5);
          sigma_d2_time_sp ~ normal(0, 1);
          d2_time_sp_raw ~ std_normal();
        // Observer
          d2_obsDE ~ normal(0, .25);
        // Traits
          d3_mass ~ normal(0, .5);
          d3_elevMedian ~ normal(0, .5);
          d3_migratory ~ normal(0, 1);
          d3_dietCarn ~ normal(0, .5);
          d3_x_time_elevMedian ~ normal(0,.5);
  // Likelihood computed via reduce_sum
    target += reduce_sum(
      // partial_sum function
        partial_sum,
      // intercept
        intercept_expanded,
      // grainsize
        grainsize,
      // main integer data
        integer_data,
      // variable sizes
        n_spCl, n_spSr, n_sp, n_fam, n_spObs,
      // max number of visits
        n_visit_max,
      // continuous covariates
          // Elevations
            relev,
            relev2,
            lowland_x_relev,
            lowland_x_relev2,
          // Biogeography and traits
            elevMedian,
            elevBreadth,
            elevMedian_x_forestPresent,
            elevMedian_x_forestSpecialist,
            mass,
            elevMedian_x_pasture,
            elevBreadth_x_pasture,
            elevMedian_x_forestPresent_x_pasture,
            elevMedian_x_forestSpecialist_x_pasture,
            mass_x_pasture,
          // Nuisance
            time,
            time_x_elev,
          // Distance to range
            distance_to_range,
      // parameters:
        // Occupancy elevation
          mu_b1_relev, b1_relev_sp, mu_b1_relev2, b1_relev2_sp,
          b1_lowland, b1_x_lowland_relev, b1_x_lowland_relev2,
        // Occupancy pasture
          mu_b2_pasture, b2_pasture_sp, b2_pasture_fam,
        // Occupancy traits
          b3_mountain_barrier, b3_valley_barrier,
          b3_elevMedian, b3_elevBreadth,
          b3_forestPresent, b3_forestSpecialist, b3_tfSpecialist,
          b3_dryForestPresent, b3_floodDrySpecialist,
          b3_aridPresent,
          b3_migratory,
          b3_x_elevMedian_forestPresent, b3_x_elevMedian_forestSpecialist,
          b3_mass,
          b3_dietInvert, b3_dietCarn, b3_dietFruitNect, b3_dietGran,
          b4_mountain_barrier, b4_valley_barrier,
          b4_elevMedian, b4_elevBreadth,
          b4_forestPresent, b4_forestSpecialist, b4_tfSpecialist, b4_dryForestPresent, b4_floodDrySpecialist, b4_aridPresent,
          b4_migratory,
          b4_x_elevMedian_forestPresent, b4_x_elevMedian_forestSpecialist,
          b4_mass,
          b4_dietInvert, b4_dietCarn, b4_dietFruitNect, b4_dietGran,
        // Occupancy distance-to-range
          mu_b5_distance_to_range, b5_distance_to_range_sp,
        // Detection intercepts
          mu_d0, d0_sp, d0_fam, d0_spObs,
        // Detection pasture
          mu_d1_pasture, d1_pasture_sp, d1_pasture_fam,
        // Detection nuisance
          mu_d2_time, d2_time_sp, d2_obsDE,
        // Detection traits
          d3_mass, d3_elevMedian, d3_migratory, d3_dietCarn, d3_x_time_elevMedian
    ); // end reduce_sum call
} // end model block

