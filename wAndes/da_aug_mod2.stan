// data-augmented Stan model for the West Andes analysis presented in 
// "Biogeographic multi-species occupancy models for large-scale survey data"
// Note: code generated with brms 2.16.1

functions {
 /* compute correlated group-level effects
  * Args: 
  *   z: matrix of unscaled group-level effects
  *   SD: vector of standard deviation parameters
  *   L: cholesky factor correlation matrix
  * Returns: 
  *   matrix of scaled group-level effects
  */ 
  matrix scale_r_cor(matrix z, vector SD, matrix L) {
    // r is stored in another dimension order than z
    return transpose(diag_pre_multiply(SD, L) * z);
  }
    real occupancy_augmented_lpmf(
    int[] y, // detection data
    vector mu, // lin pred for detection
    vector occ, // lin pred for occupancy. Only the first vint1[1] elements matter.
    real Omega, // lin pred for availability.  Only the first element matters
    int[] vint1, // # units (n_unit). Elements after 1 irrelevant.
    int[] vint2, // # sampling events per unit (n_rep). Elements after vint1[1] irrelevant.
    int[] vint3, // Indicator for > 0 detections (Q). Elements after vint1[1] irrelevant.
    
    int[] vint4, // # species (observed + augmented)
    int[] vint5, // # Indicator for species was observed.  Elements after vint4[1] irrelevant
    
    int[] vint6, // species
  
  // indices for jth repeated sampling event to each unit (elements after vint1[1] irrelevant):
    int[] vint7,
    int[] vint8,
    int[] vint9,
    int[] vint10
) {
  // Create array of the rep indices that correspond to each unit.
    int index_array[vint1[1], 4];
      index_array[,1] = vint7[1:vint1[1]];
      index_array[,2] = vint8[1:vint1[1]];
      index_array[,3] = vint9[1:vint1[1]];
      index_array[,4] = vint10[1:vint1[1]];

  // Initialize and compute log-likelihood
    real lp = 0;
    
    for (sp in 1:vint4[1]) {
      real lp_s = 0;
      if (vint5[sp] == 1) {
        for (i in 1:vint1[1]) {
          if (vint6[i] == sp) {
            int indices[vint2[i]] = index_array[i, 1:vint2[i]];
            if (vint3[i] == 1) {
              lp_s += bernoulli_logit_lpmf(1 | occ[i]);
              lp_s += bernoulli_logit_lpmf(y[indices] | mu[indices]);
            }
            if (vint3[i] == 0) {
              lp_s += log_sum_exp(bernoulli_logit_lpmf(1 | occ[i]) + 
                                    sum(log1m_inv_logit(mu[indices])), bernoulli_logit_lpmf(0 | occ[i]));
            }
          }
        }
        lp += log_inv_logit(Omega) + lp_s;
      } else {
        for (i in 1:vint1[1]) {
          if (vint6[i] == sp) {
            int indices[vint2[i]] = index_array[i, 1:vint2[i]];
            lp_s += log_sum_exp(bernoulli_logit_lpmf(1 | occ[i]) + 
                                  sum(log1m_inv_logit(mu[indices])), bernoulli_logit_lpmf(0 | occ[i]));
          }
        }
        lp += log_sum_exp(log1m_inv_logit(Omega), log_inv_logit(Omega) + lp_s);  
      }
    }
    return(lp);
  }

}
data {
  int<lower=1> N;  // total number of observations
  int<lower=1> N_unit;
  int Y[N];  // response variable
  // data for custom integer vectors
  int vint1[N];
  // data for custom integer vectors
  int vint2[N];
  // data for custom integer vectors
  int vint3[N];
  // data for custom integer vectors
  int vint4[N];
  // data for custom integer vectors
  int vint5[N];
  // data for custom integer vectors
  int vint6[N];
  // data for custom integer vectors
  int vint7[N];
  // data for custom integer vectors
  int vint8[N];
  // data for custom integer vectors
  int vint9[N];
  // data for custom integer vectors
  int vint10[N];
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  int<lower=1> K_occ;  // number of population-level effects
  matrix[N_unit, K_occ] X_occ;  // population-level design matrix
  // data for group-level effects of ID 1
  int<lower=1> N_1;  // number of grouping levels
  int<lower=1> M_1;  // number of coefficients per level
  int<lower=1> J_1[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_1_1;
  vector[N_unit] Z_1_occ_2;
  int<lower=1> NC_1;  // number of group-level correlations
  // data for group-level effects of ID 2
  int<lower=1> N_2;  // number of grouping levels
  int<lower=1> M_2;  // number of coefficients per level
  int<lower=1> J_2[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_2_1;
  vector[N] Z_2_2;
  // data for group-level effects of ID 3
  int<lower=1> N_3;  // number of grouping levels
  int<lower=1> M_3;  // number of coefficients per level
  int<lower=1> J_3[N_unit];  // grouping indicator per observation
  // group-level predictor values
  vector[N_unit] Z_3_occ_1;
  // data for group-level effects of ID 4
  int<lower=1> N_4;  // number of grouping levels
  int<lower=1> M_4;  // number of coefficients per level
  int<lower=1> J_4[N_unit];  // grouping indicator per observation
  // group-level predictor values
  vector[N_unit] Z_4_occ_1;
  vector[N_unit] Z_4_occ_2;
  vector[N_unit] Z_4_occ_3;
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  int Kc = K - 1;
  matrix[N, Kc] Xc;  // centered version of X without an intercept
  vector[Kc] means_X;  // column means of X before centering
  int Kc_occ = K_occ - 1;
  matrix[N_unit, Kc_occ] Xc_occ;  // centered version of X_occ without an intercept
  vector[Kc_occ] means_X_occ;  // column means of X_occ before centering
  for (i in 2:K) {
    means_X[i - 1] = mean(X[, i]);
    Xc[, i - 1] = X[, i] - means_X[i - 1];
  }
  for (i in 2:K_occ) {
    means_X_occ[i - 1] = mean(X_occ[, i]);
    Xc_occ[, i - 1] = X_occ[, i] - means_X_occ[i - 1];
  }
}
parameters {
  real<upper=16> Intercept_Omega;  // temporary intercept for centered predictors
  vector[Kc] b;  // population-level effects
  real Intercept;  // temporary intercept for centered predictors
  vector[Kc_occ] b_occ;  // population-level effects
  real Intercept_occ;  // temporary intercept for centered predictors
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  matrix[M_1, N_1] z_1;  // standardized group-level effects
  cholesky_factor_corr[M_1] L_1;  // cholesky factor of correlation matrix
  vector<lower=0>[M_2] sd_2;  // group-level standard deviations
  vector[N_2] z_2[M_2];  // standardized group-level effects
  vector<lower=0>[M_3] sd_3;  // group-level standard deviations
  vector[N_3] z_3[M_3];  // standardized group-level effects
  vector<lower=0>[M_4] sd_4;  // group-level standard deviations
  vector[N_4] z_4[M_4];  // standardized group-level effects
}
transformed parameters {
  matrix[N_1, M_1] r_1;  // actual group-level effects
  // using vectors speeds up indexing in loops
  vector[N_1] r_1_1;
  vector[N_1] r_1_occ_2;
  vector[N_2] r_2_1;  // actual group-level effects
  vector[N_2] r_2_2;  // actual group-level effects
  vector[N_3] r_3_occ_1;  // actual group-level effects
  vector[N_4] r_4_occ_1;  // actual group-level effects
  vector[N_4] r_4_occ_2;  // actual group-level effects
  vector[N_4] r_4_occ_3;  // actual group-level effects
  // compute actual group-level effects
  r_1 = scale_r_cor(z_1, sd_1, L_1);
  r_1_1 = r_1[, 1];
  r_1_occ_2 = r_1[, 2];
  r_2_1 = (sd_2[1] * (z_2[1]));
  r_2_2 = (sd_2[2] * (z_2[2]));
  r_3_occ_1 = (sd_3[1] * (z_3[1]));
  r_4_occ_1 = (sd_4[1] * (z_4[1]));
  r_4_occ_2 = (sd_4[2] * (z_4[2]));
  r_4_occ_3 = (sd_4[3] * (z_4[3]));
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] mu = Intercept + Xc * b;
    // initialize linear predictor term
    vector[N_unit] occ = Intercept_occ + Xc_occ * b_occ;
    // initialize linear predictor term
    real Omega = Intercept_Omega;
    for (n in 1:N) {
      // add more terms to the linear predictor
      mu[n] += r_1_1[J_1[n]] * Z_1_1[n] + r_2_1[J_2[n]] * Z_2_1[n] + r_2_2[J_2[n]] * Z_2_2[n];
    }
    for (n in 1:N_unit) {
      // add more terms to the linear predictor
      occ[n] += r_1_occ_2[J_1[n]] * Z_1_occ_2[n] + r_3_occ_1[J_3[n]] * Z_3_occ_1[n] + r_4_occ_1[J_4[n]] * Z_4_occ_1[n] + r_4_occ_2[J_4[n]] * Z_4_occ_2[n] + r_4_occ_3[J_4[n]] * Z_4_occ_3[n];
    }
    target += occupancy_augmented_lpmf(Y | mu, occ, Omega, vint1, vint2, vint3, vint4, vint5, vint6, vint7, vint8, vint9, vint10);
  }
  // priors including constants
  target += normal_lpdf(b[1] | 0, 1);
  target += normal_lpdf(b[2] | 0, 1);
  target += normal_lpdf(b[3] | 0, 1);
  target += logistic_lpdf(Intercept | 0, 1);
  target += normal_lpdf(b_occ[1] | 0, 5);
  target += normal_lpdf(b_occ[2] | 0, 5);
  target += normal_lpdf(b_occ[3] | 0, 5);
  target += logistic_lpdf(Intercept_occ | 0, 1);
  target += logistic_lpdf(Intercept_Omega | 0,1);
  target += normal_lpdf(sd_1[1] | 0, 5)
    - 1 * normal_lccdf(0 | 0, 5);
  target += normal_lpdf(sd_1[2] | 0, 5)
    - 1 * normal_lccdf(0 | 0, 5);
  target += std_normal_lpdf(to_vector(z_1));
  target += lkj_corr_cholesky_lpdf(L_1 | 1);
  target += normal_lpdf(sd_2[1] | 0, 5)
    - 1 * normal_lccdf(0 | 0, 5);
  target += normal_lpdf(sd_2[2] | 0, 5)
    - 1 * normal_lccdf(0 | 0, 5);
  target += std_normal_lpdf(z_2[1]);
  target += std_normal_lpdf(z_2[2]);
  target += normal_lpdf(sd_3 | 0, 5)
    - 1 * normal_lccdf(0 | 0, 5);
  target += std_normal_lpdf(z_3[1]);
  target += normal_lpdf(sd_4[1] | 0, 5)
    - 1 * normal_lccdf(0 | 0, 5);
  target += normal_lpdf(sd_4[2] | 0, 5)
    - 1 * normal_lccdf(0 | 0, 5);
  target += normal_lpdf(sd_4[3] | 0, 5)
    - 1 * normal_lccdf(0 | 0, 5);
  target += std_normal_lpdf(z_4[1]);
  target += std_normal_lpdf(z_4[2]);
  target += std_normal_lpdf(z_4[3]);
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept - dot_product(means_X, b);
  // actual population-level intercept
  real b_occ_Intercept = Intercept_occ - dot_product(means_X_occ, b_occ);
  // actual population-level intercept
  real b_Omega_Intercept = Intercept_Omega;
  // compute group-level correlations
  corr_matrix[M_1] Cor_1 = multiply_lower_tri_self_transpose(L_1);
  vector<lower=-1,upper=1>[NC_1] cor_1;
  // extract upper diagonal of correlation matrix
  for (k in 1:M_1) {
    for (j in 1:(k - 1)) {
      cor_1[choose(k - 1, 2) + j] = Cor_1[j, k];
    }
  }
}