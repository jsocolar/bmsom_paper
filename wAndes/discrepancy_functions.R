# Functions to compute discrepancy variables from simulated data replicates

##### Sum of Q #####
# Sum of detections, based on posterior predictive distribution that ties latent occupancy to data



##### Mackenzie-Bailey statistics #####
mb_counts <- function(histories, nv){
  dhs <- c("1111",
           "1000", "0100", "0010", "0001",
           "1100", "1010", "1001", "0110", "0101", "0011",
           "1110", "1101", "1011", "0111",
           "0000")
  dh_counts <- rep(0,16)
  for(i in 1:16){
    dh_counts[i] <- sum(nv == 4 & histories == dhs[i])
  }
  return(dh_counts)
}

##### Join count statistics on clusters #####
# function to summarize Q by cluster
cluster_q <- function(z_info, qdata){
  if(nrow(z_info) != length(qdata)){
    stop("lengths of z_info and qdata differ")
  }
  spcl <- 1:max(z_info$id_spCl)
  return(spcl %in% z_info$id_spCl[qdata == 1])
}

# function to get weights matrix from a slice of the birds dataframe (e.g. for a single species)
get_neighbor_weights <- function(bird_df, weights_function){
  n <- nrow(bird_df)
  distance_matrix <- geosphere::distm(bird_df[,c("lon", "lat")])
  weights_matrix <- weights_function(distance_matrix)
  diag(weights_matrix) <- 0
  return(weights_matrix)
}

# function to compute join-count statistic
join_counts_stat <- function(presence_vector, weights_matrix){
  if(!all(presence_vector %in% c(0,1))){
    stop("presence vector has entries that are not zero or one")
  }
  if(sum(diag(weights_matrix) != 0) != 0){
    stop("weights matrix has nonzero diagonal elements")
  }
  if(dim(weights_matrix)[1] != dim(weights_matrix)[2]){
    stop("weigts matrix is not square")
  }
  if(dim(weights_matrix)[1] != length(presence_vector)){
    stop("weights matrix dimensions not equal to presence vector length")
  }
  if(!isSymmetric(weights_matrix)){
    stop("weights matrix is not symmetric. join_counts_stat() currently requires a symmetric weights matrix")
  }
  
  W <- sum(weights_matrix)
  n <- length(presence_vector)
  nb <- sum(presence_vector)
  nw <- n - nb
  Ebw <- W*nb*nw/(n^2)
  
  p2 <- as.matrix(presence_vector)
  bw_matrix <- (!(p2 %*% t(p2))) * (!((1-p2) %*% (1-t(p2))))
  bw <- .5*sum(bw_matrix*weights_matrix)
  
  s1 <- .5*sum(weights_matrix^2)
  s2 <- sum((4*colSums(weights_matrix)^2))
  var <- .25*(2*s1*nb*nw/(n^2) + (s2 - 2*s1)*nb*nw*(nb+nw)/(n^3) + 4*(s1 - s2)*(nb^2)*(nw^2)/(n^4))
  
  return((bw - Ebw)/sqrt(var))
}



#   
# 
# # function to compute join-count statistic pa matrix
# join_counts_stat <- function(bird_df, weights_function, type = "simrep"){
#   njoins <- sum(distance_matrix > 0 & distance_matrix < threshold)/2
#   bb <- ww <- bw <- 0
#   if(type == "simrep"){
#     nb <- sum(bd_sp$cl_q)
#     nw <- n - nb
#     pb <- nb/n
#     pw <- nw/n
#     Ebb <- njoins*pb^2
#     Eww <- njoins*pw^2
#     Ebw <- njoins*2*pb*pw
#     for(i in 2:nrow(bd_sp)){
#       for(j in 1:(i-1)){
#         if(distance_matrix[i,j] < threshold){
#           if(bd_sp$cl_q[i] == 1 & bd_sp$cl_q[j] == 1){
#             bb <- bb+1
#           }else if(bd_sp$cl_q[i] == 0 & bd_sp$cl_q[j] == 0){
#             ww <- ww+1
#           }else{
#             bw <- bw + 1
#           }
#         }
#       }
#     }
#   }else if(type == "truevalue"){
#     nb <- sum(bd_sp$cl_q_real)
#     nw <- n - nb
#     pb <- nb/n
#     pw <- nw/n
#     Ebb <- njoins*pb^2
#     Eww <- njoins*pw^2
#     Ebw <- njoins*2*pb*pw
#     for(i in 2:nrow(bd_sp)){
#       for(j in 1:(i-1)){
#         if(distance_matrix[i,j] < threshold){
#           if(bd_sp$cl_q_real[i] == 1 & bd_sp$cl_q_real[j] == 1){
#             bb <- bb+1
#           }else if(bd_sp$cl_q_real[i] == 0 & bd_sp$cl_q_real[j] == 0){
#             ww <- ww+1
#           }else{
#             bw <- bw + 1
#           }
#         }
#       }
#     }
#   }
#   s1 <- 2*njoins
#   s3 <- 4*sum(apply(distance_matrix, MARGIN = 1, FUN = function(x){sum(x>0 & x < threshold)^2}))
#   var <- .25*( (s3-s1)*nb*(n-nb)/(n*(n-1)) + 4*(s1^2 - s3)*nb*(nb-1)*(n-nb)*(n-nb-1)/(n*(n-1)*(n-2)*(n-3))) - Ebw^2
#   return((bw - Ebw)/sqrt(var))
# }
# 
# species_list <- unique(birds$species)
# 
# # cluster effects included
# jcs <- jcs_forest <- jcs_pasture <- matrix(nrow=100, ncol=length(species_list))
# jcs_real <- jcs_real_forest <- jcs_real_pasture <- vector()
# for(i in 1:100){
#   print(i)
#   iter = 20*i
#   data_rep <- get_Z_probs(draws, iter, z_info, cluster_effect = "include", return_pdet = T)
#   data_rep$data_rep_margin <- rbinom(nrow(data_rep), 1, data_rep$psi*data_rep$pdet)
#   data_rep$id_spCl <- z_info$id_spCl
#   cq <- cluster_q(data_rep, colname = "data_rep_margin")
#   cq_real <- cluster_q(data_rep, colname = "Q")
#   for(j in 1:length(species_list)){
#     species <- species_list[j]
#     
#     bd2 <- data.frame(species = z_info$species, id_sp = z_info$id_sp, id_spCl = z_info$id_spCl, 
#                       lon = birds$lon, lat = birds$lat, pasture = birds$pasture,
#                       cl_q = cq[z_info$id_spCl],
#                       cl_q_real = cq_real[z_info$id_spCl])
#     bd3 <- bd2[!duplicated(z_info$id_spCl),]
#     
#     bd4_forest <- bd3[bd3$pasture == 0, ]
#     bd4_pasture <- bd3[bd3$pasture == 1, ]
#     
#     # combined
#     bd_sp <- bd3[bd3$species == species, ]
#     if(i == 1){
#       if(sum(bd_sp$cl_q_real) > 1){
#         jcs_real[j] <- join_counts_stat(bd_sp, 1e5, type = "truevalue")
#       }else{
#         jcs_real[j] <- NA
#       }
#     }
#     
#     if(sum(bd_sp$cl_q_real) > 1){
#       jcs[i,j] <- join_counts_stat(bd_sp, 1e5, type = "simrep")
#     }
#     
#     # forest
#     bd_sp <- bd4_forest[bd4_forest$species == species,]
#     if(i == 1){
#       if(sum(bd_sp$cl_q_real) > 1){
#         jcs_real_forest[j] <- join_counts_stat(bd_sp, 1e5, type = "truevalue")
#       }else{
#         jcs_real_forest[j] <- NA
#       }
#     }
#     
#     if(sum(bd_sp$cl_q_real) > 1){
#       jcs_forest[i,j] <- join_counts_stat(bd_sp, 1e5, type = "simrep")
#     }
#     
#     # pasture
#     bd_sp <- bd4_pasture[bd4_pasture$species == species,]
#     if(i == 1){
#       if(sum(bd_sp$cl_q_real) > 1){
#         jcs_real_pasture[j] <- join_counts_stat(bd_sp, 1e5, type = "truevalue")
#       }else{
#         jcs_real_pasture[j] <- NA
#       }
#     }
#     
#     if(sum(bd_sp$cl_q_real) > 1){
#       jcs_pasture[i,j] <- join_counts_stat(bd_sp, 1e5, type = "simrep")
#     }
#   }
# }