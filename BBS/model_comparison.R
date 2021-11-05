library(brms)
library(loo)
setwd("/Users/jacob/Dropbox/Work/Occupancy/biogeographicMSOM")

bbs_naive <- readRDS("stan_outputs/brms_bbs/naive.RDS")
bbs_dist <- readRDS("stan_outputs/brms_bbs/dist_final.RDS")
bbs_clip <- readRDS("stan_outputs/brms_bbs/dist_clip_400K_final.RDS")

bbs_naive_loo <- add_criterion(bbs_naive, "loo", moment_match = T, cores = 4)
bbs_dist_loo <- add_criterion(bbs_dist, "loo", moment_match = T, cores = 4)
saveRDS(bbs_naive_loo, "stan_outputs/brms_bbs/naive_loo.RDS")
saveRDS(bbs_dist_loo, "stan_outputs/brms_bbs/dist_loo.RDS")

bbs_naive_loo <- readRDS("stan_outputs/brms_bbs/naive_loo.RDS")
bbs_dist_loo <- readRDS("stan_outputs/brms_bbs/dist_loo.RDS")

species <- unique(bbs_naive$data$species)

naive_species_loo2 <- dist_species_loo2 <- list()
for (i in 1:length(species)) {
  print(i)
  naive_species_loo2[[i]] <- add_criterion(bbs_naive, "loo", moment_match = T, cores = 4, 
                                           newdata = bbs_naive$data[bbs_naive$data$species == species[i], ])
  dist_species_loo2[[i]] <- add_criterion(bbs_dist, "loo", moment_match = T, cores = 4, 
                                          newdata = bbs_dist$data[bbs_dist$data$species == species[i], ])
}

saveRDS(naive_species_loo2, "stan_outputs/brms_bbs/naive_species_loo2.RDS")
saveRDS(dist_species_loo2, "stan_outputs/brms_bbs/dist_species_loo2.RDS")

naive_species_loo2 <- readRDS("stan_outputs/brms_bbs/naive_species_loo2.RDS")
dist_species_loo2 <- readRDS("stan_outputs/brms_bbs/dist_species_loo2.RDS")

brms::loo_compare(bbs_naive_loo, bbs_dist_loo)
a <- 0
naive_dist_comp <- naive_dist_comp_se <- vector()
for (i in 1:51) {
  if(!identical(grepl("dist", row.names(brms::loo_compare(naive_species_loo2[[i]], dist_species_loo2[[i]]))), c(T,F))) {
    print(i)
    print(species[i])
    print(brms::loo_compare(naive_species_loo2[[i]], dist_species_loo2[[i]]))
    a <- a - brms::loo_compare(naive_species_loo2[[i]], dist_species_loo2[[i]])[2,1]
    naive_dist_comp[i] <- -brms::loo_compare(naive_species_loo2[[i]], dist_species_loo2[[i]])[2,1]
  } else {
    a <- a + brms::loo_compare(naive_species_loo2[[i]], dist_species_loo2[[i]])[2,1]
    naive_dist_comp[i] <- brms::loo_compare(naive_species_loo2[[i]], dist_species_loo2[[i]])[2,1]
  }
  naive_dist_comp_se[i] <- brms::loo_compare(naive_species_loo2[[i]], dist_species_loo2[[i]])[2,2]
}
a
o_ndc <- order(naive_dist_comp)
naive_dist_comp <- naive_dist_comp[o_ndc]
naive_dist_comp_se <- naive_dist_comp_se[o_ndc]

###########
bbs_dc_loo <- add_criterion(bbs_dist, "loo", moment_match = T, cores = 4, newdata = bbs_clip$data)
bbs_clip_loo <- add_criterion(bbs_clip, "loo", moment_match = T, cores = 4)
saveRDS(bbs_dc_loo, "stan_outputs/brms_bbs/naive_loo.RDS")
saveRDS(bbs_dist_loo, "stan_outputs/brms_bbs/dist_loo.RDS")

bbs_dc_loo <- readRDS("stan_outputs/brms_bbs/naive_loo.RDS")
bbs_dist_loo <- readRDS("stan_outputs/brms_bbs/dist_loo.RDS")

dc_species_loo2 <- clip_species_loo2 <- list()
for (i in 1:length(species)) {
  print(i)
  dc_species_loo2[[i]] <- add_criterion(bbs_dist, "loo", moment_match = T, cores = 4, 
                                        newdata = bbs_clip$data[bbs_clip$data$species == species[i], ])
  clip_species_loo2[[i]] <- add_criterion(bbs_clip, "loo", moment_match = T, cores = 4, 
                                          newdata = bbs_clip$data[bbs_clip$data$species == species[i], ])
}

saveRDS(dc_species_loo2, "stan_outputs/brms_bbs/dc_species_loo2.RDS")
saveRDS(clip_species_loo2, "stan_outputs/brms_bbs/clip_species_loo2.RDS")

dc_species_loo2 <- readRDS("stan_outputs/brms_bbs/dc_species_loo2.RDS")
clip_species_loo2 <- readRDS("stan_outputs/brms_bbs/clip_species_loo2.RDS")

brms::loo_compare(bbs_dc_loo, bbs_clip_loo)
a <- 0
dist_clip_comp <- dist_clip_comp_se <- vector()
for (i in 1:51) {
  if(!identical(grepl("clip", row.names(brms::loo_compare(dc_species_loo2[[i]], clip_species_loo2[[i]]))), c(T,F))) {
    print(i)
    print(species[i])
    print(brms::loo_compare(dc_species_loo2[[i]], clip_species_loo2[[i]]))
    a <- a - brms::loo_compare(dc_species_loo2[[i]], clip_species_loo2[[i]])[2,1]
    dist_clip_comp[i] <- -brms::loo_compare(dc_species_loo2[[i]], clip_species_loo2[[i]])[2,1]
  } else {
    a <- a + brms::loo_compare(dc_species_loo2[[i]], clip_species_loo2[[i]])[2,1]
    dist_clip_comp[i] <- brms::loo_compare(dc_species_loo2[[i]], clip_species_loo2[[i]])[2,1]
  }
  dist_clip_comp_se[i] <- brms::loo_compare(dc_species_loo2[[i]], clip_species_loo2[[i]])[2,2]
}

o_dcc <- order(dist_clip_comp)
dist_clip_comp <- dist_clip_comp[o_dcc]
dist_clip_comp_se <- dist_clip_comp_se[o_dcc]

par(mfrow = c(1,2))
plot(-naive_dist_comp, pch = 16, xaxt = "n", xlab = "", ylab = "",
     ylim = c(0,600))
mtext("ELPD difference", side = 2, line = 2.8, cex = 1.5)
for(i in 1:length(dist_clip_comp)){
  lines(x = c(i,i), y = c(-naive_dist_comp[i] + 2*naive_dist_comp_se[i],
                          -naive_dist_comp[i] - 2*naive_dist_comp_se[i]))
}
abline(h = 0)
title(main="bMSOM vs traditional MSOM", cex.main=1)
plot(-dist_clip_comp, pch = 16, xaxt = "n", xlab = "", ylab = "",
     ylim = c(-3,60))
for(i in 1:length(dist_clip_comp)){
  lines(x = c(i,i), y = c(-dist_clip_comp[i] + 2*dist_clip_comp_se[i],
                          -dist_clip_comp[i] - 2*dist_clip_comp_se[i]))
}
abline(h = 0)
title(main="bMSOM: clipped vs unclipped", cex.main = 1)

