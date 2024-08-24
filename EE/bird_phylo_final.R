################################################################################
#### Phylogenetic models for birds and mammals
#### Written by Hannah Correia | Johns Hopkins University | hcorrei2@jhu.edu
#### Last updated 24 August 2024
################################################################################

#### Load libraries
library(TreeTools)
library(phytools)

library(tidyverse)
library(tidybayes)
library(brms)
library(MCMCglmm) # for phylo relationship matrices only
# library(coda) # for better plots in parallel

# library(insight) # for variance extraction (not working since version >0.20.1)


### Load combined data and Amniota phylogenetic tree
load("combinedMB_data.rda")

## Read in supertree
supertree <- read.nexus("combined_tree.nex")

## Compute branch lengths 
combined_phylo <- compute.brlen(supertree, method="coalescent")
# is.ultrametric(combined_phylo)
# plot(combined_phylo, cex = 0.3)

## Manually name the nodes, to prevent warnings
combined_phylo <- makeNodeLabel(combined_phylo, "number")

## Snowy plover still considered Kentish plover based on most current phylo tree
## Therefore, have to rename "phylo_name" of snowy plover to Kentish plover phylo name
combinedMB_data[combinedMB_data$common_name=="Snowy plover",]$phylo_name <- 
  unique(combinedMB_data[combinedMB_data$common_name=="Kentish plover",]$phylo_name)

## Check all data species are now represented in phylo
phylo.names <- sort(unique(combined_phylo$tip.label))
data.names <- sort(unique(combinedMB_data$phylo_name))
identical(data.names, phylo.names) ## Do they match?
# ## If FALSE, which species in the data file are not in the phylogeny?
# data.names[which(!data.names %in% phylo.names)]  ## Species in the data file that are not on the phylo
# phylo.names[which(!phylo.names %in% data.names)]  ## Species on the phylo that are not in the data file

## If all species in the phylo are accounted for in the data:
combined_data <- combinedMB_data

## If species don't match, remove extras from data file so they do:
## (doesn't matter if there are extra species in the phylogeny file)
# drop.species <- c(data.names[which(!data.names %in% phylo.names)], 
#                   phylo.names[which(!phylo.names %in% data.names)])
# combined_data <- combinedMB_data[!combinedMB_data$phylo_name %in% drop.species,]

## Add extra column for models accounting for multiple populations per species
combined_data$phylo_name2 <- combined_data$phylo_name


## Prep Amniota tree - phylogenetic relationship matrix
inv.phylo <- MCMCglmm::inverseA(combined_phylo, nodes = "TIPS", scale = TRUE)
AC <- solve(inv.phylo$Ainv)
rownames(AC) <- rownames(inv.phylo$Ainv)


#### Load birds and mammals phylogenetic trees ####

## Load and prep bird tree and phylogenetic relationship matrix
bird_tree <- read.nexus("119-spe.nex")
consensusBTree <- Consensus(bird_tree, p = 1, check.labels = TRUE)
bird_phylo <- compute.brlen(consensusBTree, method="coalescent")
# is.ultrametric(bird_phylo)
inv.phylo <- MCMCglmm::inverseA(bird_phylo, nodes = "TIPS", scale = TRUE)
AB <- solve(inv.phylo$Ainv)
rownames(AB) <- rownames(inv.phylo$Ainv)

## Load and prep mammal tree and phylogenetic relationship matrix
mammal_tree <- read.nexus("61-Mamm.nex")
consensusMTree <- Consensus(mammal_tree, p = 1, check.labels = TRUE)
mamm_phylo <- compute.brlen(consensusMTree, method="coalescent")
# is.ultrametric(mamm_phylo)
inv.phylo <- MCMCglmm::inverseA(mamm_phylo, nodes = "TIPS", scale = TRUE)
AM <- solve(inv.phylo$Ainv)
rownames(AM) <- rownames(inv.phylo$Ainv)




#### Phylogenetic correlation of MP in birds ####

## (no avgbrood in predictors)
mp.bird.brm <- brm(round(pmult*nbrood) | trials(nbrood) ~ 
                     (1 | gr(phylo_name, cov = AB)) + (1 | phylo_name2), 
                   data = combined_data[combined_data$taxa=="bird",],
                   family = "binomial",
                   data2 = list(AB = AB), 
                   sample_prior = "yes",
                   silent = 2, refresh = 0, 
                   seed = 6765)

## Phylogenetic ICC using `insight` (not working in `insight`>0.20.1):
# var_comp <- get_variance(mp.bird.brm)
# tau2_phy <- var_comp$var.intercept[[1]] # extract first component, which is the phylogenetic variance
# tau2 <- var_comp$var.random
# sig2 <- var_comp$var.residual
# (ICC_mp_bird <- tau2_phy/(tau2 + sig2))

## Calculate a mean and CI for phylo ICC using posterior draws from model:
sig2_logistic <- (pi^2)/3  ## resid var approximation for logistic regression; see Nakagawa & Schielzeth (2010)
lambda <- as.mcmc(as_draws_matrix(mp.bird.brm, variable = "sd_phylo_name__Intercept")^2)/
  (as.mcmc(as_draws_matrix(mp.bird.brm, "sd_phylo_name__Intercept")^2) +
     as.mcmc(as_draws_matrix(mp.bird.brm, "sd_phylo_name2__Intercept")^2) +
     sig2_logistic)
mean(lambda)
# HPDinterval(lambda)
quantile(lambda, c(0.025, 0.975))

## Calculate a mean and CI for phylo ICC using hypothesis test from brms package:
## sig2 = (pi^2)/3 = 3.289868 approximation for logistic regression
hyp <- "sd_phylo_name__Intercept^2 / (sd_phylo_name__Intercept^2 + sd_phylo_name2__Intercept^2 + 3.289868) = 0"
(hyp_mp_bird <- hypothesis(mp.bird.brm, hyp, class = NULL))



#### Phylogenetic correlation of MP in mammals ####

mp.mamm.brm <- brm(round(pmult*nbrood) | trials(nbrood) ~ 
                     (1 | gr(phylo_name, cov = AM)) + (1 | phylo_name2), 
                   data = combined_data[combined_data$taxa=="mammal",],
                   family = "binomial",
                   data2 = list(AM = AM), 
                   sample_prior = "yes",
                   silent = 2, refresh = 0, 
                   seed = 6765)

## Phylogenetic ICC using `insight` (not working in `insight`>0.20.1):
# var_comp <- get_variance(mp.mamm.brm)
# tau2_phy <- var_comp$var.intercept[[1]] # extract first component, which is the phylogenetic variance
# tau2 <- var_comp$var.random
# sig2 <- var_comp$var.residual
# (ICC_mp_mamm <- tau2_phy/(tau2 + sig2))

## Calculate a mean and CI for phylo ICC using posterior draws from model:
sig2_logistic <- (pi^2)/3  ## resid var approximation for logistic regression; see Nakagawa & Schielzeth (2010)
lambda <- as.mcmc(as_draws_matrix(mp.mamm.brm, variable = "sd_phylo_name__Intercept")^2)/
  (as.mcmc(as_draws_matrix(mp.mamm.brm, "sd_phylo_name__Intercept")^2) +
     as.mcmc(as_draws_matrix(mp.mamm.brm, "sd_phylo_name2__Intercept")^2) +
     sig2_logistic)
mean(lambda)
# HPDinterval(lambda)
quantile(lambda, c(0.025, 0.975))

## Calculate a mean and CI for phylo ICC using hypothesis test from brms package:
## sig2 = (pi^2)/3 = 3.289868 approximation for logistic regression
hyp <- "sd_phylo_name__Intercept^2 / (sd_phylo_name__Intercept^2 + sd_phylo_name2__Intercept^2 + 3.289868) = 0"
(hyp_mp_mamm <- hypothesis(mp.mamm.brm, hyp, class = NULL))



#### Phylogenetic correlation of resid MP in birds ####

resid.bird.brm <- brm(resid_pmult ~ (1 | gr(phylo_name, cov = AB)) + (1 | phylo_name2), 
                      data = combined_data[combined_data$taxa=="bird",],
                      family = "gaussian",
                      data2 = list(AB = AB), 
                      iter = 3000,
                      sample_prior = "yes",
                      silent = 2, refresh = 0, 
                      seed = 6765)

## Calculate a mean and CI for phylo ICC using hypothesis test from brms package:
hyp <- "sd_phylo_name__Intercept^2 / (sd_phylo_name__Intercept^2 + sd_phylo_name2__Intercept^2 + sigma^2) = 0"
(hyp_resid_bird <- hypothesis(resid.bird.brm, hyp, class = NULL))



#### Phylogenetic correlation of resid MP in mammals ####

resid.mamm.brm <- brm(resid_pmult ~ (1 | gr(phylo_name, cov = AM)) + (1 | phylo_name2), 
                      data = combined_data[combined_data$taxa=="mammal",],
                      family = "gaussian",
                      data2 = list(AM = AM), 
                      iter = 5000,
                      control = list(adapt_delta = 0.95),
                      sample_prior = "yes",
                      silent = 2, refresh = 0, 
                      seed = 6765)

## Calculate a mean and CI for phylo ICC using hypothesis test from brms package:
hyp <- "sd_phylo_name__Intercept^2 / (sd_phylo_name__Intercept^2 + sd_phylo_name2__Intercept^2 + sigma^2) = 0"
(hyp <- hypothesis(resid.mamm.brm, hyp, class = NULL))



#### Phylogenetic correlation of MP in combined tree ####

mp.both.brm <- brm(round(pmult*nbrood) | trials(nbrood) ~ 
                     (1 | gr(phylo_name, cov = AC)) + (1 | phylo_name2), 
                   data = combined_data,
                   family = "binomial",
                   data2 = list(AC = AC), 
                   sample_prior = "yes",
                   silent = 2, refresh = 0, 
                   seed = 6765)

## Phylogenetic ICC using `insight` (not working in `insight`>0.20.1):
# var_comp <- get_variance(mp.both.brm)
# tau2_phy <- var_comp$var.intercept[[1]] # extract first component, which is the phylogenetic variance
# tau2 <- var_comp$var.random
# sig2 <- var_comp$var.residual
# (ICC_mp_both <- tau2_phy/(tau2 + sig2))

## Calculate a mean and CI for phylo ICC using hypothesis test from brms package:
## sig2 = (pi^2)/3 = 3.289868 approximation for logistic regression
hyp <- "sd_phylo_name__Intercept^2 / (sd_phylo_name__Intercept^2 + sd_phylo_name2__Intercept^2 + 3.289868) = 0"
(hyp_mp_both <- hypothesis(mp.both.brm, hyp, class = NULL))



#### Phylogenetic correlation of resid MP in combined tree ####

resid.both.brm <- brm(resid_pmult ~ (1 | gr(phylo_name, cov = AC)) + (1 | phylo_name2), 
                      data = combined_data,
                      family = "gaussian",
                      data2 = list(AC = AC), 
                      iter = 4000,
                      sample_prior = "yes",
                      silent = 2, refresh = 0, 
                      seed = 6765)

## Calculate a mean and CI for phylo ICC using hypothesis test from brms package:
hyp <- "sd_phylo_name__Intercept^2 / (sd_phylo_name__Intercept^2 + sd_phylo_name2__Intercept^2 + sigma^2) = 0"
(hyp <- hypothesis(resid.both.brm, hyp, class = NULL))



#### Phylogenetic correlation of body mass in birds ####
bmass.bird.brm <- brm(Bmass ~ (1 | gr(phylo_name, cov = AB)) + (1 | phylo_name2), 
                      data = combined_data[combined_data$taxa=="bird",],
                      family = "exponential",
                      data2 = list(AB = AB), 
                      sample_prior = "yes",
                      silent = 2, refresh = 0, 
                      seed = 6765)

## Phylogenetic ICC using `insight` (not working in `insight`>0.20.1):
# var_comp <- get_variance(bmass.bird.brm)
# tau2_phy <- var_comp$var.intercept[[1]] # extract first component, which is the phylogenetic variance
# tau2 <- var_comp$var.random
# sig2 <- var_comp$var.residual
# (ICC_bmass_bird <- tau2_phy/(tau2 + sig2))

## Calculate a mean and CI for phylo ICC using posterior draws from model:
sig2_exp <- 1  ## brms normalizes the exp variance to a standard scale; variance set to 1
lambda <- as.mcmc(as_draws_matrix(bmass.bird.brm, variable = "sd_phylo_name__Intercept")^2)/
  (as.mcmc(as_draws_matrix(bmass.bird.brm, "sd_phylo_name__Intercept")^2) +
     as.mcmc(as_draws_matrix(bmass.bird.brm, "sd_phylo_name2__Intercept")^2) +
     sig2_exp)
mean(lambda)
# HPDinterval(lambda)
quantile(lambda, c(0.025, 0.975))

## Calculate a mean and CI for phylo ICC using hypothesis test from brms package:
## sig2_exp = 1 b/c brms normalizes the exp variance to a standard scale
hyp <- "sd_phylo_name__Intercept^2 / (sd_phylo_name__Intercept^2 + sd_phylo_name2__Intercept^2 + 1) = 0"
(hyp_bmass_bird <- hypothesis(bmass.bird.brm, hyp, class = NULL))



#### Phylogenetic correlation of body mass in mammals ####
bmass.mamm.brm <- brm(Bmass ~ (1 | gr(phylo_name, cov = AM)) + (1 | phylo_name2), 
                      data = combined_data[combined_data$taxa=="mammal",],
                      family = "exponential",
                      data2 = list(AM = AM), 
                      sample_prior = "yes",
                      silent = 2, refresh = 0, 
                      seed = 6765)

## Phylogenetic ICC using `insight` (not working in `insight`>0.20.1):
# var_comp <- get_variance(bmass.mamm.brm)
# tau2_phy <- var_comp$var.intercept[[1]] # extract first component, which is the phylogenetic variance
# tau2 <- var_comp$var.random
# sig2 <- var_comp$var.residual
# (ICC_bmass_mamm <- tau2_phy/(tau2 + sig2))

## Calculate a mean and CI for phylo ICC using posterior draws from model:
sig2_exp <- 1  ## brms normalizes the exp variance to a standard scale; variance set to 1
lambda <- as.mcmc(as_draws_matrix(bmass.mamm.brm, variable = "sd_phylo_name__Intercept")^2)/
  (as.mcmc(as_draws_matrix(bmass.mamm.brm, "sd_phylo_name__Intercept")^2) +
     as.mcmc(as_draws_matrix(bmass.mamm.brm, "sd_phylo_name2__Intercept")^2) +
     sig2_exp)
mean(lambda)
# HPDinterval(lambda)
quantile(lambda, c(0.025, 0.975))

## Calculate a mean and CI for phylo ICC using hypothesis test from brms package:
## sig2_exp = 1 b/c brms normalizes the exp variance to a standard scale
hyp <- "sd_phylo_name__Intercept^2 / (sd_phylo_name__Intercept^2 + sd_phylo_name2__Intercept^2 + 1) = 0"
(hyp_bmass_mamm <- hypothesis(bmass.mamm.brm, hyp, class = NULL))




#### Compare monogamous and non-monogamous birds ####

#### Is MP different between mono and xmono birds when accounting for brood size? ####
mp.mono.sbrood.brm <- brm(round(pmult*nbrood) | trials(nbrood) ~ I(log10(avgbrood)) + social_monogamy +
                            (1 | gr(phylo_name, cov = AB)) + (1 | phylo_name2),
                          ## Remove birds with no info on monogamy
                          data = combined_data[combined_data$taxa=="bird" & 
                                                 !combined_data$social_monogamy=="no_info",], 
                          family = "binomial",
                          data2 = list(AB = AB), 
                          sample_prior = "yes",
                          control = list(adapt_delta = 0.9),
                          silent = 2, refresh = 0, 
                          seed = 6765)

## Is "social_monogamy" term significant?
summary(mp.mono.sbrood.brm)




#### Is pB-p different between mono and xmono birds when accounting for brood size? ####
resid.mono.sbrood.brm <- brm(resid_pmult ~ I(log10(avgbrood)) + social_monogamy +
                               (1 | gr(phylo_name, cov = AB)) + (1 | phylo_name2),
                             ## Remove birds with no info on monogamy
                             data = combined_data[combined_data$taxa=="bird" & 
                                                    !combined_data$social_monogamy=="no_info",], 
                             family = "gaussian",
                             data2 = list(AB = AB), 
                             sample_prior = "yes",
                             iter = 3000,
                             silent = 2, refresh = 0, 
                             seed = 6765)

## Is "social_monogamy" term significant?
summary(resid.mono.sbrood.brm)




#### Compare birds and mammals ####

#### Is MP different between birds and mammals when accounting for brood size? ####
mp.both.sbrood.brm <- brm(round(pmult*nbrood) | trials(nbrood) ~ I(log10(avgbrood)) + taxa +
                            (1 | gr(phylo_name, cov = AC)) + (1 | phylo_name2),
                          data = combined_data,
                          family = "binomial",
                          data2 = list(AC = AC), 
                          sample_prior = "yes",
                          silent = 2, refresh = 0, 
                          seed = 6765)

## Is "taxa" term significant?
summary(mp.both.sbrood.brm)




#### Is pB-p different between birds and mammals when accounting for brood size? ####
resid.both.sbrood.brm <- brm(resid_pmult ~ I(log10(avgbrood)) + taxa +
                               (1 | gr(phylo_name, cov = AC)) + (1 | phylo_name2),
                             data = combined_data,
                             family = "gaussian",
                             data2 = list(AC = AC), 
                             sample_prior = "yes",
                             iter = 4000,
                             silent = 2, refresh = 0, 
                             seed = 6765)

## Is "taxa" term significant?
summary(resid.both.sbrood.brm)




