###############################################################################
###############################################################################
##  R code for Bayesian analysis of multiple paternity in bird and in mammals
##  Written by Hannah Correia @ Johns Hopkins University: hcorrei2@jhu.edu
##
##  Must run bird-mammal_MCMC_final.R first, and have output files
##  "MCMC_resids_mammals.rda"
##  "MCMC_resids_birds.rda" for each subset analysis 
##  (see NOTE, lines 69-84, in "bird-mammal_MCMC_final.R" for details)
##
##  Requires "k_alpha.R" code from Zapf et al. 2016 
##  (https://doi.org/10.1186/s12874-016-0200-9 see "Additional file 3" 
##  in the supplementary materials)
###############################################################################
###############################################################################


library(readxl)
library(tidyverse)
library(reshape2)
library(irr)
library(clusrank)
library(vcd)


#### Compare microsat DNA with DNA fingerprinting studies in birds ####

paternity_birds <- read.csv("paternity_birds.csv") ## (source data) 
# View(paternity_birds)


#### CMH test ####
### Cochran-Mantel-Haenszel test for DNA methods

## Create category for pmult
paternity_birds$pmult_cat <- ifelse(paternity_birds$pmult==0, "azero", "positive")

## Create category for brood size 
## integer with a "7+" group
paternity_birds$brood_id <- ifelse((paternity_birds$avgbrood - trunc(paternity_birds$avgbrood)) < 0.5, 
                                   floor(paternity_birds$avgbrood), ceiling(paternity_birds$avgbrood))
paternity_birds$brood_group <- ifelse(paternity_birds$brood_id >= 7, ">6.5",
                                      ifelse(paternity_birds$brood_id <=2, "<2.5",
                                             as.character(paternity_birds$brood_id)))
paternity_birds$brood_group <- ordered(paternity_birds$brood_group, 
                                       levels = c("<2.5","3","4","5","6",">6.5"))
## What are the sizes of each of these groups?
paternity_birds %>% filter(!is.na(pmult)) %>% group_by(brood_group) %>% summarise(count = n())

## Three-dimensional table 
ptab <- paternity_birds %>%
  filter(!is.na(pmult)) %>%
  count(DNA, pmult_cat, brood_group) %>%
  group_by(DNA) 

ptab5 <- xtabs(n ~ DNA + pmult_cat + brood_group, ptab) 

## CMH test
mantelhaen.test(ptab5, exact = TRUE, correct = FALSE)

## Conditional odds ratios
oddsratio(ptab5, 3, log = FALSE)

## Conditional log odds ratio plot
lor <- oddsratio(ptab5,3)
exp(confint(lor)) ## CI
summary(lor)
plot(lor, xlab="Brood size", ylab = "LOR(DNA / pmult)", main = "")



## need estimated pB for clusWilcox from separate Bayesian analyses 
load("analysis_microsat/MCMC_resids.rda") ## results from microsat DNA birds Bayesian analysis
MCMC_resids_microsat <- MCMC_resids

load("analysis_fingerprint/MCMC_resids.rda") ## results from fingerprint DNA birds Bayesian analysis
MCMC_resids_fingerprint <- MCMC_resids

MCMC_resids_microsat$DNA <- "microsatellite"
MCMC_resids_fingerprint$DNA <- "fingerprint"

birds_resids <- rbind(MCMC_resids_microsat, MCMC_resids_fingerprint)

## Create category for brood size 
## integer with a "7+" group
birds_resids$brood_id <- ifelse((birds_resids$avgbrood - trunc(birds_resids$avgbrood)) < 0.5, 
                                floor(birds_resids$avgbrood), ceiling(birds_resids$avgbrood))
birds_resids$brood_group <- ifelse(birds_resids$brood_id >= 7, ">6.5",
                                   ifelse(birds_resids$brood_id <=2, "<2.5",
                                          as.character(birds_resids$brood_id)))
birds_resids$brood_group <- ordered(birds_resids$brood_group, 
                                    levels = c("<2.5","3","4","5","6",">6.5"))

## clusWilcox doesn't play nice with characters or factors - convert to numeric
birds_resids$DNA_num <- ifelse(birds_resids$DNA=="microsatellite", 1, 2)
birds_resids$brood_grpnum <- as.integer(factor(birds_resids$brood_group))+1

## (expect microsat to find more MP > 0, therefore larger observed pmult in microsat data)
clusWilcox.test(pmult, cluster = brood_grpnum, group = DNA_num, data = birds_resids,
                alternative = "greater", method = "ds", exact = TRUE)

## (expect microsat to find more MP > 0, therefore smaller pB-p in microsat data)
## i.e., (pB-p)_microsat - (pB-p)_fingerprint < 0
clusWilcox.test(resid_pmult, cluster = brood_grpnum, group = DNA_num, data = birds_resids,
                alternative = "less", method = "ds", exact = TRUE)

## Difference in nsires between methods 
sire_ms <- birds_resids[birds_resids$DNA=="microsatellite",]$avgsire
sire_fp <- birds_resids[birds_resids$DNA=="fingerprint",]$avgsire

t.test(sire_ms, sire_fp, alternative = "two.sided", var.equal = FALSE)

## Difference in brood sizes between methods 
broodsz_ms <- birds_resids[birds_resids$DNA=="microsatellite",]$avgbrood
broodsz_fp <- birds_resids[birds_resids$DNA=="fingerprint",]$avgbrood

t.test(broodsz_ms, broodsz_fp, alternative = "two.sided", var.equal = FALSE)

## Difference in brood samples between methods (need source data)
samples_ms <- paternity_birds[paternity_birds$DNA=="microsatellite",]$nbrood
samples_fp <- paternity_birds[paternity_birds$DNA=="fingerprint",]$nbrood

t.test(samples_ms, samples_fp, alternative = "two.sided", var.equal = FALSE)




###  BASED ON THE ABOVE ANALYSES, WE CONCLUDE TO ONLY USE MICROSATELLITE STUDIES
###  FOR BIRDS AND MOVE ON TO COMPARING SOCIALLY MONOGAMOUS BIRDS WITH THOSE THAT
###  ARE NOT SOCIALLY MONOGAMOUS


#### Compare socially monogamous birds with birds that are NOT socially monogamous ####
paternity_birds <- read.csv("paternity_birds.csv") ## (source data)
# View(paternity_birds)


#### CMH test ####
### Cochran-Mantel-Haenszel test for Social Monogamy

## Create category for pmult
paternity_birds$pmult_cat <- ifelse(paternity_birds$pmult==0, "azero", "positive")

## Create category for brood size 
## integer with a "7+" group
paternity_birds$brood_id <- ifelse((paternity_birds$avgbrood - trunc(paternity_birds$avgbrood)) < 0.5, 
                                   floor(paternity_birds$avgbrood), ceiling(paternity_birds$avgbrood))
paternity_birds$brood_group <- ifelse(paternity_birds$brood_id >= 7, ">6.5",
                                      ifelse(paternity_birds$brood_id <=2, "<2.5",
                                             as.character(paternity_birds$brood_id)))
paternity_birds$brood_group <- ordered(paternity_birds$brood_group, 
                                       levels = c("<2.5","3","4","5","6",">6.5"))
## What are the sizes of each of these groups?
paternity_birds %>% filter(!is.na(pmult)) %>% group_by(brood_group) %>% summarise(count = n())

## SocMon is currently a 3-category factor - change to 2-category (mono/xmono only)
paternity_birds$SocMon_bin <- ifelse(paternity_birds$SocMon %in% c("mono","xmono"), 
                                     paternity_birds$SocMon, NA)

## Three-dimensional table 
ptab <- paternity_birds %>%
  filter(!is.na(pmult)) %>%
  count(SocMon_bin, pmult_cat, brood_group) %>%
  group_by(SocMon_bin) 

ptab5b <- xtabs(n ~ SocMon_bin + pmult_cat + brood_group, ptab) 

## CMH test
mantelhaen.test(ptab5b, exact = TRUE, correct = FALSE)

## Conditional odds ratios
oddsratio(ptab5b, 3, log = FALSE)

## Conditional log odds ratio plot
lor <- oddsratio(ptab5b,3)
exp(confint(lor)) ## CI
summary(lor)
plot(lor, xlab="Brood size", ylab = "LOR(SocMon / pmult)", main = "")



## need estimated pB for clusWilcox from separate Bayesian analyses 
load("analysis_mono/MCMC_resids.rda") ## results from monogamous birds Bayesian analysis
MCMC_resids_mono <- MCMC_resids

load("analysis_xmono/MCMC_resids.rda") ## results from non-monogamous birds Bayesian analysis
MCMC_resids_xmono <- MCMC_resids

MCMC_resids_mono$SocMon <- "mono"
MCMC_resids_xmono$SocMon <- "xmono"

birds_resids <- rbind(MCMC_resids_mono, MCMC_resids_xmono)


## Clustered Wilcoxon rank sum test for p and pB-p b/w Socially Monog AND NON SocMono
birds_resids$mono_num <- ifelse(birds_resids$SocMon=="xmono", 1, 2)
## Clusters matched to same as CMH test for socmono comparison
birds_resids$brood_id <- ifelse((birds_resids$avgbrood - trunc(birds_resids$avgbrood)) < 0.5,
                                floor(birds_resids$avgbrood), ceiling(birds_resids$avgbrood))
birds_resids$brood_group <- ifelse(birds_resids$brood_id >= 7, ">6.5",
                                   ifelse(birds_resids$brood_id <=2, "<2.5",
                                          as.character(birds_resids$brood_id)))
## clusWilcox doesn't play nice with characters or factors - convert to numeric
birds_resids$brood_group <- factor(birds_resids$brood_group)
birds_resids$brood_grpnum <- as.integer(factor(birds_resids$brood_group))+1

## (expect xmono to find more MP > 0, therefore larger observed pmult in xmono data)
clusWilcox.test(pmult, cluster = brood_grpnum, group = mono_num, data = birds_resids,
                alternative = "greater", method = "ds", exact = TRUE)

## (expect xmono to find more MP > 0, therefore smaller pB-p in xmono data)
## i.e., (pB-p)_xmono - (pB-p)_mono < 0
clusWilcox.test(resid_pmult, cluster = brood_grpnum, group = mono_num, data = birds_resids,
                alternative = "less", method = "ds", exact = TRUE)





#### Compare birds with mammals ####

###  BASED ON THE ABOVE ANALYSES, WE CONCLUDE TO ONLY USE MICROSATELLITE STUDIES
###  FOR BIRDS AND MOVE ON TO THE MCMC ANALYSIS OF BIRDS AND MAMMALS USING
###  POPULATIONS AS UNITS OF INTEREST, IN ORDER TO THEN COMPARE THE TWO TAXA.




#### Load results files from Bayesian MCMC analysis:

## birds (source data)
paternity_birds <- read.csv("paternity_birds.csv") 

## mammals (source data)
paternity_mammals <- read.csv("paternity_mammals.csv") 


#### Krippendorff's alpha ####
## Thanks to Zapf et al. 2016 (https://doi.org/10.1186/s12874-016-0200-9)
## for R code (see "Additional file 3" in the supplementary materials)
source("../k_alpha.R")

## Create unique species ID number for each species and obs number for each pmult within species
bird_dat <- paternity_birds %>%
  mutate(species2 = factor(species), species_ID2 = factor(paste0("S",factor(as.numeric(species2))))) %>%
  group_by(species_ID2) %>%
  mutate(obs_ID2 = factor(paste0("P", 1:length(pmult)))) %>%
  ungroup() 


#### Krippendorf's alpha for MP of birds
## Keep only species with more than one observation
bird_mp <- bird_dat %>% select(species_ID2, obs_ID2, pmult) %>%
  group_by(species_ID2) %>% filter(length(pmult)>1) %>% ungroup() %>%
  pivot_wider(names_from = species_ID2, values_from = pmult, id_cols = obs_ID2, 
              values_fill = NA)

mp_mat <- as.matrix(bird_mp[,-1])

k_alpha(t(mp_mat), alpha_q = 0.05, nboot = 1000, scaling = "ratio")


#### Krippendorff's alpha for brood size of birds
bird_bs <- bird_dat %>% select(species_ID2, obs_ID2, avgbrood) %>%
  group_by(species_ID2) %>% filter(length(avgbrood)>1) %>% ungroup() %>%
  # mutate(avgbrood = avgbrood/10) %>%
  pivot_wider(names_from = species_ID2, values_from = avgbrood, id_cols = obs_ID2, 
              values_fill = NA)

bs_mat <- as.matrix(bird_bs[,-1])

k_alpha(t(bs_mat), alpha_q = 0.05, nboot = 1000, scaling = "ratio")


#### Krippendorff's alpha for MP of mammals

## Create unique species ID number for each species and obs number for each pmult within species
mammal_dat <- paternity_mammals %>%
  mutate(species2 = factor(species), species_ID2 = factor(paste0("S",factor(as.numeric(species2))))) %>%
  group_by(species_ID2) %>%
  mutate(obs_ID2 = factor(paste0("P", 1:length(pmult)))) %>%
  ungroup() 

## Keep only species with more than one observation
mammal_mp <- mammal_dat %>% select(species_ID2, obs_ID2, pmult) %>%
  group_by(species_ID2) %>% filter(length(pmult)>1) %>% ungroup() %>%
  pivot_wider(names_from = species_ID2, values_from = pmult, id_cols = obs_ID2, 
              values_fill = NA)

mp_mat <- as.matrix(mammal_mp[,-1])

k_alpha(t(mp_mat), alpha_q = 0.05, nboot = 1000, scaling = "ratio")


#### Krippendorff's alpha for brood size of mammals
mammal_bs <- mammal_dat %>% select(species_ID2, obs_ID2, avgbrood) %>%
  group_by(species_ID2) %>% filter(length(avgbrood)>1) %>% ungroup() %>%
  # mutate(avgbrood = avgbrood/10) %>%
  pivot_wider(names_from = species_ID2, values_from = avgbrood, id_cols = obs_ID2, 
              values_fill = NA)

bs_mat <- as.matrix(mammal_bs[,-1])

k_alpha(t(bs_mat), alpha_q = 0.05, nboot = 1000, scaling = "ratio")


## birds (results data, microsat-only)
load("MCMC_resids_birds.rda")
bird_resids <- MCMC_resids ## results from microsat DNA birds Bayesian analysis

## mammals (results data)
load("MCMC_resids_mammals.rda")
mammal_resids <- MCMC_resids ## results from Bayesian analysis

bird_resids$group <- "Birds"
mammal_resids$group <- "Mammals"

## Combine mammals and birds results
resid_dat <- rbind(mammal_resids, bird_resids)


#### CMH test ####
### Cochran-Mantel-Haenszel test for DNA methods

## Create category for pmult
resid_dat$pmult_cat <- ifelse(resid_dat$pmult==0, "azero", "positive")

## integer with a "7+" group
resid_dat$brood_id <- ifelse((resid_dat$avgbrood - trunc(resid_dat$avgbrood)) < 0.5, 
                             floor(resid_dat$avgbrood), ceiling(resid_dat$avgbrood))
resid_dat$brood_group <- ifelse(resid_dat$brood_id >= 7, ">6.5",
                                ifelse(resid_dat$brood_id <=2, "<2.5",
                                       as.character(resid_dat$brood_id)))
resid_dat$brood_group <- ordered(resid_dat$brood_group, 
                                 levels = c("<2.5","3","4","5","6",">6.5"))


## Three-dimensional table 
ptab <- resid_dat %>%
  filter(!is.na(pmult)) %>%
  count(group, pmult_cat, brood_group) %>%
  group_by(group) 

ptab6 <- xtabs(n ~ group + pmult_cat + brood_group, ptab) 

## CMH test
mantelhaen.test(ptab6, exact = TRUE, correct = FALSE)

## Marginal odds ratios
AC <- margin.table(ptab6, c(2,1))

oddsratio(ptab6, 3, log = FALSE)

lor <- oddsratio(ptab6,3)
exp(confint(lor)) ## CI
summary(lor)
plot(lor, xlab="Brood/litter size", ylab = "LOR(taxa / pmult)", main = "")



#### Non-parametric comparison of mammals and birds in binned broodsizes ####
#### Rank-sum test for clustered data, Datta and Satten 2005 

## "ds" method (i.e., Datta and Satten, informative clusters) does not behave with
## groups that are factors or characters; need numeric groups
resid_dat$group_num <- ifelse(resid_dat$group=="Birds", 1, 2) ## bird = 1, mammals = 2

## Clusters matched to same as CMH test for DNA method comparison
resid_dat$brood_id <- ifelse((resid_dat$avgbrood - trunc(resid_dat$avgbrood)) < 0.5, 
                             floor(resid_dat$avgbrood), ceiling(resid_dat$avgbrood))
resid_dat$brood_group <- ifelse(resid_dat$brood_id >= 7, ">6.5",
                                ifelse(resid_dat$brood_id <=2, "<=2",
                                       as.character(resid_dat$brood_id)))
## clusWilcox doesn't play nice with characters or factors - convert to numeric
resid_dat$brood_group <- factor(resid_dat$brood_group)
resid_dat$brood_grpnum <- as.integer(factor(resid_dat$brood_group))+1

## "Are values of pmult for birds (group_num=1) signif smaller than those of mammals (group_num=2)?"
## i.e. pmult_birds - pmult_mammals < 0
clusWilcox.test(pmult, cluster = brood_grpnum, group = group_num, data = resid_dat, 
                alternative = "less", method = "ds", exact  = TRUE)

## "Are values of pB-p for birds (group_num=1) signif larger than those of mammals (group_num=2)?"
## i.e. (pB-p)_birds - (pB-p)_mammals > 0
clusWilcox.test(resid_pmult, cluster = brood_grpnum, group = group_num, data = resid_dat, 
                alternative = "greater", method = "ds", exact  = TRUE)


## Difference in nsires between methods 
sires_b <- resid_dat[resid_dat$group=="Birds",]$avgsire
sires_m <- resid_dat[resid_dat$group=="Mammals",]$avgsire

t.test(sires_b, sires_m, alternative = "two.sided", var.equal = FALSE)

## Difference in brood sizes between methods 
broodsz_b <- resid_dat[resid_dat$group=="Birds",]$avgbrood
broodsz_m <- resid_dat[resid_dat$group=="Mammals",]$avgbrood

t.test(broodsz_b, broodsz_m, alternative = "two.sided", var.equal = FALSE)

## Difference in brood samples between methods (need source data)
samples_b <- mammal_dat$nbrood
samples_m <- bird_dat[bird_dat$DNA=="microsatellite",]$nbrood

t.test(samples_b, samples_m, alternative = "two.sided", var.equal = FALSE)



