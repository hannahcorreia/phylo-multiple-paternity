###############################################################################
###############################################################################
##  R code for Bayesian analysis of multiple paternity in bird and in mammals
##  Written by Hannah Correia @ Johns Hopkins University: hcorrei2@jhu.edu
###############################################################################
###############################################################################

#### **Bayesian analysis of birds** ####

## For large nbroods, will likely run into "vector memory exhausted" error
## On smaller Mac systems, 
## Open terminal,and type in the following three lines:
## cd ~
## touch .Renviron
## open .Renviron
## Save the following as the first line of .Renviron:
## R_MAX_VSIZE=100Gb 
## Note: This limit includes both physical and virtual memory.



# clear the workspace
rm(list = ls(all = T))

library(rjags)
library(R2OpenBUGS)
library(coda)
library(extraDistr)
library(ggplot2)

##### post.summ function courtesy of Ben Staton
post.summ = function(post, var) {
  
  # coerce to matrix for easy subsetting
  post.samp = as.matrix(post)
  
  # if parameter is indexed
  if(substr(var, nchar(var), nchar(var)) == "[") {
    # extract columns with headers equal to the desired variable
    post.sub = post.samp[,substr(colnames(post.samp), 1, nchar(var)) == var]
    # calculate desired quantities
    summ = apply(post.sub, 2, function(x) c(mean = mean(x), sd = sd(x), quantile(x, c(0.5, 0.025, 0.975))))
    return(summ)
  }
  
  # if parameter is not indexed
  if(substr(var, nchar(var), nchar(var)) != "[") {
    # extract the column with the same header as the desired variable
    post.sub = post.samp[,substr(colnames(post.samp), 1, nchar(var)) == var]
    # calculate the desired quantities
    summ = c(mean = mean(post.sub), sd = sd(post.sub), quantile(post.sub, c(0.5, 0.025, 0.975)))
    return(summ)
  }
}

#########################################################
#### ESTIMATION OF Q (RATE OF SUCCESS OF EACH SIRE) 
#### AND CALCULATION OF M (# MATES) FROM R (# SIRES)
#########################################################

##### read and prepare data #####
dat1 <- read.csv("paternity_birds.csv")

## remove birds with nbrood = NA and avgbrood = NA (can't run without avgbrood)
dat <- dat1[!is.na(dat1$nbrood) & !is.na(dat1$avgbrood),]


###############################################################################
##  NOTE:
##  The following code runs through ALL bird species together
##  In order to replicate analyses correctly, the following code must
##  be run on each subsample of bird species. 
##  E.g., for comparing socially monogamous birds to socially non-monogamous 
##  birds, must subset data to only species that are socially monogamous, 
##  run the code, then subset the data to only species that are socially
##  non-monogamous, then run the code again.
##
##  Separate files/folders for each subset analysis are recommended.
###############################################################################

## Uncomment one of the following depending on subset of birds being analyzed
# dat <- dat[dat$DNA=="microsatellite",]  ## for microsatellite DNA populations
# dat <- dat[dat$DNA=="fingerprint",]  ## for fingerprinting populations
# dat <- dat[dat$SocMon=="mono",]  ## for socially monogamous populations (microsat DNA only)
# dat <- dat[dat$SocMon=="xmono",]  ## for socially non-monogamous populations (microsat DNA only)




# expand data to create multiple broods within each species (sample size * 10)
for(i in 1:nrow(dat)){
  temp.sp <- expand.grid(species = rep(dat[i,]$species, dat[i,]$nbrood*10))
  if(i == 1){
    temp <- temp.sp
  } else {
    temp <- rbind(temp, temp.sp)
  } 
}
mates.df <- merge(temp, dat, all.x = TRUE, all.y = FALSE)


##### specify model code #####
model_byspecies_q.file <- "model_byspecies_q.txt"
jagsscript.byspecies_q <- cat("
  model{
  # define prior: probability of q (prob of success for each sire)
  q ~ dbeta(1,1)
  
  # likelihood: stochastic relationship
  # sires ~ (truncated binomial)
  for(i in 1:N){
    nsires[i] ~ dbinom(q, littersize[i]) T(0,)
  }
  
  # DERIVED QUANTITIES
  # Calculate the number of mates (deterministic relation?)
  # q is probability of success here
  # for(i in 1:N){
  #   nmates[i] <- (q * nsires[i]) / (1-q)
  # }
  
  # mates ~ (negative binomial); 
  ##
  # NOTE: in R, representation of NB dist is same as Casella-Berger stat theory book
  # whereas, Ash uses the alternative NB pmf seen in the note for NB dist in Casella-Berger
  # Therefore, P(X=x|r,p) used in R is P(X=M-r|r,q) in terms of number of failures M-r
  # mean of NB is μ = n(1-p)/p and variance n(1-p)/p^2, where p is the prob of SUCCESS
  # R is calculating the number of failures M-r, and we generated the number of successes (sires, R)
  # therefore the number of mates M = failures + successes = failed.m + nsires
  ##
  for(i in 1:N){
    failed.m[i] ~ dnbinom(q, nsires[i]) T(nsires[i],)
    nmates[i] <- failed.m[i] + nsires[i]
  }
  for(i in 1:N){
    est.M[i] <- nsires[i]/q
  }
  
  # Calculated quantities from original data
  for(j in 1:J){
    avg.R[j] <- avgbrood[j]*q/(1-(1-q)^avgbrood[j])     # avg sires, given avgbrood and estimated q
  }
  for(j in 1:J){
    avg.M[j] <- avgsire[j]/q     # avg mates, given avgsire and estimated q
  }
  for(j in 1:J){
    est.pmult.kq[j] <- 1 - avgbrood[j] * q * (1 - q)^(avgbrood[j] - 1) / (1 - (1 - q)^(avgbrood[j]))   
    # estimated prob. of mult. paternity, given avgbrood and estimated q
  }
}
",file = model_byspecies_q.file)



#### BEGINNING OF REGENERATION OF DATA LOOP ####
# Remove any datapoints with NA for avgsire, since Poisson data gen will fail
mates.dat <- mates.df[!is.na(mates.df$avgsire),]

set.seed(36849)

#### (1) generate broods and sires from Poisson for each species ####
mates <- cbind(mates.dat, littersize = NA, nsires = NA)
for(i in 1:nrow(mates)){
  # generate broods
  nsample <- mates[i,]$nbrood   # number of broods
  avg.bs <- mates[i,]$avgbrood  # average brood size
  gen.k <- rtpois(1, avg.bs, a = 0, b = Inf)
  mates[i,]$littersize <- gen.k
  # generate sires
  avg.sire <- mates[i,]$avgsire  # average number of sires
  sire.max <- mates[i,]$littersize # max number of sires per litter is litter size
  gen.r <- rtpois(1, avg.sire, a = 0, b = sire.max)
  mates[i,]$nsires <- gen.r
}
## remove multiple paternity proportion of 1 or 0
#mates <- mates[mates$pmult != 0 & mates$pmult != 1,]
## remove an average sire value of NA
mates <- mates[!is.na(mates$avgsire), ]

## save data generation (b/c it takes too long to redo each time)
# save(mates, file = "mates.rda")
# load("mates.rda")



##### (2) MCMC dimensions #####
ni = 10000
nb = 1000
nc = 2 # needs to match number of initial values proposed
nt = 2
n.iter = ni + nb

##### (3) parameters to monitor #####
params = c("q", "nsires", "failed.m", "nmates", "est.M", "avg.R", "avg.M", "est.pmult.kq")


##### CREATE LIST TO HOLD ALL RESULTS FROM JAGS
jags.results <- list() 

##### (5) run the model in JAGS #####
# Initial conditions for loop
# Likelihood needs to change for each nsire value, so q calculated for each # nsire
mates.p <- mates
species <- sort(unique(mates.p$species))
b <- 1
# Start JAGS
#starttime <- Sys.time()
for(r in species){
  # data containing only relevant values
  m.dat <- mates.p[mates.p$species==r,]
  orig.dat <- dat[dat$species==r,]
  jags.dat <- list(nsires = m.dat$nsires, 
                   littersize = m.dat$littersize, 
                   N = nrow(m.dat), 
                   avgbrood = orig.dat$avgbrood,
                   avgsire = orig.dat$avgsire,
                   J = nrow(orig.dat))
  
  # run JAGS (not using initial values for q)
  jmod <- jags.model(file = model_byspecies_q.file, data = jags.dat, n.chains = nc, n.adapt = 1000)  
  update(jmod, n.iter = nb, by = 1, progress.bar = 'text')
  post <- coda.samples(jmod, params, n.iter = ni, thin = nt) 
  
  # save current JAGS output
  name <- paste0("species=",r)
  jags.results[[name]] <- post
  
  # Save quantities from data
  MCMC_temp <- cbind(
    as.character(orig.dat$species),
    as.numeric(orig.dat$avgbrood),
    as.numeric(orig.dat$avgsire),
    as.numeric(orig.dat$pmult),
    post.summ(post, "q")[[1]],  # mean of the est. param. q
    post.summ(post, "q")[[2]],  # sd of the est. param. q
    post.summ(post, "q")[[4]],  # 2.5% est. param. q
    post.summ(post, "q")[[5]],  # 97.5% est. param. q
    post.summ(post, "nsires")[[1]],  # mean of nsires
    post.summ(post, "nsires")[[2]],  # sd of nsires
    post.summ(post, "nsires")[[4]],  # 2.5% nsires
    post.summ(post, "nsires")[[5]],  # 97.5% nsires
    post.summ(post, "nmates")[[1]],  # mean of nmates
    post.summ(post, "nmates")[[2]],  # sd of nmates
    post.summ(post, "nmates")[[4]],  # 2.5% nmates
    post.summ(post, "nmates")[[5]],  # 97.5% nmates
    post.summ(post, "est.M")[[1]],  # mean of est.M
    post.summ(post, "est.M")[[2]],  # sd of est.M
    post.summ(post, "est.M")[[4]],  # 2.5% est.M
    post.summ(post, "est.M")[[5]],  # 97.5% est.M
    post.summ(post, "est.pmult.kq")[[1]],  # mean of the est.pmult.kq
    post.summ(post, "est.pmult.kq")[[2]],  # sd of the est.pmult.kq
    post.summ(post, "est.pmult.kq")[[4]],  # 2.5% est.pmult.kq
    post.summ(post, "est.pmult.kq")[[5]],  # 97.5% est.pmult.kq
    post.summ(post, "avg.R")[[1]],  # mean of avg.R
    post.summ(post, "avg.R")[[2]],  # sd of avg.R
    post.summ(post, "avg.R")[[4]],  # 2.5% avg.R
    post.summ(post, "avg.R")[[5]],  # 97.5% avg.R
    post.summ(post, "avg.M")[[1]],  # mean of avg.M
    post.summ(post, "avg.M")[[2]],  # sd of avg.M
    post.summ(post, "avg.M")[[4]],  # 2.5% avg.M
    post.summ(post, "avg.M")[[5]]  # 97.5% avg.M
  )
  
  if(b==1) MCMC_sumtemp <- MCMC_temp else MCMC_sumtemp <- rbind(MCMC_sumtemp, MCMC_temp)
  
  # move to next step
  b <- b + 1
}

MCMC_summary <- data.frame(species = as.character(MCMC_sumtemp[,1]),
                           avgbrood = as.numeric(MCMC_sumtemp[,2]),
                           avgsire = as.numeric(MCMC_sumtemp[,3]),
                           pmult = as.numeric(MCMC_sumtemp[,4]),
                           mean.q = as.numeric(MCMC_sumtemp[,5]),
                           sd.q = as.numeric(MCMC_sumtemp[,6]),
                           q2.5 = as.numeric(MCMC_sumtemp[,7]),
                           q97.5 = as.numeric(MCMC_sumtemp[,8]),
                           mean.nsires = as.numeric(MCMC_sumtemp[,9]),
                           sd.nsires = as.numeric(MCMC_sumtemp[,10]),
                           nsires2.5 = as.numeric(MCMC_sumtemp[,11]),
                           nsires97.5 = as.numeric(MCMC_sumtemp[,12]),
                           mean.nmates = as.numeric(MCMC_sumtemp[,13]),
                           sd.nmates = as.numeric(MCMC_sumtemp[,14]),
                           nmates2.5 = as.numeric(MCMC_sumtemp[,15]),
                           nmates97.5 = as.numeric(MCMC_sumtemp[,16]),
                           mean.estM = as.numeric(MCMC_sumtemp[,17]),
                           sd.estM = as.numeric(MCMC_sumtemp[,18]),
                           estM2.5 = as.numeric(MCMC_sumtemp[,19]),
                           estM97.5 = as.numeric(MCMC_sumtemp[,20]),
                           mean.est.pmult = as.numeric(MCMC_sumtemp[,21]),
                           sd.est.pmult = as.numeric(MCMC_sumtemp[,22]),
                           est.pmult2.5 = as.numeric(MCMC_sumtemp[,23]),
                           est.pmult97.5 = as.numeric(MCMC_sumtemp[,24]),
                           mean.avgR = as.numeric(MCMC_sumtemp[,25]),
                           sd.avgR = as.numeric(MCMC_sumtemp[,26]),
                           avgR2.5 = as.numeric(MCMC_sumtemp[,27]),
                           avgR97.5 = as.numeric(MCMC_sumtemp[,28]),
                           mean.avgM = as.numeric(MCMC_sumtemp[,29]),
                           sd.avgM = as.numeric(MCMC_sumtemp[,30]),
                           avgM2.5 = as.numeric(MCMC_sumtemp[,31]),
                           avgM97.5 = as.numeric(MCMC_sumtemp[,32]))

save(MCMC_summary, file = paste0("MCMC_summary_birds.rda"))


# ################ DIAGNOSTICS BELOW ###################
# current <- 1 #change for different species (61 total)
# nam.current <- names(jags.results[current])
# jags.post <- jags.results[[nam.current]]
# gelman.diag(jags.post, multivariate = F) # values below 1.1 should be okay
# gelman.plot(jags.post, multivariate = F)
# n.eff <- effectiveSize(jags.post)
# 
# # visualize trace and posterior plots
# par(mar=c(2,2,1,1))
# plot(jags.post)
# 
# ##### make inference #####
# post.summ(jags.post, "q") # probability of success for each sire for a given data-gen run
# post.summ(jags.post, "nsires")
# post.summ(jags.post, "nmates")
# post.summ(jags.post, "avg.R") # avg number of sires given the est prob of success for each sire, using avg litter size
# post.summ(jags.post, "avg.M")
# post.summ(jags.post, "est.pmult.kq")

# #### Compare the est.sires (calculated with est. q) to true.sires
# sires.results <- sires.results[order(sires.results$species),]
# head(sires.results)
# for(i in unique(sires.results$true.sires)){
#   print(apply(sires.results[sires.results$true.sires==i, 2:5], 2, mean))
# }


################ CALCULATE RESIDUALS ###################
MCMC_resids <- MCMC_summary
loess_pmult <- loess(mean.est.pmult ~ avgbrood, data = MCMC_resids)
MCMC_resids$pred_pmult <- predict(loess_pmult, newdata = MCMC_resids$avgbrood)
# lower limit
loess_pmult_lcl <- loess(est.pmult2.5 ~ avgbrood, data = MCMC_resids)
MCMC_resids$pred_pmult_lcl <- predict(loess_pmult_lcl, newdata = MCMC_resids$avgbrood)
# upper limit
loess_pmult_ucl <- loess(est.pmult97.5 ~ avgbrood, data = MCMC_resids)
MCMC_resids$pred_pmult_ucl <- predict(loess_pmult_ucl, newdata = MCMC_resids$avgbrood)
# residual calculations
MCMC_resids$resid_pmult <- MCMC_resids$pred_pmult - MCMC_resids$pmult
MCMC_resids$resid_pmult_lcl <- MCMC_resids$pred_pmult_lcl - MCMC_resids$pmult
MCMC_resids$resid_pmult_ucl <- MCMC_resids$pred_pmult_ucl - MCMC_resids$pmult

# save dataframe with all residual calculations
save(MCMC_resids, file = "MCMC_resids_birds.rda")
# write.csv(MCMC_resids, file = "MCMC_resids_birds.csv")


################ CALCULATE EFFECT SIZES ###################
birds_1a <- read.csv("paternity_birds.csv", header = TRUE)
## remove birds with nbrood = NA
birds_1b <- birds_1a[!is.na(birds_1a$nbrood),]
birds_1 <- birds_1b[!is.na(birds_1b$avgsire),]

# yi = Bayes pmult - true pmult (i.e. MCMC_total_resids$resid_pmult)
# vi = pmult*(1-pmult)/nbrood

birds_singlerun1 <- merge(birds_1[,c(1:3)], MCMC_resids, by = c("species", "avgbrood"))
# order by brood size
birds_singlerun <- birds_singlerun1[order(birds_singlerun1$avgbrood, birds_singlerun1$pmult, birds_singlerun1$nbrood),]

### variance p(1-p)/n for fixed p assumption
birds_singlerun$vi <- birds_singlerun$pmult*(1-birds_singlerun$pmult)/birds_singlerun$nbrood

## labels
bird_species <- birds_singlerun[order(birds_singlerun$avgbrood),]$species

library(metafor)
paternity.meta.bayes <- rma(resid_pmult, vi, data = birds_singlerun, slab = bird_species)
paternity.meta.bayes


#### Need to plot forest plot over three pages (way too long otherwise)
## Part 1
res1 <- paternity.meta.bayes
res1$vi.f <- res1$vi.f[1:68]
res1$yi.f <- res1$yi.f[1:68]
res1$slab <- bird_species[1:68]


## Species should already be ordered by avgbrood and labelled correctly - see above
pdf(file = paste0("forest_birds1-w-labels.pdf"), width = 9, height = 15)
print(
  forest(res1, #slab = part1.species, 
         #order = order(birds_singlerun[birds_singlerun$species %in% part1.species,]$avgbrood), 
         cex = 1,
         xlim = c(-2.5,2),
         xlab = expression(paste(p[B], " - p"))),
  quote = FALSE,
  text(-2.5, 70, "Species",  pos=4),
  text(2, 70, expression(paste(p[B], " - p [95% CI]")), pos=2)
)
dev.off()

## Part 2
res2 <- paternity.meta.bayes
res2$vi.f <- res2$vi.f[69:136]
res2$yi.f <- res2$yi.f[69:136]
res2$slab <- bird_species[69:136]


pdf(file = paste0("forest_birds2-w-labels.pdf"), width = 9, height = 15)
print(
  forest(res2, #slab = part2.species, 
         #order = order(birds_singlerun[birds_singlerun$species %in% part2.species,]$avgbrood), 
         cex = 1,
         xlim = c(-2.5,2),
         xlab = expression(paste(p[B], " - p"))),
  quote = FALSE,
  text(-2.5, 70, "Species",  pos=4),
  text(2, 70, expression(paste(p[B], " - p [95% CI]")), pos=2)
)
dev.off()





#### **Bayesian analysis of mammals** ####

# clear the workspace
rm(list = ls(all = T))

library(rjags)
library(R2OpenBUGS)
library(coda)
library(extraDistr)
library(ggplot2)

##### specify Ben's post.summ function
post.summ = function(post, var) {
  
  # coerce to matrix for easy subsetting
  post.samp = as.matrix(post)
  
  # if parameter is indexed
  if(substr(var, nchar(var), nchar(var)) == "[") {
    # extract columns with headers equal to the desired variable
    post.sub = post.samp[,substr(colnames(post.samp), 1, nchar(var)) == var]
    # calculate desired quantities
    summ = apply(post.sub, 2, function(x) c(mean = mean(x), sd = sd(x), quantile(x, c(0.5, 0.025, 0.975))))
    return(summ)
  }
  
  # if parameter is not indexed
  if(substr(var, nchar(var), nchar(var)) != "[") {
    # extract the column with the same header as the desired variable
    post.sub = post.samp[,substr(colnames(post.samp), 1, nchar(var)) == var]
    # calculate the desired quantities
    summ = c(mean = mean(post.sub), sd = sd(post.sub), quantile(post.sub, c(0.5, 0.025, 0.975)))
    return(summ)
  }
}

#########################################################
#### ESTIMATION OF Q (RATE OF SUCCESS OF EACH SIRE) 
#### AND CALCULATION OF M (# MATES) FROM R (# SIRES)
#########################################################

##### read and prepare data #####
dat <- read.csv("paternity_mammals.csv")

# expand data to create multiple broods within each species (sample size * 10)
for(i in 1:nrow(dat)){
  temp.sp <- expand.grid(species = rep(dat[i,]$species, dat[i,]$nbrood*10))
  if(i == 1){
    temp <- temp.sp
  } else {
    temp <- rbind(temp, temp.sp)
  } 
}
mates.df <- merge(temp, dat, all.x = TRUE, all.y = FALSE)


##### specify model code #####
model_byspecies_q.file <- "model_byspecies_q.txt"
jagsscript.byspecies_q <- cat("
  model{
  # define prior: probability of q (prob of success for each sire)
  q ~ dbeta(1,1)
  
  # likelihood: stochastic relationship
  # sires ~ (truncated binomial)
  for(i in 1:N){
    nsires[i] ~ dbinom(q, littersize[i]) T(0,)
  }
  
  # DERIVED QUANTITIES
  # Calculate the number of mates (deterministic relation?)
  # q is probability of success here
  # for(i in 1:N){
  #   nmates[i] <- (q * nsires[i]) / (1-q)
  # }
  
  # mates ~ (negative binomial); 
  ##
  # NOTE: in R, representation of NB dist is same as Casella-Berger stat theory book
  # whereas, Ash uses the alternative NB pmf seen in the note for NB dist in Casella-Berger
  # Therefore, P(X=x|r,p) used in R is P(X=M-r|r,q) in terms of number of failures M-r
  # mean of NB is μ = n(1-p)/p and variance n(1-p)/p^2, where p is the prob of SUCCESS
  # R is calculating the number of failures M-r, and we generated the number of successes (sires, R)
  # therefore the number of mates M = failures + successes = failed.m + nsires
  ##
  for(i in 1:N){
    failed.m[i] ~ dnbinom(q, nsires[i]) T(nsires[i],)
    nmates[i] <- failed.m[i] + nsires[i]
  }
  for(i in 1:N){
    est.M[i] <- nsires[i]/q
  }
  
  # Calculated quantities from original data
  for(j in 1:J){
    avg.R[j] <- avgbrood[j]*q/(1-(1-q)^avgbrood[j])     # avg sires, given avgbrood and estimated q
  }
  for(j in 1:J){
    avg.M[j] <- avgsire[j]/q     # avg mates, given avgsire and estimated q
  }
  for(j in 1:J){
    est.pmult.kq[j] <- 1 - avgbrood[j] * q * (1 - q)^(avgbrood[j] - 1) / (1 - (1 - q)^(avgbrood[j]))   
    # estimated prob. of mult. paternity, given avgbrood and estimated q
  }
}
",file = model_byspecies_q.file)



#### BEGINNING OF REGENERATION OF DATA LOOP ####
# Remove any datapoints with NA for avgsire and avgbrood, since Poisson data gen will fail
mates.dat <- mates.df[!is.na(mates.df$avgsire),]

set.seed(36849)

#### (1) generate broods and sires from Poisson for each species ####
mates <- cbind(mates.dat, littersize = NA, nsires = NA)
for(i in 1:nrow(mates)){
  # generate broods
  nsample <- mates[i,]$nbrood   # number of broods
  avg.bs <- mates[i,]$avgbrood  # average brood size
  gen.k <- rtpois(1, avg.bs, a = 0, b = Inf)
  mates[i,]$littersize <- gen.k
  # generate sires
  avg.sire <- mates[i,]$avgsire  # average number of sires
  sire.max <- mates[i,]$littersize # max number of sires per litter is litter size
  gen.r <- rtpois(1, avg.sire, a = 0, b = sire.max)
  mates[i,]$nsires <- gen.r
}
## remove multiple paternity proportion of 1 or 0
#mates <- mates[mates$pmult != 0 & mates$pmult != 1,]
## remove an average sire value of NA
mates <- mates[!is.na(mates$avgsire), ]

##### (2) MCMC dimensions #####
ni = 10000
nb = 1000
nc = 2 # needs to match number of initial values proposed
nt = 2
n.iter = ni + nb

##### (3) parameters to monitor #####
params = c("q", "nsires", "failed.m", "nmates", "est.M", "avg.R", "avg.M", "est.pmult.kq")


##### CREATE LIST TO HOLD ALL RESULTS FROM JAGS
# # max number of iterations through all combinations of sires and litters
# maxlitter.for.minsire <- length(min(mates.dat$nsires):max(mates.dat$littersize))
# maxlitter.for.maxsire <- length(max(mates.dat$nsires):max(mates.dat$littersize))
# nlist <- sum(maxlitter.for.minsire:maxlitter.for.maxsire) 
# # slightly more slots in list than necessary, e.g. sires 5 & 6 only have max of 13 littersize
# # could change with random generation of litter sizes though
jags.results <- list() 

##### (5) run the model in JAGS #####
# Initial conditions for loop
# Likelihood needs to change for each nsire value, so q calculated for each # nsire
mates.p <- mates
species <- sort(unique(mates.p$species))
b <- 1
# Start JAGS
#starttime <- Sys.time()
for(r in species){
  # data containing only relevant values
  m.dat <- mates.p[mates.p$species==r,]
  orig.dat <- dat[dat$species==r,]
  jags.dat <- list(nsires = m.dat$nsires, 
                   littersize = m.dat$littersize, 
                   N = nrow(m.dat), 
                   avgbrood = orig.dat$avgbrood,
                   avgsire = orig.dat$avgsire,
                   J = nrow(orig.dat))
  
  # run JAGS (not using initial values for q)
  jmod <- jags.model(file = model_byspecies_q.file, data = jags.dat, n.chains = nc, n.adapt = 1000)  
  update(jmod, n.iter = nb, by = 1, progress.bar = 'text')
  post <- coda.samples(jmod, params, n.iter = ni, thin = nt) 
  
  # save current JAGS output
  name <- paste0("species=",r)
  jags.results[[name]] <- post
  
  # Save quantities from data
  MCMC_temp <- cbind(
    as.character(orig.dat$species),
    as.numeric(orig.dat$avgbrood),
    as.numeric(orig.dat$avgsire),
    as.numeric(orig.dat$pmult),
    post.summ(post, "q")[[1]],  # mean of the est. param. q
    post.summ(post, "q")[[2]],  # sd of the est. param. q
    post.summ(post, "q")[[4]],  # 2.5% est. param. q
    post.summ(post, "q")[[5]],  # 97.5% est. param. q
    post.summ(post, "nsires")[[1]],  # mean of nsires
    post.summ(post, "nsires")[[2]],  # sd of nsires
    post.summ(post, "nsires")[[4]],  # 2.5% nsires
    post.summ(post, "nsires")[[5]],  # 97.5% nsires
    post.summ(post, "nmates")[[1]],  # mean of nmates
    post.summ(post, "nmates")[[2]],  # sd of nmates
    post.summ(post, "nmates")[[4]],  # 2.5% nmates
    post.summ(post, "nmates")[[5]],  # 97.5% nmates
    post.summ(post, "est.M")[[1]],  # mean of est.M
    post.summ(post, "est.M")[[2]],  # sd of est.M
    post.summ(post, "est.M")[[4]],  # 2.5% est.M
    post.summ(post, "est.M")[[5]],  # 97.5% est.M
    post.summ(post, "est.pmult.kq")[[1]],  # mean of the est.pmult.kq
    post.summ(post, "est.pmult.kq")[[2]],  # sd of the est.pmult.kq
    post.summ(post, "est.pmult.kq")[[4]],  # 2.5% est.pmult.kq
    post.summ(post, "est.pmult.kq")[[5]],  # 97.5% est.pmult.kq
    post.summ(post, "avg.R")[[1]],  # mean of avg.R
    post.summ(post, "avg.R")[[2]],  # sd of avg.R
    post.summ(post, "avg.R")[[4]],  # 2.5% avg.R
    post.summ(post, "avg.R")[[5]],  # 97.5% avg.R
    post.summ(post, "avg.M")[[1]],  # mean of avg.M
    post.summ(post, "avg.M")[[2]],  # sd of avg.M
    post.summ(post, "avg.M")[[4]],  # 2.5% avg.M
    post.summ(post, "avg.M")[[5]]  # 97.5% avg.M
  )
  
  if(b==1) MCMC_sumtemp <- MCMC_temp else MCMC_sumtemp <- rbind(MCMC_sumtemp, MCMC_temp)
  
  # move to next step
  b <- b + 1
}

MCMC_summary <- data.frame(species = as.character(MCMC_sumtemp[,1]),
                           avgbrood = as.numeric(MCMC_sumtemp[,2]),
                           avgsire = as.numeric(MCMC_sumtemp[,3]),
                           pmult = as.numeric(MCMC_sumtemp[,4]),
                           mean.q = as.numeric(MCMC_sumtemp[,5]),
                           sd.q = as.numeric(MCMC_sumtemp[,6]),
                           q2.5 = as.numeric(MCMC_sumtemp[,7]),
                           q97.5 = as.numeric(MCMC_sumtemp[,8]),
                           mean.nsires = as.numeric(MCMC_sumtemp[,9]),
                           sd.nsires = as.numeric(MCMC_sumtemp[,10]),
                           nsires2.5 = as.numeric(MCMC_sumtemp[,11]),
                           nsires97.5 = as.numeric(MCMC_sumtemp[,12]),
                           mean.nmates = as.numeric(MCMC_sumtemp[,13]),
                           sd.nmates = as.numeric(MCMC_sumtemp[,14]),
                           nmates2.5 = as.numeric(MCMC_sumtemp[,15]),
                           nmates97.5 = as.numeric(MCMC_sumtemp[,16]),
                           mean.estM = as.numeric(MCMC_sumtemp[,17]),
                           sd.estM = as.numeric(MCMC_sumtemp[,18]),
                           estM2.5 = as.numeric(MCMC_sumtemp[,19]),
                           estM97.5 = as.numeric(MCMC_sumtemp[,20]),
                           mean.est.pmult = as.numeric(MCMC_sumtemp[,21]),
                           sd.est.pmult = as.numeric(MCMC_sumtemp[,22]),
                           est.pmult2.5 = as.numeric(MCMC_sumtemp[,23]),
                           est.pmult97.5 = as.numeric(MCMC_sumtemp[,24]),
                           mean.avgR = as.numeric(MCMC_sumtemp[,25]),
                           sd.avgR = as.numeric(MCMC_sumtemp[,26]),
                           avgR2.5 = as.numeric(MCMC_sumtemp[,27]),
                           avgR97.5 = as.numeric(MCMC_sumtemp[,28]),
                           mean.avgM = as.numeric(MCMC_sumtemp[,29]),
                           sd.avgM = as.numeric(MCMC_sumtemp[,30]),
                           avgM2.5 = as.numeric(MCMC_sumtemp[,31]),
                           avgM97.5 = as.numeric(MCMC_sumtemp[,32]))

save(MCMC_summary, file = paste0("MCMC_summary_mammals.rda"))


# ################ DIAGNOSTICS BELOW ###################
# current <- 1 #change for different species (61 total)
# nam.current <- names(jags.results[current])
# jags.post <- jags.results[[nam.current]]
# gelman.diag(jags.post, multivariate = F) # values below 1.1 should be okay
# gelman.plot(jags.post, multivariate = F)
# n.eff <- effectiveSize(jags.post)
# 
# # visualize trace and posterior plots
# par(mar=c(2,2,1,1))
# plot(jags.post)
# 
# ##### make inference #####
# post.summ(jags.post, "q") # probability of success for each sire for a given data-gen run
# post.summ(jags.post, "nsires")
# post.summ(jags.post, "nmates")
# post.summ(jags.post, "avg.R") # avg number of sires given the est prob of success for each sire, using avg litter size
# post.summ(jags.post, "avg.M")
# post.summ(jags.post, "est.pmult.kq")

# #### Compare the est.sires (calculated with est. q) to true.sires
# sires.results <- sires.results[order(sires.results$species),]
# head(sires.results)
# for(i in unique(sires.results$true.sires)){
#   print(apply(sires.results[sires.results$true.sires==i, 2:5], 2, mean))
# }


################ CALCULATE RESIDUALS ###################
MCMC_resids <- MCMC_summary
loess_pmult <- loess(mean.est.pmult ~ avgbrood, data = MCMC_resids)
MCMC_resids$pred_pmult <- predict(loess_pmult, newdata = MCMC_resids$avgbrood)
# lower limit
loess_pmult_lcl <- loess(est.pmult2.5 ~ avgbrood, data = MCMC_resids)
MCMC_resids$pred_pmult_lcl <- predict(loess_pmult_lcl, newdata = MCMC_resids$avgbrood)
# upper limit
loess_pmult_ucl <- loess(est.pmult97.5 ~ avgbrood, data = MCMC_resids)
MCMC_resids$pred_pmult_ucl <- predict(loess_pmult_ucl, newdata = MCMC_resids$avgbrood)
# residual calculations
MCMC_resids$resid_pmult <- MCMC_resids$pred_pmult - MCMC_resids$pmult
MCMC_resids$resid_pmult_lcl <- MCMC_resids$pred_pmult_lcl - MCMC_resids$pmult
MCMC_resids$resid_pmult_ucl <- MCMC_resids$pred_pmult_ucl - MCMC_resids$pmult

# save dataframe with all residual calculations
save(MCMC_resids, file = "MCMC_resids_mammals.rda")
# write.csv(MCMC_resids, file = "MCMC_resids_mammals.csv")



################ CALCULATE EFFECT SIZES ###################
mammals_1 <- read.csv("paternity_mammals.csv", header = TRUE)
mammals_1 <- mammals_1[!is.na(mammals_1$avgsire),]

# yi = Bayes pmult - true pmult (i.e. MCMC_total_resids$resid_pmult)
# vi = pmult*(1-pmult)/nbrood

mammals_singlerun <- merge(mammals_1[,c(1:3)], MCMC_resids, by = c("species", "avgbrood"))
# order by brood size
mammals_singlerun <- mammals_singlerun[order(mammals_singlerun$avgbrood, mammals_singlerun$pmult, mammals_singlerun$nbrood),]

### variance p(1-p)/n for fixed p assumption
mammals_singlerun$vi <- mammals_singlerun$pmult*(1-mammals_singlerun$pmult)/mammals_singlerun$nbrood

library(metafor)
paternity.meta.bayes <- rma(resid_pmult, vi, data = mammals_singlerun)
paternity.meta.bayes
pdf(file = paste0("forest_mammals-w-labels.pdf"), width = 11.5, height = 13)
print(
  forest(paternity.meta.bayes, slab = mammals_singlerun$species, 
         order = order(mammals_singlerun$avgbrood), cex = 1,
         xlim = c(-2.5,2),
         xlab = expression(paste(p[B], " - p"))),
  quote = FALSE,
  text(-2.5, 65, "Species",  pos=4),
  text(2, 65, expression(paste(p[B], " - p [95% CI]")), pos=2)
)
dev.off()

