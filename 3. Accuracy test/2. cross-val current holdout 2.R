### R code for data input for Cross validation test (create two datasets; one half = predict, one half = holdout)
### This script is to obtain predictive results from cross-validation
### Script from Ardiantiono et al. 2024. Improved occupancy and cost-effectiveness of species monitoring programs through data integration
### Jags adaptation from Nimble script by Jeffrey W. Doser (doserjef@msu.edu)
#==================================================================================================================


rm(list = ls())
library(coda)
library(tidyverse)

# Set working directory (where the all data required to run the models are stored)
# Please amend the filepath to ensure compatibility with your folder structure
setwd("D:\\Data analysis\\Chapter 2_Integrated\\Ardiantiono et al. 2024. Data integration\\3. Accuracy test")
getwd()

# Load in raw data (USE INPUT CROSSVAL 2) --------------------------------------------------------
source("D:\\Data analysis\\Chapter 2_Integrated\\Ardiantiono et al. 2024. Data integration\\3. Accuracy test\\0. Input file integrated_Global_CrossVal_2_holdout.R") 
load("0. cross-val-indices.R")

# Functions ---------------------------------------------------------------
logit <- function(theta, a = 0, b = 1) {log((theta-a)/(b-theta))}
logit.inv <- function(z, a = 0, b = 1) {b-(b-a)/(1+exp(z))}

# CALL THIS CODE
# These functions compute log pointwise-posterior predictive densities that...
# are used to assess predictive performance. 
source("2. lpd-sampler.R")

# Current hold-out set ----------------------------------------------------
# Takes value 1 or 2 depending on which set is the current hold out set
# Provided as user input when running script from terminal.
curr.set <- 2
pred.set <- ifelse(curr.set == 1, 2, 1)


# Read in all data sets ---------------------------------------------------
# SMART --------------------------------
samples.sm <- list()
jags.out <- load("Results\\icom-SM-results-4e+05-iterations-1-holdout-.R")
out.list <- as.matrix(as.data.frame(out[[1]]))
out.list <- out.list[,-c(1:7)] #we remove community parameters to extract other parameter easily
samples.sm[[1]] <- out.list
jags.out <- load("Results\\icom-SM-results-4e+05-iterations-2-holdout-.R")
out.list <- as.matrix(as.data.frame(out[[1]]))
out.list <- out.list[,-c(1:7)]
samples.sm[[2]] <- out.list

# CT-TR-SM (Integrated) --------------------------------
samples.ct.tr.sm <- list()
jags.out <- load("Results\\icom-CT-TR-SM-results-4e+05-iterations-1-holdout-.R")
out.list <- as.matrix(as.data.frame(out[[1]]))
out.list <- out.list[,-c(1:7)] #we remove community parameters to extract other par easily
samples.ct.tr.sm[[1]] <- out.list
jags.out <- load("Results\\icom-CT-TR-SM-results-4e+05-iterations-2-holdout-.R")
out.list <- as.matrix(as.data.frame(out[[1]]))
out.list <- out.list[,-c(1:7)]
samples.ct.tr.sm[[2]] <- out.list


# Define constants --------------------------------------------------------
I <- 3
sp <- c('SP1', 'SP2', 'SP3')
sp.names <- c('SP1', 'SP2', 'SP3')
sp.indx <- 1:3
# Initiate arrays to store lpd values. Adjust the structure based on combination models
# For example if we have seven models as predictors and we want to test how accurate they predict
# species occurence in three datasets (i.e. sites with patrol, transect, camera data)
# then we set it as array(NA, dim = c(I, 3, 7))
elpd.ct.vals <- array(NA, dim = c(I, 2, 2))
elpd.tr.vals <- array(NA, dim = c(I, 2, 2))
elpd.sm.vals <- array(NA, dim = c(I, 2, 2))


# The remainder of the code computes the lpd according to the description 
# provided in Doser et al. (2021). A single log-pointwide predictive density
# measure is calculated for each model, for each species, for each data set location included
# in that model, and for each model that could have potentially generated
# those data. This approach accounts for model uncertainty and enables 
# assessment of predictive performance for individual species, the entire 
# community, and at different parts of the entire region of interest. 

###########################
# Patrol as predicting model -------------------------------------------------

##### Patrol dataset 
samples.fit <- samples.sm[[curr.set]]
samples.pred <- samples.sm[[pred.set]]

access.fit <- c(covars.sm$Access[-sm.x.indices[[curr.set]]])
mean.access.fit <- mean(access.fit)
sd.access.fit <- sd(access.fit)
access.fit <- (access.fit - mean(access.fit)) / sd(access.fit)
access.pred <- c(covars.sm$Access[sm.x.indices[[curr.set]]])
access.pred <- (access.pred - mean.access.fit) / sd.access.fit 

tree.fit <- c(covars.sm$Tree[-sm.x.indices[[curr.set]]])
mean.tree.fit <- mean(tree.fit)
sd.tree.fit <- sd(tree.fit)
tree.fit <- (tree.fit - mean(tree.fit)) / sd(tree.fit)
tree.pred <- c(covars.sm$Tree[sm.x.indices[[curr.set]]])
tree.pred <- (tree.pred - mean.tree.fit) / sd.tree.fit 

canopy.fit <- c(covars.sm$Tree_var[-sm.x.indices[[curr.set]]])
mean.canopy.fit <- mean(canopy.fit)
sd.canopy.fit <- sd(canopy.fit)
canopy.fit <- (canopy.fit - mean(canopy.fit)) / sd(canopy.fit)
canopy.pred <- c(covars.sm$Tree_var[sm.x.indices[[curr.set]]])
canopy.pred <- (canopy.pred - mean.canopy.fit) / sd.canopy.fit 

n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.sm <- length(sm.x.indices[[curr.set]])
elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                    samples.pred = samples.pred, 
                    access.pred = access.pred, 
                    tree.pred = tree.pred, 
                    canopy.pred = canopy.pred,
                    my.iter = my.iter, 
                    n.iter = n.iter, 
                    I = I, type = 'sm')
elpd.vals.1a <- elpd.vals.1[, !colSums(is.na(elpd.vals.1))] #OPTIONAL if there's columns with NAs
elpd.sm.vals[, 1, 1] <- apply(log(elpd.vals.1a), 1, sum)


#### CT + TR + SM (Integrated)
samples.fit <- samples.sm[[curr.set]]
samples.pred <- samples.ct.tr.sm[[pred.set]]

access.fit <- c(covars.sm$Access[-sm.x.indices[[curr.set]]])
mean.access.fit <- mean(access.fit)
sd.access.fit <- sd(access.fit)
access.fit <- (access.fit - mean(access.fit)) / sd(access.fit)
access.pred <- c(covars.ct$Access[ct.x.indices[[curr.set]]],
                 covars.tr$Access[tr.x.indices[[curr.set]]],
                 covars.sm$Access[sm.x.indices[[curr.set]]])
access.pred <- (access.pred - mean.access.fit) / sd.access.fit 

tree.fit <- c(covars.sm$Tree[-sm.x.indices[[curr.set]]])
mean.tree.fit <- mean(tree.fit)
sd.tree.fit <- sd(tree.fit)
tree.fit <- (tree.fit - mean(tree.fit)) / sd(tree.fit)
tree.pred <- c(covars.ct$Tree[ct.x.indices[[curr.set]]],
               covars.tr$Tree[tr.x.indices[[curr.set]]],
               covars.sm$Tree[sm.x.indices[[curr.set]]])
tree.pred <- (tree.pred - mean.tree.fit) / sd.tree.fit 

canopy.fit <- c(covars.sm$Tree_var[-sm.x.indices[[curr.set]]])
mean.canopy.fit <- mean(canopy.fit)
sd.canopy.fit <- sd(canopy.fit)
canopy.fit <- (canopy.fit - mean(canopy.fit)) / sd(canopy.fit)
canopy.pred <- c(covars.ct$Tree_var[ct.x.indices[[curr.set]]],
                 covars.tr$Tree_var[tr.x.indices[[curr.set]]],
                 covars.sm$Tree_var[sm.x.indices[[curr.set]]])
canopy.pred <- (canopy.pred - mean.canopy.fit) / sd.canopy.fit 

n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.ct <- length(ct.x.indices[[curr.set]])
J.tr <- length(tr.x.indices[[curr.set]])
J.sm <- length(sm.x.indices[[curr.set]])
access.ct.pred <- access.pred[1:J.ct]
access.tr.pred <- access.pred[1:J.tr + J.ct]
access.sm.pred <- access.pred[1:J.sm + J.tr +J.ct]
tree.ct.pred <- tree.pred[1:J.ct]
tree.tr.pred <- tree.pred[1:J.tr + J.ct]
tree.sm.pred <- tree.pred[1:J.sm + J.tr +J.ct]
canopy.ct.pred <- canopy.pred[1:J.ct]
canopy.tr.pred <- canopy.pred[1:J.tr + J.ct]
canopy.sm.pred <- canopy.pred[1:J.sm + J.tr +J.ct]

#This is to test the performance of patrol data in predicting species occurrence in camera sites
ct.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                       samples.pred = samples.pred, 
                       access.pred = access.ct.pred, 
                       tree.pred = tree.ct.pred, 
                       canopy.pred = canopy.ct.pred,
                       my.iter = my.iter, 
                       n.iter = n.iter, 
                       I = I, type = 'ct')
ct.elpd.vals.1a <- ct.elpd.vals.1[, !colSums(is.na(ct.elpd.vals.1))] #OPTIONAL if there's columns with NAs
elpd.ct.vals[, 1, 1] <- apply(log(ct.elpd.vals.1a), 1, sum)

#This is to test the performance of patrol data in predicting species occurrence in transect sites
tr.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                       samples.pred = samples.pred, 
                       access.pred = access.tr.pred, 
                       tree.pred = tree.tr.pred, 
                       canopy.pred = canopy.tr.pred,
                       my.iter = my.iter, 
                       n.iter = n.iter, 
                       I = I, type = 'tr')
tr.elpd.vals.1a <- tr.elpd.vals.1[, !colSums(is.na(tr.elpd.vals.1))] #OPTIONAL if there's columns with NAs
elpd.tr.vals[, 1, 1] <- apply(log(tr.elpd.vals.1a), 1, sum)

#This is to test the performance of patrol data in predicting species occurrence in patrol sites
sm.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                       samples.pred = samples.pred, 
                       access.pred = access.sm.pred, 
                       tree.pred = tree.sm.pred, 
                       canopy.pred = canopy.sm.pred,
                       my.iter = my.iter, 
                       n.iter = n.iter, 
                       I = I, type = 'sm')
sm.elpd.vals.1a <- sm.elpd.vals.1[, !colSums(is.na(sm.elpd.vals.1))] #OPTIONAL if there's columns with NAs
elpd.sm.vals[, 2, 1] <- apply(log(sm.elpd.vals.1a), 1, sum)

################################################
#################################
#############
# Integrated) as predicting model -------------------------------------------------

##### Patrol 
samples.fit <- samples.sm[[curr.set]]
samples.pred <- samples.sm[[pred.set]]

access.fit <- c(covars.ct$Access[-ct.x.indices[[curr.set]]],
                covars.tr$Access[-tr.x.indices[[curr.set]]],
                covars.sm$Access[-sm.x.indices[[curr.set]]])
mean.access.fit <- mean(access.fit)
sd.access.fit <- sd(access.fit)
access.fit <- (access.fit - mean(access.fit)) / sd(access.fit)
access.pred <- c(covars.sm$Access[sm.x.indices[[curr.set]]])
access.pred <- (access.pred - mean.access.fit) / sd.access.fit 

tree.fit <- c(covars.ct$Tree[-ct.x.indices[[curr.set]]],
              covars.tr$Tree[-tr.x.indices[[curr.set]]],
              covars.sm$Tree[-sm.x.indices[[curr.set]]])
mean.tree.fit <- mean(tree.fit)
sd.tree.fit <- sd(tree.fit)
tree.fit <- (tree.fit - mean(tree.fit)) / sd(tree.fit)
tree.pred <- c(covars.sm$Tree[sm.x.indices[[curr.set]]])
tree.pred <- (tree.pred - mean.tree.fit) / sd.tree.fit 

canopy.fit <- c(covars.ct$Tree_var[-ct.x.indices[[curr.set]]],
                covars.tr$Tree_var[-tr.x.indices[[curr.set]]],
                covars.sm$Tree_var[-sm.x.indices[[curr.set]]])
mean.canopy.fit <- mean(canopy.fit)
sd.canopy.fit <- sd(canopy.fit)
canopy.fit <- (canopy.fit - mean(canopy.fit)) / sd(canopy.fit)
canopy.pred <- c(covars.sm$Tree_var[sm.x.indices[[curr.set]]])
canopy.pred <- (canopy.pred - mean.canopy.fit) / sd.canopy.fit 

n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.sm <- length(sm.x.indices[[curr.set]])
elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                    samples.pred = samples.pred, 
                    access.pred = access.pred, 
                    tree.pred = tree.pred, 
                    canopy.pred = canopy.pred,
                    my.iter = my.iter, 
                    n.iter = n.iter, 
                    I = I, type = 'sm')
elpd.vals.1a <- elpd.vals.1[, !colSums(is.na(elpd.vals.1))] #OPTIONAL if there's columns with NAs
elpd.sm.vals[, 1, 2] <- apply(log(elpd.vals.1a), 1, sum)


#### CT + TR + SM 
samples.fit <- samples.sm[[curr.set]]
samples.pred <- samples.ct.tr.sm[[pred.set]]

access.fit <- c(covars.ct$Access[-ct.x.indices[[curr.set]]],
                covars.tr$Access[-tr.x.indices[[curr.set]]],
                covars.sm$Access[-sm.x.indices[[curr.set]]])
mean.access.fit <- mean(access.fit)
sd.access.fit <- sd(access.fit)
access.fit <- (access.fit - mean(access.fit)) / sd(access.fit)
access.pred <- c(covars.ct$Access[ct.x.indices[[curr.set]]],
                 covars.tr$Access[tr.x.indices[[curr.set]]],
                 covars.sm$Access[sm.x.indices[[curr.set]]])
access.pred <- (access.pred - mean.access.fit) / sd.access.fit 

tree.fit <- c(covars.ct$Tree[-ct.x.indices[[curr.set]]],
              covars.tr$Tree[-tr.x.indices[[curr.set]]],
              covars.sm$Tree[-sm.x.indices[[curr.set]]])
mean.tree.fit <- mean(tree.fit)
sd.tree.fit <- sd(tree.fit)
tree.fit <- (tree.fit - mean(tree.fit)) / sd(tree.fit)
tree.pred <- c(covars.ct$Tree[ct.x.indices[[curr.set]]],
               covars.tr$Tree[tr.x.indices[[curr.set]]],
               covars.sm$Tree[sm.x.indices[[curr.set]]])
tree.pred <- (tree.pred - mean.tree.fit) / sd.tree.fit 

canopy.fit <- c(covars.ct$Tree_var[-ct.x.indices[[curr.set]]],
                covars.tr$Tree_var[-tr.x.indices[[curr.set]]],
                covars.sm$Tree_var[-sm.x.indices[[curr.set]]])
mean.canopy.fit <- mean(canopy.fit)
sd.canopy.fit <- sd(canopy.fit)
canopy.fit <- (canopy.fit - mean(canopy.fit)) / sd(canopy.fit)
canopy.pred <- c(covars.ct$Tree_var[ct.x.indices[[curr.set]]],
                 covars.tr$Tree_var[tr.x.indices[[curr.set]]],
                 covars.sm$Tree_var[sm.x.indices[[curr.set]]])
canopy.pred <- (canopy.pred - mean.canopy.fit) / sd.canopy.fit 

n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.ct <- length(ct.x.indices[[curr.set]])
J.tr <- length(tr.x.indices[[curr.set]])
J.sm <- length(sm.x.indices[[curr.set]])
access.ct.pred <- access.pred[1:J.ct]
access.tr.pred <- access.pred[1:J.tr + J.ct]
access.sm.pred <- access.pred[1:J.sm + J.tr +J.ct]
tree.ct.pred <- tree.pred[1:J.ct]
tree.tr.pred <- tree.pred[1:J.tr + J.ct]
tree.sm.pred <- tree.pred[1:J.sm + J.tr +J.ct]
canopy.ct.pred <- canopy.pred[1:J.ct]
canopy.tr.pred <- canopy.pred[1:J.tr + J.ct]
canopy.sm.pred <- canopy.pred[1:J.sm + J.tr +J.ct]

ct.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                       samples.pred = samples.pred, 
                       access.pred = access.ct.pred, 
                       tree.pred = tree.ct.pred, 
                       canopy.pred = canopy.ct.pred,
                       my.iter = my.iter, 
                       n.iter = n.iter, 
                       I = I, type = 'ct')
ct.elpd.vals.1a <- ct.elpd.vals.1[, !colSums(is.na(ct.elpd.vals.1))] #OPTIONAL if there's columns with NAs
elpd.ct.vals[, 1, 2] <- apply(log(ct.elpd.vals.1a), 1, sum)

tr.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                       samples.pred = samples.pred, 
                       access.pred = access.tr.pred, 
                       tree.pred = tree.tr.pred, 
                       canopy.pred = canopy.tr.pred,
                       my.iter = my.iter, 
                       n.iter = n.iter, 
                       I = I, type = 'tr')
tr.elpd.vals.1a <- tr.elpd.vals.1[, !colSums(is.na(tr.elpd.vals.1))] #OPTIONAL if there's columns with NAs
elpd.tr.vals[, 1, 2] <- apply(log(tr.elpd.vals.1a), 1, sum)

sm.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                       samples.pred = samples.pred, 
                       access.pred = access.sm.pred, 
                       tree.pred = tree.sm.pred, 
                       canopy.pred = canopy.sm.pred,
                       my.iter = my.iter, 
                       n.iter = n.iter, 
                       I = I, type = 'sm')
sm.elpd.vals.1a <- sm.elpd.vals.1[, !colSums(is.na(sm.elpd.vals.1))] #OPTIONAL if there's columns with NAs
elpd.sm.vals[, 2, 2] <- apply(log(sm.elpd.vals.1a), 1, sum)

## Save all!
save(elpd.ct.vals, elpd.tr.vals, 
     elpd.sm.vals, file = paste('cross-val-results-', 
                                curr.set, '.R', sep = ''))
