# cross-val.R: code to obtain predictive results from cross-validation
#              for the Leuser tiger-prey study
# Adapted from Doser et al. 2021
# Citation: 

rm(list = ls())
library(coda)
library(tidyverse)

# Set working directory
setwd("C:\\Users\\User\\OneDrive - University of Kent\\Data Analysis\\Chapter 2_Integrated modeling\\8. Accuracy test")
getwd()

# Functions ---------------------------------------------------------------
logit <- function(theta, a = 0, b = 1) {log((theta-a)/(b-theta))}
logit.inv <- function(z, a = 0, b = 1) {b-(b-a)/(1+exp(z))}
# These functions compute log pointwise-posterior predictive densities that
# are used to assess predictive performance. 
source("lpd-sampler V1.R")

# Current hold-out set ----------------------------------------------------
# Takes value 1 or 2 depending on which set is the current hold out set
# Provided as user input when running script from terminal. 
# Select 2, as this is our hold-out set
curr.set <- 1 #could use whichever half, doesnt matter bcs we will run the other half too
pred.set <- ifelse(curr.set == 1, 2, 1)

# Read in all data sets ---------------------------------------------------
# CT --------------------------------
samples.ct <- list()
jags.out <- load("icom-CT-results-4e+05-iterations-1-holdout-.R")
out.list <- as.matrix(as.data.frame(out[[1]]))
out.list <- out.list[,-c(1:7)] #we remove community parameters to extract other par easily
samples.ct[[1]] <- out.list
jags.out <- load("icom-CT-results-4e+05-iterations-2-holdout-.R")
out.list <- as.matrix(as.data.frame(out[[1]]))
out.list <- out.list[,-c(1:7)]
samples.ct[[2]] <- out.list
#you can see z.sm has no dimension info (z.sm.1 - z.sm.975, not like z.sm.1[1] etc)
#this is fine as we know we have three species so just divide it

# TR --------------------------------
samples.tr <- list()
jags.out <- load("icom-TR-results-4e+05-iterations-1-holdout-.R")
out.list <- as.matrix(as.data.frame(out[[1]]))
out.list <- out.list[,-c(1:7)] #we remove community parameters to extract other par easily
samples.tr[[1]] <- out.list
jags.out <- load("icom-TR-results-4e+05-iterations-2-holdout-.R")
out.list <- as.matrix(as.data.frame(out[[1]]))
out.list <- out.list[,-c(1:7)]
samples.tr[[2]] <- out.list

# SMART --------------------------------
samples.sm <- list()
jags.out <- load("icom-SM-results-4e+05-iterations-1-holdout-.R")
out.list <- as.matrix(as.data.frame(out[[1]]))
out.list <- out.list[,-c(1:7)] #we remove community parameters to extract other par easily
samples.sm[[1]] <- out.list
jags.out <- load("icom-SM-results-4e+05-iterations-2-holdout-.R")
out.list <- as.matrix(as.data.frame(out[[1]]))
out.list <- out.list[,-c(1:7)]
samples.sm[[2]] <- out.list

# CT-TR --------------------------------
samples.ct.tr <- list()
jags.out <- load("icom-CT-TR-results-4e+05-iterations-1-holdout-.R")
out.list <- as.matrix(as.data.frame(out[[1]]))
out.list <- out.list[,-c(1:7)] #we remove community parameters to extract other par easily
samples.ct.tr[[1]] <- out.list
jags.out <- load("icom-CT-TR-results-4e+05-iterations-2-holdout-.R")
out.list <- as.matrix(as.data.frame(out[[1]]))
out.list <- out.list[,-c(1:7)]
samples.ct.tr[[2]] <- out.list

# CT-SM --------------------------------
samples.ct.sm <- list()
jags.out <- load("icom-CT-SM-results-4e+05-iterations-1-holdout-.R")
out.list <- as.matrix(as.data.frame(out[[1]]))
out.list <- out.list[,-c(1:7)] #we remove community parameters to extract other par easily
samples.ct.sm[[1]] <- out.list
jags.out <- load("icom-CT-SM-results-4e+05-iterations-2-holdout-.R")
out.list <- as.matrix(as.data.frame(out[[1]]))
out.list <- out.list[,-c(1:7)]
samples.ct.sm[[2]] <- out.list

# TR-SM --------------------------------
samples.tr.sm <- list()
jags.out <- load("icom-TR-SM-results-4e+05-iterations-1-holdout-.R")
out.list <- as.matrix(as.data.frame(out[[1]]))
out.list <- out.list[,-c(1:7)] #we remove community parameters to extract other par easily
samples.tr.sm[[1]] <- out.list
jags.out <- load("icom-TR-SM-results-4e+05-iterations-2-holdout-.R")
out.list <- as.matrix(as.data.frame(out[[1]]))
out.list <- out.list[,-c(1:7)]
samples.tr.sm[[2]] <- out.list

# CT-TR-SM --------------------------------
samples.ct.tr.sm <- list()
jags.out <- load("icom-CT-TR-SM-results-4e+05-iterations-1-holdout-.R")
out.list <- as.matrix(as.data.frame(out[[1]]))
out.list <- out.list[,-c(1:7)] #we remove community parameters to extract other par easily
samples.ct.tr.sm[[1]] <- out.list
jags.out <- load("icom-CT-TR-SM-results-4e+05-iterations-2-holdout-.R")
out.list <- as.matrix(as.data.frame(out[[1]]))
out.list <- out.list[,-c(1:7)]
samples.ct.tr.sm[[2]] <- out.list


# Define constants --------------------------------------------------------
I <- 3
sp <- c('PAT', 'RUU', 'SUS')
sp.names <- c('Sumatran tiger', 'Sambar deer', 'Wild pig')
sp.indx <- 1:3
# Initiate arrays to store lpd values. 
elpd.ct.vals <- array(NA, dim = c(I, 4, 7))
elpd.tr.vals <- array(NA, dim = c(I, 4, 7))
elpd.sm.vals <- array(NA, dim = c(I, 4, 7))

# Load in raw data --------------------------------------------------------
source("C:\\Users\\User\\OneDrive - University of Kent\\Data Analysis\\Chapter 2_Integrated modeling\\8. Accuracy test\\Input file integrated_Global_CrossVal_1_holdout.R") 
load("cross-val-indices.R")


# The remainder of the code computes the lpd according to the description 
# provided in Doser et al. (2021). A single log-pointwide predictive density
# measure is calculated for each model, for each species, for each data set location included
# in that model, and for each model that could have potentially generated
# those data. This approach accounts for model uncertainty and enables 
# assessment of predictive performance for individual species, the entire 
# community, and at different parts of the entire region of interest. 

# CT as predicting model -------------------------------------------------
##### CT 
samples.fit <- samples.ct[[curr.set]]
samples.pred <- samples.ct[[pred.set]]

access.fit <- c(covars.ct$Access[-ct.x.indices[[curr.set]]])
mean.access.fit <- mean(access.fit)
sd.access.fit <- sd(access.fit)
access.fit <- (access.fit - mean(access.fit)) / sd(access.fit)
access.pred <- c(covars.ct$Access[ct.x.indices[[curr.set]]])
access.pred <- (access.pred - mean.access.fit) / sd.access.fit 

tree.fit <- c(covars.ct$Tree[-ct.x.indices[[curr.set]]])
mean.tree.fit <- mean(tree.fit)
sd.tree.fit <- sd(tree.fit)
tree.fit <- (tree.fit - mean(tree.fit)) / sd(tree.fit)
tree.pred <- c(covars.ct$Tree[ct.x.indices[[curr.set]]])
tree.pred <- (tree.pred - mean.tree.fit) / sd.tree.fit 

canopy.fit <- c(covars.ct$Canopy[-ct.x.indices[[curr.set]]])
mean.canopy.fit <- mean(canopy.fit)
sd.canopy.fit <- sd(canopy.fit)
canopy.fit <- (canopy.fit - mean(canopy.fit)) / sd(canopy.fit)
canopy.pred <- c(covars.ct$Canopy[ct.x.indices[[curr.set]]])
canopy.pred <- (canopy.pred - mean.canopy.fit) / sd.canopy.fit 

n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.ct <- length(ct.x.indices[[curr.set]])
elpd.vals.1 <- elpd(samples.fit = samples.fit, 
		    samples.pred = samples.pred, 
		    access.pred = access.pred, 
		    tree.pred = tree.pred, 
		    canopy.pred = canopy.pred,
		    my.iter = my.iter, 
		    n.iter = n.iter, 
		    I = I, type = 'ct')
elpd.ct.vals[, 1, 1] <- apply(log(elpd.vals.1), 1, sum)

##### TR 
samples.fit <- samples.ct[[curr.set]]
samples.pred <- samples.tr[[pred.set]]

access.fit <- c(covars.ct$Access[-ct.x.indices[[curr.set]]])
mean.access.fit <- mean(access.fit)
sd.access.fit <- sd(access.fit)
access.fit <- (access.fit - mean(access.fit)) / sd(access.fit)
access.pred <- c(covars.tr$Access[tr.x.indices[[curr.set]]])
access.pred <- (access.pred - mean.access.fit) / sd.access.fit 

tree.fit <- c(covars.ct$Tree[-ct.x.indices[[curr.set]]])
mean.tree.fit <- mean(tree.fit)
sd.tree.fit <- sd(tree.fit)
tree.fit <- (tree.fit - mean(tree.fit)) / sd(tree.fit)
tree.pred <- c(covars.tr$Tree[tr.x.indices[[curr.set]]])
tree.pred <- (tree.pred - mean.tree.fit) / sd.tree.fit 

canopy.fit <- c(covars.ct$Canopy[-ct.x.indices[[curr.set]]])
mean.canopy.fit <- mean(canopy.fit)
sd.canopy.fit <- sd(canopy.fit)
canopy.fit <- (canopy.fit - mean(canopy.fit)) / sd(canopy.fit)
canopy.pred <- c(covars.tr$Canopy[tr.x.indices[[curr.set]]])
canopy.pred <- (canopy.pred - mean.canopy.fit) / sd.canopy.fit 

n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.tr <- length(tr.x.indices[[curr.set]])
elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                    samples.pred = samples.pred, 
                    access.pred = access.pred, 
                    tree.pred = tree.pred, 
                    canopy.pred = canopy.pred,
                    my.iter = my.iter, 
                    n.iter = n.iter, 
                    I = I, type = 'tr')
elpd.tr.vals[, 1, 1] <- apply(log(elpd.vals.1), 1, sum)

##### SM 
samples.fit <- samples.ct[[curr.set]]
samples.pred <- samples.sm[[pred.set]]

access.fit <- c(covars.ct$Access[-ct.x.indices[[curr.set]]])
mean.access.fit <- mean(access.fit)
sd.access.fit <- sd(access.fit)
access.fit <- (access.fit - mean(access.fit)) / sd(access.fit)
access.pred <- c(covars.sm$Access[sm.x.indices[[curr.set]]])
access.pred <- (access.pred - mean.access.fit) / sd.access.fit 

tree.fit <- c(covars.ct$Tree[-ct.x.indices[[curr.set]]])
mean.tree.fit <- mean(tree.fit)
sd.tree.fit <- sd(tree.fit)
tree.fit <- (tree.fit - mean(tree.fit)) / sd(tree.fit)
tree.pred <- c(covars.sm$Tree[sm.x.indices[[curr.set]]])
tree.pred <- (tree.pred - mean.tree.fit) / sd.tree.fit 

canopy.fit <- c(covars.ct$Canopy[-ct.x.indices[[curr.set]]])
mean.canopy.fit <- mean(canopy.fit)
sd.canopy.fit <- sd(canopy.fit)
canopy.fit <- (canopy.fit - mean(canopy.fit)) / sd(canopy.fit)
canopy.pred <- c(covars.sm$Canopy[sm.x.indices[[curr.set]]])
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
elpd.sm.vals[, 1, 1] <- apply(log(elpd.vals.1), 1, sum)

##### CT + TR 
samples.fit <- samples.ct[[curr.set]]
samples.pred <- samples.ct.tr[[pred.set]]

access.fit <- c(covars.ct$Access[-ct.x.indices[[curr.set]]])
mean.access.fit <- mean(access.fit)
sd.access.fit <- sd(access.fit)
access.fit <- (access.fit - mean(access.fit)) / sd(access.fit)
access.pred <- c(covars.ct$Access[ct.x.indices[[curr.set]]],
                 covars.tr$Access[tr.x.indices[[curr.set]]])
access.pred <- (access.pred - mean.access.fit) / sd.access.fit 

tree.fit <- c(covars.ct$Tree[-ct.x.indices[[curr.set]]])
mean.tree.fit <- mean(tree.fit)
sd.tree.fit <- sd(tree.fit)
tree.fit <- (tree.fit - mean(tree.fit)) / sd(tree.fit)
tree.pred <- c(covars.ct$Tree[ct.x.indices[[curr.set]]],
               covars.tr$Tree[tr.x.indices[[curr.set]]])
tree.pred <- (tree.pred - mean.tree.fit) / sd.tree.fit 

canopy.fit <- c(covars.ct$Canopy[-ct.x.indices[[curr.set]]])
mean.canopy.fit <- mean(canopy.fit)
sd.canopy.fit <- sd(canopy.fit)
canopy.fit <- (canopy.fit - mean(canopy.fit)) / sd(canopy.fit)
canopy.pred <- c(covars.ct$Canopy[ct.x.indices[[curr.set]]],
                 covars.tr$Canopy[tr.x.indices[[curr.set]]])
canopy.pred <- (canopy.pred - mean.canopy.fit) / sd.canopy.fit 

n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.ct <- length(ct.x.indices[[curr.set]])
J.tr <- length(tr.x.indices[[curr.set]])
access.ct.pred <- access.pred[1:J.ct]
access.tr.pred <- access.pred[1:J.tr + J.ct]
tree.ct.pred <- tree.pred[1:J.ct]
tree.tr.pred <- tree.pred[1:J.tr + J.ct]
canopy.ct.pred <- canopy.pred[1:J.ct]
canopy.tr.pred <- canopy.pred[1:J.tr + J.ct]

ct.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                    samples.pred = samples.pred, 
                    access.pred = access.ct.pred, 
                    tree.pred = tree.ct.pred, 
                    canopy.pred = canopy.ct.pred,
                    my.iter = my.iter, 
                    n.iter = n.iter, 
                    I = I, type = 'ct')
elpd.ct.vals[, 2, 1] <- apply(log(ct.elpd.vals.1), 1, sum)

tr.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                       samples.pred = samples.pred, 
                       access.pred = access.tr.pred, 
                       tree.pred = tree.tr.pred, 
                       canopy.pred = canopy.tr.pred,
                       my.iter = my.iter, 
                       n.iter = n.iter, 
                       I = I, type = 'tr')
elpd.tr.vals[, 2, 1] <- apply(log(tr.elpd.vals.1), 1, sum)

##### CT + SM 
samples.fit <- samples.ct[[curr.set]]
samples.pred <- samples.ct.sm[[pred.set]]

access.fit <- c(covars.ct$Access[-ct.x.indices[[curr.set]]])
mean.access.fit <- mean(access.fit)
sd.access.fit <- sd(access.fit)
access.fit <- (access.fit - mean(access.fit)) / sd(access.fit)
access.pred <- c(covars.ct$Access[ct.x.indices[[curr.set]]],
                 covars.sm$Access[sm.x.indices[[curr.set]]])
access.pred <- (access.pred - mean.access.fit) / sd.access.fit 

tree.fit <- c(covars.ct$Tree[-ct.x.indices[[curr.set]]])
mean.tree.fit <- mean(tree.fit)
sd.tree.fit <- sd(tree.fit)
tree.fit <- (tree.fit - mean(tree.fit)) / sd(tree.fit)
tree.pred <- c(covars.ct$Tree[ct.x.indices[[curr.set]]],
               covars.sm$Tree[sm.x.indices[[curr.set]]])
tree.pred <- (tree.pred - mean.tree.fit) / sd.tree.fit 

canopy.fit <- c(covars.ct$Canopy[-ct.x.indices[[curr.set]]])
mean.canopy.fit <- mean(canopy.fit)
sd.canopy.fit <- sd(canopy.fit)
canopy.fit <- (canopy.fit - mean(canopy.fit)) / sd(canopy.fit)
canopy.pred <- c(covars.ct$Canopy[ct.x.indices[[curr.set]]],
                 covars.sm$Canopy[sm.x.indices[[curr.set]]])
canopy.pred <- (canopy.pred - mean.canopy.fit) / sd.canopy.fit 

n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.ct <- length(ct.x.indices[[curr.set]])
J.sm <- length(sm.x.indices[[curr.set]])
access.ct.pred <- access.pred[1:J.ct]
access.sm.pred <- access.pred[1:J.sm + J.ct]
tree.ct.pred <- tree.pred[1:J.ct]
tree.sm.pred <- tree.pred[1:J.sm + J.ct]
canopy.ct.pred <- canopy.pred[1:J.ct]
canopy.sm.pred <- canopy.pred[1:J.sm + J.ct]

ct.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                       samples.pred = samples.pred, 
                       access.pred = access.ct.pred, 
                       tree.pred = tree.ct.pred, 
                       canopy.pred = canopy.ct.pred,
                       my.iter = my.iter, 
                       n.iter = n.iter, 
                       I = I, type = 'ct')
elpd.ct.vals[, 3, 1] <- apply(log(ct.elpd.vals.1), 1, sum)

sm.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                       samples.pred = samples.pred, 
                       access.pred = access.sm.pred, 
                       tree.pred = tree.sm.pred, 
                       canopy.pred = canopy.sm.pred,
                       my.iter = my.iter, 
                       n.iter = n.iter, 
                       I = I, type = 'sm')
elpd.sm.vals[, 2, 1] <- apply(log(sm.elpd.vals.1), 1, sum)

#### TR + SM 
samples.fit <- samples.ct[[curr.set]]
samples.pred <- samples.tr.sm[[pred.set]]

access.fit <- c(covars.ct$Access[-ct.x.indices[[curr.set]]])
mean.access.fit <- mean(access.fit)
sd.access.fit <- sd(access.fit)
access.fit <- (access.fit - mean(access.fit)) / sd(access.fit)
access.pred <- c(covars.tr$Access[tr.x.indices[[curr.set]]],
                 covars.sm$Access[sm.x.indices[[curr.set]]])
access.pred <- (access.pred - mean.access.fit) / sd.access.fit 

tree.fit <- c(covars.ct$Tree[-ct.x.indices[[curr.set]]])
mean.tree.fit <- mean(tree.fit)
sd.tree.fit <- sd(tree.fit)
tree.fit <- (tree.fit - mean(tree.fit)) / sd(tree.fit)
tree.pred <- c(covars.tr$Tree[tr.x.indices[[curr.set]]],
               covars.sm$Tree[sm.x.indices[[curr.set]]])
tree.pred <- (tree.pred - mean.tree.fit) / sd.tree.fit 

canopy.fit <- c(covars.ct$Canopy[-ct.x.indices[[curr.set]]])
mean.canopy.fit <- mean(canopy.fit)
sd.canopy.fit <- sd(canopy.fit)
canopy.fit <- (canopy.fit - mean(canopy.fit)) / sd(canopy.fit)
canopy.pred <- c(covars.tr$Canopy[tr.x.indices[[curr.set]]],
                 covars.sm$Canopy[sm.x.indices[[curr.set]]])
canopy.pred <- (canopy.pred - mean.canopy.fit) / sd.canopy.fit 

n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.tr <- length(tr.x.indices[[curr.set]])
J.sm <- length(sm.x.indices[[curr.set]])
access.tr.pred <- access.pred[1:J.tr]
access.sm.pred <- access.pred[1:J.sm + J.tr]
tree.tr.pred <- tree.pred[1:J.tr]
tree.sm.pred <- tree.pred[1:J.sm + J.tr]
canopy.tr.pred <- canopy.pred[1:J.tr]
canopy.sm.pred <- canopy.pred[1:J.sm + J.tr]

tr.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                       samples.pred = samples.pred, 
                       access.pred = access.tr.pred, 
                       tree.pred = tree.tr.pred, 
                       canopy.pred = canopy.tr.pred,
                       my.iter = my.iter, 
                       n.iter = n.iter, 
                       I = I, type = 'tr')
elpd.tr.vals[, 3, 1] <- apply(log(tr.elpd.vals.1), 1, sum)

sm.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                       samples.pred = samples.pred, 
                       access.pred = access.sm.pred, 
                       tree.pred = tree.sm.pred, 
                       canopy.pred = canopy.sm.pred,
                       my.iter = my.iter, 
                       n.iter = n.iter, 
                       I = I, type = 'sm')
elpd.sm.vals[, 3, 1] <- apply(log(sm.elpd.vals.1), 1, sum)

#### CT + TR + SM 
samples.fit <- samples.ct[[curr.set]]
samples.pred <- samples.ct.tr.sm[[pred.set]]

access.fit <- c(covars.ct$Access[-ct.x.indices[[curr.set]]])
mean.access.fit <- mean(access.fit)
sd.access.fit <- sd(access.fit)
access.fit <- (access.fit - mean(access.fit)) / sd(access.fit)
access.pred <- c(covars.ct$Access[ct.x.indices[[curr.set]]],
                 covars.tr$Access[tr.x.indices[[curr.set]]],
                 covars.sm$Access[sm.x.indices[[curr.set]]])
access.pred <- (access.pred - mean.access.fit) / sd.access.fit 

tree.fit <- c(covars.ct$Tree[-ct.x.indices[[curr.set]]])
mean.tree.fit <- mean(tree.fit)
sd.tree.fit <- sd(tree.fit)
tree.fit <- (tree.fit - mean(tree.fit)) / sd(tree.fit)
tree.pred <- c(covars.ct$Tree[ct.x.indices[[curr.set]]],
               covars.tr$Tree[tr.x.indices[[curr.set]]],
               covars.sm$Tree[sm.x.indices[[curr.set]]])
tree.pred <- (tree.pred - mean.tree.fit) / sd.tree.fit 

canopy.fit <- c(covars.ct$Canopy[-ct.x.indices[[curr.set]]])
mean.canopy.fit <- mean(canopy.fit)
sd.canopy.fit <- sd(canopy.fit)
canopy.fit <- (canopy.fit - mean(canopy.fit)) / sd(canopy.fit)
canopy.pred <- c(covars.ct$Canopy[ct.x.indices[[curr.set]]],
                 covars.tr$Canopy[tr.x.indices[[curr.set]]],
                 covars.sm$Canopy[sm.x.indices[[curr.set]]])
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
elpd.ct.vals[, 4, 1] <- apply(log(ct.elpd.vals.1), 1, sum)

tr.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                       samples.pred = samples.pred, 
                       access.pred = access.tr.pred, 
                       tree.pred = tree.tr.pred, 
                       canopy.pred = canopy.tr.pred,
                       my.iter = my.iter, 
                       n.iter = n.iter, 
                       I = I, type = 'tr')
elpd.tr.vals[, 4, 1] <- apply(log(tr.elpd.vals.1), 1, sum)

sm.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                       samples.pred = samples.pred, 
                       access.pred = access.sm.pred, 
                       tree.pred = tree.sm.pred, 
                       canopy.pred = canopy.sm.pred,
                       my.iter = my.iter, 
                       n.iter = n.iter, 
                       I = I, type = 'sm')
elpd.sm.vals[, 4, 1] <- apply(log(sm.elpd.vals.1), 1, sum)




#############
# TR as predicting model -------------------------------------------------
##### CT 
samples.fit <- samples.tr[[curr.set]]
samples.pred <- samples.ct[[pred.set]]

access.fit <- c(covars.tr$Access[-tr.x.indices[[curr.set]]])
mean.access.fit <- mean(access.fit)
sd.access.fit <- sd(access.fit)
access.fit <- (access.fit - mean(access.fit)) / sd(access.fit)
access.pred <- c(covars.ct$Access[ct.x.indices[[curr.set]]])
access.pred <- (access.pred - mean.access.fit) / sd.access.fit 

tree.fit <- c(covars.tr$Tree[-tr.x.indices[[curr.set]]])
mean.tree.fit <- mean(tree.fit)
sd.tree.fit <- sd(tree.fit)
tree.fit <- (tree.fit - mean(tree.fit)) / sd(tree.fit)
tree.pred <- c(covars.ct$Tree[ct.x.indices[[curr.set]]])
tree.pred <- (tree.pred - mean.tree.fit) / sd.tree.fit 

canopy.fit <- c(covars.tr$Canopy[-tr.x.indices[[curr.set]]])
mean.canopy.fit <- mean(canopy.fit)
sd.canopy.fit <- sd(canopy.fit)
canopy.fit <- (canopy.fit - mean(canopy.fit)) / sd(canopy.fit)
canopy.pred <- c(covars.ct$Canopy[ct.x.indices[[curr.set]]])
canopy.pred <- (canopy.pred - mean.canopy.fit) / sd.canopy.fit 

n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.ct <- length(ct.x.indices[[curr.set]])
elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                    samples.pred = samples.pred, 
                    access.pred = access.pred, 
                    tree.pred = tree.pred, 
                    canopy.pred = canopy.pred,
                    my.iter = my.iter, 
                    n.iter = n.iter, 
                    I = I, type = 'ct')
elpd.ct.vals[, 1, 2] <- apply(log(elpd.vals.1), 1, sum)

##### TR 
samples.fit <- samples.tr[[curr.set]]
samples.pred <- samples.tr[[pred.set]]

access.fit <- c(covars.tr$Access[-tr.x.indices[[curr.set]]])
mean.access.fit <- mean(access.fit)
sd.access.fit <- sd(access.fit)
access.fit <- (access.fit - mean(access.fit)) / sd(access.fit)
access.pred <- c(covars.tr$Access[tr.x.indices[[curr.set]]])
access.pred <- (access.pred - mean.access.fit) / sd.access.fit 

tree.fit <- c(covars.tr$Tree[-tr.x.indices[[curr.set]]])
mean.tree.fit <- mean(tree.fit)
sd.tree.fit <- sd(tree.fit)
tree.fit <- (tree.fit - mean(tree.fit)) / sd(tree.fit)
tree.pred <- c(covars.tr$Tree[tr.x.indices[[curr.set]]])
tree.pred <- (tree.pred - mean.tree.fit) / sd.tree.fit 

canopy.fit <- c(covars.tr$Canopy[-tr.x.indices[[curr.set]]])
mean.canopy.fit <- mean(canopy.fit)
sd.canopy.fit <- sd(canopy.fit)
canopy.fit <- (canopy.fit - mean(canopy.fit)) / sd(canopy.fit)
canopy.pred <- c(covars.tr$Canopy[tr.x.indices[[curr.set]]])
canopy.pred <- (canopy.pred - mean.canopy.fit) / sd.canopy.fit 

n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.tr <- length(tr.x.indices[[curr.set]])
elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                    samples.pred = samples.pred, 
                    access.pred = access.pred, 
                    tree.pred = tree.pred, 
                    canopy.pred = canopy.pred,
                    my.iter = my.iter, 
                    n.iter = n.iter, 
                    I = I, type = 'tr')
elpd.tr.vals[, 1, 2] <- apply(log(elpd.vals.1), 1, sum)

##### SM 
samples.fit <- samples.tr[[curr.set]]
samples.pred <- samples.sm[[pred.set]]

access.fit <- c(covars.tr$Access[-tr.x.indices[[curr.set]]])
mean.access.fit <- mean(access.fit)
sd.access.fit <- sd(access.fit)
access.fit <- (access.fit - mean(access.fit)) / sd(access.fit)
access.pred <- c(covars.sm$Access[sm.x.indices[[curr.set]]])
access.pred <- (access.pred - mean.access.fit) / sd.access.fit 

tree.fit <- c(covars.tr$Tree[-tr.x.indices[[curr.set]]])
mean.tree.fit <- mean(tree.fit)
sd.tree.fit <- sd(tree.fit)
tree.fit <- (tree.fit - mean(tree.fit)) / sd(tree.fit)
tree.pred <- c(covars.sm$Tree[sm.x.indices[[curr.set]]])
tree.pred <- (tree.pred - mean.tree.fit) / sd.tree.fit 

canopy.fit <- c(covars.tr$Canopy[-tr.x.indices[[curr.set]]])
mean.canopy.fit <- mean(canopy.fit)
sd.canopy.fit <- sd(canopy.fit)
canopy.fit <- (canopy.fit - mean(canopy.fit)) / sd(canopy.fit)
canopy.pred <- c(covars.sm$Canopy[sm.x.indices[[curr.set]]])
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
elpd.sm.vals[, 1, 2] <- apply(log(elpd.vals.1), 1, sum)

##### CT + TR 
samples.fit <- samples.tr[[curr.set]]
samples.pred <- samples.ct.tr[[pred.set]]

access.fit <- c(covars.tr$Access[-tr.x.indices[[curr.set]]])
mean.access.fit <- mean(access.fit)
sd.access.fit <- sd(access.fit)
access.fit <- (access.fit - mean(access.fit)) / sd(access.fit)
access.pred <- c(covars.ct$Access[ct.x.indices[[curr.set]]],
                 covars.tr$Access[tr.x.indices[[curr.set]]])
access.pred <- (access.pred - mean.access.fit) / sd.access.fit 

tree.fit <- c(covars.tr$Tree[-tr.x.indices[[curr.set]]])
mean.tree.fit <- mean(tree.fit)
sd.tree.fit <- sd(tree.fit)
tree.fit <- (tree.fit - mean(tree.fit)) / sd(tree.fit)
tree.pred <- c(covars.ct$Tree[ct.x.indices[[curr.set]]],
               covars.tr$Tree[tr.x.indices[[curr.set]]])
tree.pred <- (tree.pred - mean.tree.fit) / sd.tree.fit 

canopy.fit <- c(covars.tr$Canopy[-tr.x.indices[[curr.set]]])
mean.canopy.fit <- mean(canopy.fit)
sd.canopy.fit <- sd(canopy.fit)
canopy.fit <- (canopy.fit - mean(canopy.fit)) / sd(canopy.fit)
canopy.pred <- c(covars.ct$Canopy[ct.x.indices[[curr.set]]],
                 covars.tr$Canopy[tr.x.indices[[curr.set]]])
canopy.pred <- (canopy.pred - mean.canopy.fit) / sd.canopy.fit 

n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.ct <- length(ct.x.indices[[curr.set]])
J.tr <- length(tr.x.indices[[curr.set]])
access.ct.pred <- access.pred[1:J.ct]
access.tr.pred <- access.pred[1:J.tr + J.ct]
tree.ct.pred <- tree.pred[1:J.ct]
tree.tr.pred <- tree.pred[1:J.tr + J.ct]
canopy.ct.pred <- canopy.pred[1:J.ct]
canopy.tr.pred <- canopy.pred[1:J.tr + J.ct]

ct.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                       samples.pred = samples.pred, 
                       access.pred = access.ct.pred, 
                       tree.pred = tree.ct.pred, 
                       canopy.pred = canopy.ct.pred,
                       my.iter = my.iter, 
                       n.iter = n.iter, 
                       I = I, type = 'ct')
elpd.ct.vals[, 2, 2] <- apply(log(ct.elpd.vals.1), 1, sum)

tr.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                       samples.pred = samples.pred, 
                       access.pred = access.tr.pred, 
                       tree.pred = tree.tr.pred, 
                       canopy.pred = canopy.tr.pred,
                       my.iter = my.iter, 
                       n.iter = n.iter, 
                       I = I, type = 'tr')
elpd.tr.vals[, 2, 2] <- apply(log(tr.elpd.vals.1), 1, sum)

##### CT + SM 
samples.fit <- samples.tr[[curr.set]]
samples.pred <- samples.ct.sm[[pred.set]]

access.fit <- c(covars.tr$Access[-tr.x.indices[[curr.set]]])
mean.access.fit <- mean(access.fit)
sd.access.fit <- sd(access.fit)
access.fit <- (access.fit - mean(access.fit)) / sd(access.fit)
access.pred <- c(covars.ct$Access[ct.x.indices[[curr.set]]],
                 covars.sm$Access[sm.x.indices[[curr.set]]])
access.pred <- (access.pred - mean.access.fit) / sd.access.fit 

tree.fit <- c(covars.tr$Tree[-tr.x.indices[[curr.set]]])
mean.tree.fit <- mean(tree.fit)
sd.tree.fit <- sd(tree.fit)
tree.fit <- (tree.fit - mean(tree.fit)) / sd(tree.fit)
tree.pred <- c(covars.ct$Tree[ct.x.indices[[curr.set]]],
               covars.sm$Tree[sm.x.indices[[curr.set]]])
tree.pred <- (tree.pred - mean.tree.fit) / sd.tree.fit 

canopy.fit <- c(covars.tr$Canopy[-tr.x.indices[[curr.set]]])
mean.canopy.fit <- mean(canopy.fit)
sd.canopy.fit <- sd(canopy.fit)
canopy.fit <- (canopy.fit - mean(canopy.fit)) / sd(canopy.fit)
canopy.pred <- c(covars.ct$Canopy[ct.x.indices[[curr.set]]],
                 covars.sm$Canopy[sm.x.indices[[curr.set]]])
canopy.pred <- (canopy.pred - mean.canopy.fit) / sd.canopy.fit 

n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.ct <- length(ct.x.indices[[curr.set]])
J.sm <- length(sm.x.indices[[curr.set]])
access.ct.pred <- access.pred[1:J.ct]
access.sm.pred <- access.pred[1:J.sm + J.ct]
tree.ct.pred <- tree.pred[1:J.ct]
tree.sm.pred <- tree.pred[1:J.sm + J.ct]
canopy.ct.pred <- canopy.pred[1:J.ct]
canopy.sm.pred <- canopy.pred[1:J.sm + J.ct]

ct.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                       samples.pred = samples.pred, 
                       access.pred = access.ct.pred, 
                       tree.pred = tree.ct.pred, 
                       canopy.pred = canopy.ct.pred,
                       my.iter = my.iter, 
                       n.iter = n.iter, 
                       I = I, type = 'ct')
elpd.ct.vals[, 3, 2] <- apply(log(ct.elpd.vals.1), 1, sum)

sm.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                       samples.pred = samples.pred, 
                       access.pred = access.sm.pred, 
                       tree.pred = tree.sm.pred, 
                       canopy.pred = canopy.sm.pred,
                       my.iter = my.iter, 
                       n.iter = n.iter, 
                       I = I, type = 'sm')
elpd.sm.vals[, 2, 2] <- apply(log(sm.elpd.vals.1), 1, sum)

#### TR + SM 
samples.fit <- samples.tr[[curr.set]]
samples.pred <- samples.tr.sm[[pred.set]]

access.fit <- c(covars.tr$Access[-tr.x.indices[[curr.set]]])
mean.access.fit <- mean(access.fit)
sd.access.fit <- sd(access.fit)
access.fit <- (access.fit - mean(access.fit)) / sd(access.fit)
access.pred <- c(covars.tr$Access[tr.x.indices[[curr.set]]],
                 covars.sm$Access[sm.x.indices[[curr.set]]])
access.pred <- (access.pred - mean.access.fit) / sd.access.fit 

tree.fit <- c(covars.tr$Tree[-tr.x.indices[[curr.set]]])
mean.tree.fit <- mean(tree.fit)
sd.tree.fit <- sd(tree.fit)
tree.fit <- (tree.fit - mean(tree.fit)) / sd(tree.fit)
tree.pred <- c(covars.tr$Tree[tr.x.indices[[curr.set]]],
               covars.sm$Tree[sm.x.indices[[curr.set]]])
tree.pred <- (tree.pred - mean.tree.fit) / sd.tree.fit 

canopy.fit <- c(covars.tr$Canopy[-tr.x.indices[[curr.set]]])
mean.canopy.fit <- mean(canopy.fit)
sd.canopy.fit <- sd(canopy.fit)
canopy.fit <- (canopy.fit - mean(canopy.fit)) / sd(canopy.fit)
canopy.pred <- c(covars.tr$Canopy[tr.x.indices[[curr.set]]],
                 covars.sm$Canopy[sm.x.indices[[curr.set]]])
canopy.pred <- (canopy.pred - mean.canopy.fit) / sd.canopy.fit 

n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.tr <- length(tr.x.indices[[curr.set]])
J.sm <- length(sm.x.indices[[curr.set]])
access.tr.pred <- access.pred[1:J.tr]
access.sm.pred <- access.pred[1:J.sm + J.tr]
tree.tr.pred <- tree.pred[1:J.tr]
tree.sm.pred <- tree.pred[1:J.sm + J.tr]
canopy.tr.pred <- canopy.pred[1:J.tr]
canopy.sm.pred <- canopy.pred[1:J.sm + J.tr]

tr.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                       samples.pred = samples.pred, 
                       access.pred = access.tr.pred, 
                       tree.pred = tree.tr.pred, 
                       canopy.pred = canopy.tr.pred,
                       my.iter = my.iter, 
                       n.iter = n.iter, 
                       I = I, type = 'tr')
elpd.tr.vals[, 3, 2] <- apply(log(tr.elpd.vals.1), 1, sum)

sm.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                       samples.pred = samples.pred, 
                       access.pred = access.sm.pred, 
                       tree.pred = tree.sm.pred, 
                       canopy.pred = canopy.sm.pred,
                       my.iter = my.iter, 
                       n.iter = n.iter, 
                       I = I, type = 'sm')
elpd.sm.vals[, 3, 2] <- apply(log(sm.elpd.vals.1), 1, sum)

#### CT + TR + SM 
samples.fit <- samples.tr[[curr.set]]
samples.pred <- samples.ct.tr.sm[[pred.set]]

access.fit <- c(covars.tr$Access[-tr.x.indices[[curr.set]]])
mean.access.fit <- mean(access.fit)
sd.access.fit <- sd(access.fit)
access.fit <- (access.fit - mean(access.fit)) / sd(access.fit)
access.pred <- c(covars.ct$Access[ct.x.indices[[curr.set]]],
                 covars.tr$Access[tr.x.indices[[curr.set]]],
                 covars.sm$Access[sm.x.indices[[curr.set]]])
access.pred <- (access.pred - mean.access.fit) / sd.access.fit 

tree.fit <- c(covars.tr$Tree[-tr.x.indices[[curr.set]]])
mean.tree.fit <- mean(tree.fit)
sd.tree.fit <- sd(tree.fit)
tree.fit <- (tree.fit - mean(tree.fit)) / sd(tree.fit)
tree.pred <- c(covars.ct$Tree[ct.x.indices[[curr.set]]],
               covars.tr$Tree[tr.x.indices[[curr.set]]],
               covars.sm$Tree[sm.x.indices[[curr.set]]])
tree.pred <- (tree.pred - mean.tree.fit) / sd.tree.fit 

canopy.fit <- c(covars.tr$Canopy[-tr.x.indices[[curr.set]]])
mean.canopy.fit <- mean(canopy.fit)
sd.canopy.fit <- sd(canopy.fit)
canopy.fit <- (canopy.fit - mean(canopy.fit)) / sd(canopy.fit)
canopy.pred <- c(covars.ct$Canopy[ct.x.indices[[curr.set]]],
                 covars.tr$Canopy[tr.x.indices[[curr.set]]],
                 covars.sm$Canopy[sm.x.indices[[curr.set]]])
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
elpd.ct.vals[, 4, 2] <- apply(log(ct.elpd.vals.1), 1, sum)

tr.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                       samples.pred = samples.pred, 
                       access.pred = access.tr.pred, 
                       tree.pred = tree.tr.pred, 
                       canopy.pred = canopy.tr.pred,
                       my.iter = my.iter, 
                       n.iter = n.iter, 
                       I = I, type = 'tr')
elpd.tr.vals[, 4, 2] <- apply(log(tr.elpd.vals.1), 1, sum)

sm.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                       samples.pred = samples.pred, 
                       access.pred = access.sm.pred, 
                       tree.pred = tree.sm.pred, 
                       canopy.pred = canopy.sm.pred,
                       my.iter = my.iter, 
                       n.iter = n.iter, 
                       I = I, type = 'sm')
elpd.sm.vals[, 4, 2] <- apply(log(sm.elpd.vals.1), 1, sum)



#############
# SM as predicting model -------------------------------------------------
##### CT 
samples.fit <- samples.sm[[curr.set]]
samples.pred <- samples.ct[[pred.set]]

access.fit <- c(covars.sm$Access[-sm.x.indices[[curr.set]]])
mean.access.fit <- mean(access.fit)
sd.access.fit <- sd(access.fit)
access.fit <- (access.fit - mean(access.fit)) / sd(access.fit)
access.pred <- c(covars.ct$Access[ct.x.indices[[curr.set]]])
access.pred <- (access.pred - mean.access.fit) / sd.access.fit 

tree.fit <- c(covars.sm$Tree[-sm.x.indices[[curr.set]]])
mean.tree.fit <- mean(tree.fit)
sd.tree.fit <- sd(tree.fit)
tree.fit <- (tree.fit - mean(tree.fit)) / sd(tree.fit)
tree.pred <- c(covars.ct$Tree[ct.x.indices[[curr.set]]])
tree.pred <- (tree.pred - mean.tree.fit) / sd.tree.fit 

canopy.fit <- c(covars.sm$Canopy[-sm.x.indices[[curr.set]]])
mean.canopy.fit <- mean(canopy.fit)
sd.canopy.fit <- sd(canopy.fit)
canopy.fit <- (canopy.fit - mean(canopy.fit)) / sd(canopy.fit)
canopy.pred <- c(covars.ct$Canopy[ct.x.indices[[curr.set]]])
canopy.pred <- (canopy.pred - mean.canopy.fit) / sd.canopy.fit 

n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.ct <- length(ct.x.indices[[curr.set]])
elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                    samples.pred = samples.pred, 
                    access.pred = access.pred, 
                    tree.pred = tree.pred, 
                    canopy.pred = canopy.pred,
                    my.iter = my.iter, 
                    n.iter = n.iter, 
                    I = I, type = 'ct')
elpd.ct.vals[, 1, 3] <- apply(log(elpd.vals.1), 1, sum)

##### TR 
samples.fit <- samples.sm[[curr.set]]
samples.pred <- samples.tr[[pred.set]]

access.fit <- c(covars.sm$Access[-sm.x.indices[[curr.set]]])
mean.access.fit <- mean(access.fit)
sd.access.fit <- sd(access.fit)
access.fit <- (access.fit - mean(access.fit)) / sd(access.fit)
access.pred <- c(covars.tr$Access[tr.x.indices[[curr.set]]])
access.pred <- (access.pred - mean.access.fit) / sd.access.fit 

tree.fit <- c(covars.sm$Tree[-sm.x.indices[[curr.set]]])
mean.tree.fit <- mean(tree.fit)
sd.tree.fit <- sd(tree.fit)
tree.fit <- (tree.fit - mean(tree.fit)) / sd(tree.fit)
tree.pred <- c(covars.tr$Tree[tr.x.indices[[curr.set]]])
tree.pred <- (tree.pred - mean.tree.fit) / sd.tree.fit 

canopy.fit <- c(covars.sm$Canopy[-sm.x.indices[[curr.set]]])
mean.canopy.fit <- mean(canopy.fit)
sd.canopy.fit <- sd(canopy.fit)
canopy.fit <- (canopy.fit - mean(canopy.fit)) / sd(canopy.fit)
canopy.pred <- c(covars.tr$Canopy[tr.x.indices[[curr.set]]])
canopy.pred <- (canopy.pred - mean.canopy.fit) / sd.canopy.fit 

n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.tr <- length(tr.x.indices[[curr.set]])
elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                    samples.pred = samples.pred, 
                    access.pred = access.pred, 
                    tree.pred = tree.pred, 
                    canopy.pred = canopy.pred,
                    my.iter = my.iter, 
                    n.iter = n.iter, 
                    I = I, type = 'tr')
elpd.tr.vals[, 1, 3] <- apply(log(elpd.vals.1), 1, sum)

##### SM 
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

canopy.fit <- c(covars.sm$Canopy[-sm.x.indices[[curr.set]]])
mean.canopy.fit <- mean(canopy.fit)
sd.canopy.fit <- sd(canopy.fit)
canopy.fit <- (canopy.fit - mean(canopy.fit)) / sd(canopy.fit)
canopy.pred <- c(covars.sm$Canopy[sm.x.indices[[curr.set]]])
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
elpd.sm.vals[, 1, 3] <- apply(log(elpd.vals.1), 1, sum)

##### CT + TR 
samples.fit <- samples.sm[[curr.set]]
samples.pred <- samples.ct.tr[[pred.set]]

access.fit <- c(covars.sm$Access[-sm.x.indices[[curr.set]]])
mean.access.fit <- mean(access.fit)
sd.access.fit <- sd(access.fit)
access.fit <- (access.fit - mean(access.fit)) / sd(access.fit)
access.pred <- c(covars.ct$Access[ct.x.indices[[curr.set]]],
                 covars.tr$Access[tr.x.indices[[curr.set]]])
access.pred <- (access.pred - mean.access.fit) / sd.access.fit 

tree.fit <- c(covars.sm$Tree[-sm.x.indices[[curr.set]]])
mean.tree.fit <- mean(tree.fit)
sd.tree.fit <- sd(tree.fit)
tree.fit <- (tree.fit - mean(tree.fit)) / sd(tree.fit)
tree.pred <- c(covars.ct$Tree[ct.x.indices[[curr.set]]],
               covars.tr$Tree[tr.x.indices[[curr.set]]])
tree.pred <- (tree.pred - mean.tree.fit) / sd.tree.fit 

canopy.fit <- c(covars.sm$Canopy[-sm.x.indices[[curr.set]]])
mean.canopy.fit <- mean(canopy.fit)
sd.canopy.fit <- sd(canopy.fit)
canopy.fit <- (canopy.fit - mean(canopy.fit)) / sd(canopy.fit)
canopy.pred <- c(covars.ct$Canopy[ct.x.indices[[curr.set]]],
                 covars.tr$Canopy[tr.x.indices[[curr.set]]])
canopy.pred <- (canopy.pred - mean.canopy.fit) / sd.canopy.fit 

n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.ct <- length(ct.x.indices[[curr.set]])
J.tr <- length(tr.x.indices[[curr.set]])
access.ct.pred <- access.pred[1:J.ct]
access.tr.pred <- access.pred[1:J.tr + J.ct]
tree.ct.pred <- tree.pred[1:J.ct]
tree.tr.pred <- tree.pred[1:J.tr + J.ct]
canopy.ct.pred <- canopy.pred[1:J.ct]
canopy.tr.pred <- canopy.pred[1:J.tr + J.ct]

ct.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                       samples.pred = samples.pred, 
                       access.pred = access.ct.pred, 
                       tree.pred = tree.ct.pred, 
                       canopy.pred = canopy.ct.pred,
                       my.iter = my.iter, 
                       n.iter = n.iter, 
                       I = I, type = 'ct')
elpd.ct.vals[, 2, 3] <- apply(log(ct.elpd.vals.1), 1, sum)

tr.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                       samples.pred = samples.pred, 
                       access.pred = access.tr.pred, 
                       tree.pred = tree.tr.pred, 
                       canopy.pred = canopy.tr.pred,
                       my.iter = my.iter, 
                       n.iter = n.iter, 
                       I = I, type = 'tr')
elpd.tr.vals[, 2, 3] <- apply(log(tr.elpd.vals.1), 1, sum)

##### CT + SM 
samples.fit <- samples.sm[[curr.set]]
samples.pred <- samples.ct.sm[[pred.set]]

access.fit <- c(covars.sm$Access[-sm.x.indices[[curr.set]]])
mean.access.fit <- mean(access.fit)
sd.access.fit <- sd(access.fit)
access.fit <- (access.fit - mean(access.fit)) / sd(access.fit)
access.pred <- c(covars.ct$Access[ct.x.indices[[curr.set]]],
                 covars.sm$Access[sm.x.indices[[curr.set]]])
access.pred <- (access.pred - mean.access.fit) / sd.access.fit 

tree.fit <- c(covars.sm$Tree[-sm.x.indices[[curr.set]]])
mean.tree.fit <- mean(tree.fit)
sd.tree.fit <- sd(tree.fit)
tree.fit <- (tree.fit - mean(tree.fit)) / sd(tree.fit)
tree.pred <- c(covars.ct$Tree[ct.x.indices[[curr.set]]],
               covars.sm$Tree[sm.x.indices[[curr.set]]])
tree.pred <- (tree.pred - mean.tree.fit) / sd.tree.fit 

canopy.fit <- c(covars.sm$Canopy[-sm.x.indices[[curr.set]]])
mean.canopy.fit <- mean(canopy.fit)
sd.canopy.fit <- sd(canopy.fit)
canopy.fit <- (canopy.fit - mean(canopy.fit)) / sd(canopy.fit)
canopy.pred <- c(covars.ct$Canopy[ct.x.indices[[curr.set]]],
                 covars.sm$Canopy[sm.x.indices[[curr.set]]])
canopy.pred <- (canopy.pred - mean.canopy.fit) / sd.canopy.fit 

n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.ct <- length(ct.x.indices[[curr.set]])
J.sm <- length(sm.x.indices[[curr.set]])
access.ct.pred <- access.pred[1:J.ct]
access.sm.pred <- access.pred[1:J.sm + J.ct]
tree.ct.pred <- tree.pred[1:J.ct]
tree.sm.pred <- tree.pred[1:J.sm + J.ct]
canopy.ct.pred <- canopy.pred[1:J.ct]
canopy.sm.pred <- canopy.pred[1:J.sm + J.ct]

ct.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                       samples.pred = samples.pred, 
                       access.pred = access.ct.pred, 
                       tree.pred = tree.ct.pred, 
                       canopy.pred = canopy.ct.pred,
                       my.iter = my.iter, 
                       n.iter = n.iter, 
                       I = I, type = 'ct')
elpd.ct.vals[, 3, 3] <- apply(log(ct.elpd.vals.1), 1, sum)

sm.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                       samples.pred = samples.pred, 
                       access.pred = access.sm.pred, 
                       tree.pred = tree.sm.pred, 
                       canopy.pred = canopy.sm.pred,
                       my.iter = my.iter, 
                       n.iter = n.iter, 
                       I = I, type = 'sm')
elpd.sm.vals[, 2, 3] <- apply(log(sm.elpd.vals.1), 1, sum)

#### TR + SM 
samples.fit <- samples.sm[[curr.set]]
samples.pred <- samples.tr.sm[[pred.set]]

access.fit <- c(covars.sm$Access[-sm.x.indices[[curr.set]]])
mean.access.fit <- mean(access.fit)
sd.access.fit <- sd(access.fit)
access.fit <- (access.fit - mean(access.fit)) / sd(access.fit)
access.pred <- c(covars.tr$Access[tr.x.indices[[curr.set]]],
                 covars.sm$Access[sm.x.indices[[curr.set]]])
access.pred <- (access.pred - mean.access.fit) / sd.access.fit 

tree.fit <- c(covars.sm$Tree[-sm.x.indices[[curr.set]]])
mean.tree.fit <- mean(tree.fit)
sd.tree.fit <- sd(tree.fit)
tree.fit <- (tree.fit - mean(tree.fit)) / sd(tree.fit)
tree.pred <- c(covars.tr$Tree[tr.x.indices[[curr.set]]],
               covars.sm$Tree[sm.x.indices[[curr.set]]])
tree.pred <- (tree.pred - mean.tree.fit) / sd.tree.fit 

canopy.fit <- c(covars.sm$Canopy[-sm.x.indices[[curr.set]]])
mean.canopy.fit <- mean(canopy.fit)
sd.canopy.fit <- sd(canopy.fit)
canopy.fit <- (canopy.fit - mean(canopy.fit)) / sd(canopy.fit)
canopy.pred <- c(covars.tr$Canopy[tr.x.indices[[curr.set]]],
                 covars.sm$Canopy[sm.x.indices[[curr.set]]])
canopy.pred <- (canopy.pred - mean.canopy.fit) / sd.canopy.fit 

n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.tr <- length(tr.x.indices[[curr.set]])
J.sm <- length(sm.x.indices[[curr.set]])
access.tr.pred <- access.pred[1:J.tr]
access.sm.pred <- access.pred[1:J.sm + J.tr]
tree.tr.pred <- tree.pred[1:J.tr]
tree.sm.pred <- tree.pred[1:J.sm + J.tr]
canopy.tr.pred <- canopy.pred[1:J.tr]
canopy.sm.pred <- canopy.pred[1:J.sm + J.tr]

tr.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                       samples.pred = samples.pred, 
                       access.pred = access.tr.pred, 
                       tree.pred = tree.tr.pred, 
                       canopy.pred = canopy.tr.pred,
                       my.iter = my.iter, 
                       n.iter = n.iter, 
                       I = I, type = 'tr')
elpd.tr.vals[, 3, 3] <- apply(log(tr.elpd.vals.1), 1, sum)

sm.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                       samples.pred = samples.pred, 
                       access.pred = access.sm.pred, 
                       tree.pred = tree.sm.pred, 
                       canopy.pred = canopy.sm.pred,
                       my.iter = my.iter, 
                       n.iter = n.iter, 
                       I = I, type = 'sm')
elpd.sm.vals[, 3, 3] <- apply(log(sm.elpd.vals.1), 1, sum)

#### CT + TR + SM 
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

canopy.fit <- c(covars.sm$Canopy[-sm.x.indices[[curr.set]]])
mean.canopy.fit <- mean(canopy.fit)
sd.canopy.fit <- sd(canopy.fit)
canopy.fit <- (canopy.fit - mean(canopy.fit)) / sd(canopy.fit)
canopy.pred <- c(covars.ct$Canopy[ct.x.indices[[curr.set]]],
                 covars.tr$Canopy[tr.x.indices[[curr.set]]],
                 covars.sm$Canopy[sm.x.indices[[curr.set]]])
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
elpd.ct.vals[, 4, 3] <- apply(log(ct.elpd.vals.1), 1, sum)

tr.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                       samples.pred = samples.pred, 
                       access.pred = access.tr.pred, 
                       tree.pred = tree.tr.pred, 
                       canopy.pred = canopy.tr.pred,
                       my.iter = my.iter, 
                       n.iter = n.iter, 
                       I = I, type = 'tr')
elpd.tr.vals[, 4, 3] <- apply(log(tr.elpd.vals.1), 1, sum)

sm.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                       samples.pred = samples.pred, 
                       access.pred = access.sm.pred, 
                       tree.pred = tree.sm.pred, 
                       canopy.pred = canopy.sm.pred,
                       my.iter = my.iter, 
                       n.iter = n.iter, 
                       I = I, type = 'sm')
elpd.sm.vals[, 4, 3] <- apply(log(sm.elpd.vals.1), 1, sum)



#############
# CT-TR as predicting model -------------------------------------------------
##### CT 
samples.fit <- samples.sm[[curr.set]]
samples.pred <- samples.ct[[pred.set]]

access.fit <- c(covars.ct$Access[-ct.x.indices[[curr.set]]],
                covars.tr$Access[-tr.x.indices[[curr.set]]])
mean.access.fit <- mean(access.fit)
sd.access.fit <- sd(access.fit)
access.fit <- (access.fit - mean(access.fit)) / sd(access.fit)
access.pred <- c(covars.ct$Access[ct.x.indices[[curr.set]]])
access.pred <- (access.pred - mean.access.fit) / sd.access.fit 

tree.fit <- c(covars.ct$Tree[-ct.x.indices[[curr.set]]],
              covars.tr$Tree[-tr.x.indices[[curr.set]]])
mean.tree.fit <- mean(tree.fit)
sd.tree.fit <- sd(tree.fit)
tree.fit <- (tree.fit - mean(tree.fit)) / sd(tree.fit)
tree.pred <- c(covars.ct$Tree[ct.x.indices[[curr.set]]])
tree.pred <- (tree.pred - mean.tree.fit) / sd.tree.fit 

canopy.fit <- c(covars.ct$Canopy[-ct.x.indices[[curr.set]]],
                covars.tr$Canopy[-tr.x.indices[[curr.set]]])
mean.canopy.fit <- mean(canopy.fit)
sd.canopy.fit <- sd(canopy.fit)
canopy.fit <- (canopy.fit - mean(canopy.fit)) / sd(canopy.fit)
canopy.pred <- c(covars.ct$Canopy[ct.x.indices[[curr.set]]])
canopy.pred <- (canopy.pred - mean.canopy.fit) / sd.canopy.fit 

n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.ct <- length(ct.x.indices[[curr.set]])
elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                    samples.pred = samples.pred, 
                    access.pred = access.pred, 
                    tree.pred = tree.pred, 
                    canopy.pred = canopy.pred,
                    my.iter = my.iter, 
                    n.iter = n.iter, 
                    I = I, type = 'ct')
elpd.ct.vals[, 1, 4] <- apply(log(elpd.vals.1), 1, sum)

##### TR 
samples.fit <- samples.sm[[curr.set]]
samples.pred <- samples.tr[[pred.set]]

access.fit <- c(covars.ct$Access[-ct.x.indices[[curr.set]]],
                covars.tr$Access[-tr.x.indices[[curr.set]]])
mean.access.fit <- mean(access.fit)
sd.access.fit <- sd(access.fit)
access.fit <- (access.fit - mean(access.fit)) / sd(access.fit)
access.pred <- c(covars.tr$Access[tr.x.indices[[curr.set]]])
access.pred <- (access.pred - mean.access.fit) / sd.access.fit 

tree.fit <- c(covars.ct$Tree[-ct.x.indices[[curr.set]]],
              covars.tr$Tree[-tr.x.indices[[curr.set]]])
mean.tree.fit <- mean(tree.fit)
sd.tree.fit <- sd(tree.fit)
tree.fit <- (tree.fit - mean(tree.fit)) / sd(tree.fit)
tree.pred <- c(covars.tr$Tree[tr.x.indices[[curr.set]]])
tree.pred <- (tree.pred - mean.tree.fit) / sd.tree.fit 

canopy.fit <- c(covars.ct$Canopy[-ct.x.indices[[curr.set]]],
                covars.tr$Canopy[-tr.x.indices[[curr.set]]])
mean.canopy.fit <- mean(canopy.fit)
sd.canopy.fit <- sd(canopy.fit)
canopy.fit <- (canopy.fit - mean(canopy.fit)) / sd(canopy.fit)
canopy.pred <- c(covars.tr$Canopy[tr.x.indices[[curr.set]]])
canopy.pred <- (canopy.pred - mean.canopy.fit) / sd.canopy.fit 

n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.tr <- length(tr.x.indices[[curr.set]])
elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                    samples.pred = samples.pred, 
                    access.pred = access.pred, 
                    tree.pred = tree.pred, 
                    canopy.pred = canopy.pred,
                    my.iter = my.iter, 
                    n.iter = n.iter, 
                    I = I, type = 'tr')
elpd.tr.vals[, 1, 4] <- apply(log(elpd.vals.1), 1, sum)

##### SM 
samples.fit <- samples.sm[[curr.set]]
samples.pred <- samples.sm[[pred.set]]

access.fit <- c(covars.ct$Access[-ct.x.indices[[curr.set]]],
                covars.tr$Access[-tr.x.indices[[curr.set]]])
mean.access.fit <- mean(access.fit)
sd.access.fit <- sd(access.fit)
access.fit <- (access.fit - mean(access.fit)) / sd(access.fit)
access.pred <- c(covars.sm$Access[sm.x.indices[[curr.set]]])
access.pred <- (access.pred - mean.access.fit) / sd.access.fit 

tree.fit <- c(covars.ct$Tree[-ct.x.indices[[curr.set]]],
              covars.tr$Tree[-tr.x.indices[[curr.set]]])
mean.tree.fit <- mean(tree.fit)
sd.tree.fit <- sd(tree.fit)
tree.fit <- (tree.fit - mean(tree.fit)) / sd(tree.fit)
tree.pred <- c(covars.sm$Tree[sm.x.indices[[curr.set]]])
tree.pred <- (tree.pred - mean.tree.fit) / sd.tree.fit 

canopy.fit <- c(covars.ct$Canopy[-ct.x.indices[[curr.set]]],
                covars.tr$Canopy[-tr.x.indices[[curr.set]]])
mean.canopy.fit <- mean(canopy.fit)
sd.canopy.fit <- sd(canopy.fit)
canopy.fit <- (canopy.fit - mean(canopy.fit)) / sd(canopy.fit)
canopy.pred <- c(covars.sm$Canopy[sm.x.indices[[curr.set]]])
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
elpd.sm.vals[, 1, 4] <- apply(log(elpd.vals.1), 1, sum)

##### CT + TR 
samples.fit <- samples.sm[[curr.set]]
samples.pred <- samples.ct.tr[[pred.set]]

access.fit <- c(covars.ct$Access[-ct.x.indices[[curr.set]]],
                covars.tr$Access[-tr.x.indices[[curr.set]]])
mean.access.fit <- mean(access.fit)
sd.access.fit <- sd(access.fit)
access.fit <- (access.fit - mean(access.fit)) / sd(access.fit)
access.pred <- c(covars.ct$Access[ct.x.indices[[curr.set]]],
                 covars.tr$Access[tr.x.indices[[curr.set]]])
access.pred <- (access.pred - mean.access.fit) / sd.access.fit 

tree.fit <- c(covars.ct$Tree[-ct.x.indices[[curr.set]]],
              covars.tr$Tree[-tr.x.indices[[curr.set]]])
mean.tree.fit <- mean(tree.fit)
sd.tree.fit <- sd(tree.fit)
tree.fit <- (tree.fit - mean(tree.fit)) / sd(tree.fit)
tree.pred <- c(covars.ct$Tree[ct.x.indices[[curr.set]]],
               covars.tr$Tree[tr.x.indices[[curr.set]]])
tree.pred <- (tree.pred - mean.tree.fit) / sd.tree.fit 

canopy.fit <- c(covars.ct$Canopy[-ct.x.indices[[curr.set]]],
                covars.tr$Canopy[-tr.x.indices[[curr.set]]])
mean.canopy.fit <- mean(canopy.fit)
sd.canopy.fit <- sd(canopy.fit)
canopy.fit <- (canopy.fit - mean(canopy.fit)) / sd(canopy.fit)
canopy.pred <- c(covars.ct$Canopy[ct.x.indices[[curr.set]]],
                 covars.tr$Canopy[tr.x.indices[[curr.set]]])
canopy.pred <- (canopy.pred - mean.canopy.fit) / sd.canopy.fit 

n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.ct <- length(ct.x.indices[[curr.set]])
J.tr <- length(tr.x.indices[[curr.set]])
access.ct.pred <- access.pred[1:J.ct]
access.tr.pred <- access.pred[1:J.tr + J.ct]
tree.ct.pred <- tree.pred[1:J.ct]
tree.tr.pred <- tree.pred[1:J.tr + J.ct]
canopy.ct.pred <- canopy.pred[1:J.ct]
canopy.tr.pred <- canopy.pred[1:J.tr + J.ct]

ct.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                       samples.pred = samples.pred, 
                       access.pred = access.ct.pred, 
                       tree.pred = tree.ct.pred, 
                       canopy.pred = canopy.ct.pred,
                       my.iter = my.iter, 
                       n.iter = n.iter, 
                       I = I, type = 'ct')
elpd.ct.vals[, 2, 4] <- apply(log(ct.elpd.vals.1), 1, sum)

tr.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                       samples.pred = samples.pred, 
                       access.pred = access.tr.pred, 
                       tree.pred = tree.tr.pred, 
                       canopy.pred = canopy.tr.pred,
                       my.iter = my.iter, 
                       n.iter = n.iter, 
                       I = I, type = 'tr')
elpd.tr.vals[, 2, 4] <- apply(log(tr.elpd.vals.1), 1, sum)

##### CT + SM 
samples.fit <- samples.sm[[curr.set]]
samples.pred <- samples.ct.sm[[pred.set]]

access.fit <- c(covars.ct$Access[-ct.x.indices[[curr.set]]],
                covars.tr$Access[-tr.x.indices[[curr.set]]])
mean.access.fit <- mean(access.fit)
sd.access.fit <- sd(access.fit)
access.fit <- (access.fit - mean(access.fit)) / sd(access.fit)
access.pred <- c(covars.ct$Access[ct.x.indices[[curr.set]]],
                 covars.sm$Access[sm.x.indices[[curr.set]]])
access.pred <- (access.pred - mean.access.fit) / sd.access.fit 

tree.fit <- c(covars.ct$Tree[-ct.x.indices[[curr.set]]],
              covars.tr$Tree[-tr.x.indices[[curr.set]]])
mean.tree.fit <- mean(tree.fit)
sd.tree.fit <- sd(tree.fit)
tree.fit <- (tree.fit - mean(tree.fit)) / sd(tree.fit)
tree.pred <- c(covars.ct$Tree[ct.x.indices[[curr.set]]],
               covars.sm$Tree[sm.x.indices[[curr.set]]])
tree.pred <- (tree.pred - mean.tree.fit) / sd.tree.fit 

canopy.fit <- c(covars.ct$Canopy[-ct.x.indices[[curr.set]]],
                covars.tr$Canopy[-tr.x.indices[[curr.set]]])
mean.canopy.fit <- mean(canopy.fit)
sd.canopy.fit <- sd(canopy.fit)
canopy.fit <- (canopy.fit - mean(canopy.fit)) / sd(canopy.fit)
canopy.pred <- c(covars.ct$Canopy[ct.x.indices[[curr.set]]],
                 covars.sm$Canopy[sm.x.indices[[curr.set]]])
canopy.pred <- (canopy.pred - mean.canopy.fit) / sd.canopy.fit 

n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.ct <- length(ct.x.indices[[curr.set]])
J.sm <- length(sm.x.indices[[curr.set]])
access.ct.pred <- access.pred[1:J.ct]
access.sm.pred <- access.pred[1:J.sm + J.ct]
tree.ct.pred <- tree.pred[1:J.ct]
tree.sm.pred <- tree.pred[1:J.sm + J.ct]
canopy.ct.pred <- canopy.pred[1:J.ct]
canopy.sm.pred <- canopy.pred[1:J.sm + J.ct]

ct.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                       samples.pred = samples.pred, 
                       access.pred = access.ct.pred, 
                       tree.pred = tree.ct.pred, 
                       canopy.pred = canopy.ct.pred,
                       my.iter = my.iter, 
                       n.iter = n.iter, 
                       I = I, type = 'ct')
elpd.ct.vals[, 3, 4] <- apply(log(ct.elpd.vals.1), 1, sum)

sm.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                       samples.pred = samples.pred, 
                       access.pred = access.sm.pred, 
                       tree.pred = tree.sm.pred, 
                       canopy.pred = canopy.sm.pred,
                       my.iter = my.iter, 
                       n.iter = n.iter, 
                       I = I, type = 'sm')
elpd.sm.vals[, 2, 4] <- apply(log(sm.elpd.vals.1), 1, sum)

#### TR + SM 
samples.fit <- samples.sm[[curr.set]]
samples.pred <- samples.tr.sm[[pred.set]]

access.fit <- c(covars.ct$Access[-ct.x.indices[[curr.set]]],
                covars.tr$Access[-tr.x.indices[[curr.set]]])
mean.access.fit <- mean(access.fit)
sd.access.fit <- sd(access.fit)
access.fit <- (access.fit - mean(access.fit)) / sd(access.fit)
access.pred <- c(covars.tr$Access[tr.x.indices[[curr.set]]],
                 covars.sm$Access[sm.x.indices[[curr.set]]])
access.pred <- (access.pred - mean.access.fit) / sd.access.fit 

tree.fit <- c(covars.ct$Tree[-ct.x.indices[[curr.set]]],
              covars.tr$Tree[-tr.x.indices[[curr.set]]])
mean.tree.fit <- mean(tree.fit)
sd.tree.fit <- sd(tree.fit)
tree.fit <- (tree.fit - mean(tree.fit)) / sd(tree.fit)
tree.pred <- c(covars.tr$Tree[tr.x.indices[[curr.set]]],
               covars.sm$Tree[sm.x.indices[[curr.set]]])
tree.pred <- (tree.pred - mean.tree.fit) / sd.tree.fit 

canopy.fit <- c(covars.ct$Canopy[-ct.x.indices[[curr.set]]],
                covars.tr$Canopy[-tr.x.indices[[curr.set]]])
mean.canopy.fit <- mean(canopy.fit)
sd.canopy.fit <- sd(canopy.fit)
canopy.fit <- (canopy.fit - mean(canopy.fit)) / sd(canopy.fit)
canopy.pred <- c(covars.tr$Canopy[tr.x.indices[[curr.set]]],
                 covars.sm$Canopy[sm.x.indices[[curr.set]]])
canopy.pred <- (canopy.pred - mean.canopy.fit) / sd.canopy.fit 

n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.tr <- length(tr.x.indices[[curr.set]])
J.sm <- length(sm.x.indices[[curr.set]])
access.tr.pred <- access.pred[1:J.tr]
access.sm.pred <- access.pred[1:J.sm + J.tr]
tree.tr.pred <- tree.pred[1:J.tr]
tree.sm.pred <- tree.pred[1:J.sm + J.tr]
canopy.tr.pred <- canopy.pred[1:J.tr]
canopy.sm.pred <- canopy.pred[1:J.sm + J.tr]

tr.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                       samples.pred = samples.pred, 
                       access.pred = access.tr.pred, 
                       tree.pred = tree.tr.pred, 
                       canopy.pred = canopy.tr.pred,
                       my.iter = my.iter, 
                       n.iter = n.iter, 
                       I = I, type = 'tr')
elpd.tr.vals[, 3, 4] <- apply(log(tr.elpd.vals.1), 1, sum)

sm.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                       samples.pred = samples.pred, 
                       access.pred = access.sm.pred, 
                       tree.pred = tree.sm.pred, 
                       canopy.pred = canopy.sm.pred,
                       my.iter = my.iter, 
                       n.iter = n.iter, 
                       I = I, type = 'sm')
elpd.sm.vals[, 3, 4] <- apply(log(sm.elpd.vals.1), 1, sum)

#### CT + TR + SM 
samples.fit <- samples.sm[[curr.set]]
samples.pred <- samples.ct.tr.sm[[pred.set]]

access.fit <- c(covars.ct$Access[-ct.x.indices[[curr.set]]],
                covars.tr$Access[-tr.x.indices[[curr.set]]])
mean.access.fit <- mean(access.fit)
sd.access.fit <- sd(access.fit)
access.fit <- (access.fit - mean(access.fit)) / sd(access.fit)
access.pred <- c(covars.ct$Access[ct.x.indices[[curr.set]]],
                 covars.tr$Access[tr.x.indices[[curr.set]]],
                 covars.sm$Access[sm.x.indices[[curr.set]]])
access.pred <- (access.pred - mean.access.fit) / sd.access.fit 

tree.fit <- c(covars.ct$Tree[-ct.x.indices[[curr.set]]],
              covars.tr$Tree[-tr.x.indices[[curr.set]]])
mean.tree.fit <- mean(tree.fit)
sd.tree.fit <- sd(tree.fit)
tree.fit <- (tree.fit - mean(tree.fit)) / sd(tree.fit)
tree.pred <- c(covars.ct$Tree[ct.x.indices[[curr.set]]],
               covars.tr$Tree[tr.x.indices[[curr.set]]],
               covars.sm$Tree[sm.x.indices[[curr.set]]])
tree.pred <- (tree.pred - mean.tree.fit) / sd.tree.fit 

canopy.fit <- c(covars.ct$Canopy[-ct.x.indices[[curr.set]]],
                covars.tr$Canopy[-tr.x.indices[[curr.set]]])
mean.canopy.fit <- mean(canopy.fit)
sd.canopy.fit <- sd(canopy.fit)
canopy.fit <- (canopy.fit - mean(canopy.fit)) / sd(canopy.fit)
canopy.pred <- c(covars.ct$Canopy[ct.x.indices[[curr.set]]],
                 covars.tr$Canopy[tr.x.indices[[curr.set]]],
                 covars.sm$Canopy[sm.x.indices[[curr.set]]])
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
elpd.ct.vals[, 4, 4] <- apply(log(ct.elpd.vals.1), 1, sum)

tr.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                       samples.pred = samples.pred, 
                       access.pred = access.tr.pred, 
                       tree.pred = tree.tr.pred, 
                       canopy.pred = canopy.tr.pred,
                       my.iter = my.iter, 
                       n.iter = n.iter, 
                       I = I, type = 'tr')
elpd.tr.vals[, 4, 4] <- apply(log(tr.elpd.vals.1), 1, sum)

sm.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                       samples.pred = samples.pred, 
                       access.pred = access.sm.pred, 
                       tree.pred = tree.sm.pred, 
                       canopy.pred = canopy.sm.pred,
                       my.iter = my.iter, 
                       n.iter = n.iter, 
                       I = I, type = 'sm')
elpd.sm.vals[, 4, 4] <- apply(log(sm.elpd.vals.1), 1, sum)


#############
# CT-SM as predicting model -------------------------------------------------
##### CT 
samples.fit <- samples.sm[[curr.set]]
samples.pred <- samples.ct[[pred.set]]

access.fit <- c(covars.ct$Access[-ct.x.indices[[curr.set]]],
                covars.sm$Access[-sm.x.indices[[curr.set]]])
mean.access.fit <- mean(access.fit)
sd.access.fit <- sd(access.fit)
access.fit <- (access.fit - mean(access.fit)) / sd(access.fit)
access.pred <- c(covars.ct$Access[ct.x.indices[[curr.set]]])
access.pred <- (access.pred - mean.access.fit) / sd.access.fit 

tree.fit <- c(covars.ct$Tree[-ct.x.indices[[curr.set]]],
              covars.sm$Tree[-sm.x.indices[[curr.set]]])
mean.tree.fit <- mean(tree.fit)
sd.tree.fit <- sd(tree.fit)
tree.fit <- (tree.fit - mean(tree.fit)) / sd(tree.fit)
tree.pred <- c(covars.ct$Tree[ct.x.indices[[curr.set]]])
tree.pred <- (tree.pred - mean.tree.fit) / sd.tree.fit 

canopy.fit <- c(covars.ct$Canopy[-ct.x.indices[[curr.set]]],
                covars.sm$Canopy[-sm.x.indices[[curr.set]]])
mean.canopy.fit <- mean(canopy.fit)
sd.canopy.fit <- sd(canopy.fit)
canopy.fit <- (canopy.fit - mean(canopy.fit)) / sd(canopy.fit)
canopy.pred <- c(covars.ct$Canopy[ct.x.indices[[curr.set]]])
canopy.pred <- (canopy.pred - mean.canopy.fit) / sd.canopy.fit 

n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.ct <- length(ct.x.indices[[curr.set]])
elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                    samples.pred = samples.pred, 
                    access.pred = access.pred, 
                    tree.pred = tree.pred, 
                    canopy.pred = canopy.pred,
                    my.iter = my.iter, 
                    n.iter = n.iter, 
                    I = I, type = 'ct')
elpd.ct.vals[, 1, 5] <- apply(log(elpd.vals.1), 1, sum)

##### TR 
samples.fit <- samples.sm[[curr.set]]
samples.pred <- samples.tr[[pred.set]]

access.fit <- c(covars.ct$Access[-ct.x.indices[[curr.set]]],
                covars.sm$Access[-sm.x.indices[[curr.set]]])
mean.access.fit <- mean(access.fit)
sd.access.fit <- sd(access.fit)
access.fit <- (access.fit - mean(access.fit)) / sd(access.fit)
access.pred <- c(covars.tr$Access[tr.x.indices[[curr.set]]])
access.pred <- (access.pred - mean.access.fit) / sd.access.fit 

tree.fit <- c(covars.ct$Tree[-ct.x.indices[[curr.set]]],
              covars.sm$Tree[-sm.x.indices[[curr.set]]])
mean.tree.fit <- mean(tree.fit)
sd.tree.fit <- sd(tree.fit)
tree.fit <- (tree.fit - mean(tree.fit)) / sd(tree.fit)
tree.pred <- c(covars.tr$Tree[tr.x.indices[[curr.set]]])
tree.pred <- (tree.pred - mean.tree.fit) / sd.tree.fit 

canopy.fit <- c(covars.ct$Canopy[-ct.x.indices[[curr.set]]],
                covars.sm$Canopy[-sm.x.indices[[curr.set]]])
mean.canopy.fit <- mean(canopy.fit)
sd.canopy.fit <- sd(canopy.fit)
canopy.fit <- (canopy.fit - mean(canopy.fit)) / sd(canopy.fit)
canopy.pred <- c(covars.tr$Canopy[tr.x.indices[[curr.set]]])
canopy.pred <- (canopy.pred - mean.canopy.fit) / sd.canopy.fit 

n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.tr <- length(tr.x.indices[[curr.set]])
elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                    samples.pred = samples.pred, 
                    access.pred = access.pred, 
                    tree.pred = tree.pred, 
                    canopy.pred = canopy.pred,
                    my.iter = my.iter, 
                    n.iter = n.iter, 
                    I = I, type = 'tr')
elpd.tr.vals[, 1, 5] <- apply(log(elpd.vals.1), 1, sum)

##### SM 
samples.fit <- samples.sm[[curr.set]]
samples.pred <- samples.sm[[pred.set]]

access.fit <- c(covars.ct$Access[-ct.x.indices[[curr.set]]],
                covars.sm$Access[-sm.x.indices[[curr.set]]])
mean.access.fit <- mean(access.fit)
sd.access.fit <- sd(access.fit)
access.fit <- (access.fit - mean(access.fit)) / sd(access.fit)
access.pred <- c(covars.sm$Access[sm.x.indices[[curr.set]]])
access.pred <- (access.pred - mean.access.fit) / sd.access.fit 

tree.fit <- c(covars.ct$Tree[-ct.x.indices[[curr.set]]],
              covars.sm$Tree[-sm.x.indices[[curr.set]]])
mean.tree.fit <- mean(tree.fit)
sd.tree.fit <- sd(tree.fit)
tree.fit <- (tree.fit - mean(tree.fit)) / sd(tree.fit)
tree.pred <- c(covars.sm$Tree[sm.x.indices[[curr.set]]])
tree.pred <- (tree.pred - mean.tree.fit) / sd.tree.fit 

canopy.fit <- c(covars.ct$Canopy[-ct.x.indices[[curr.set]]],
                covars.sm$Canopy[-sm.x.indices[[curr.set]]])
mean.canopy.fit <- mean(canopy.fit)
sd.canopy.fit <- sd(canopy.fit)
canopy.fit <- (canopy.fit - mean(canopy.fit)) / sd(canopy.fit)
canopy.pred <- c(covars.sm$Canopy[sm.x.indices[[curr.set]]])
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
elpd.sm.vals[, 1, 5] <- apply(log(elpd.vals.1), 1, sum)

##### CT + TR 
samples.fit <- samples.sm[[curr.set]]
samples.pred <- samples.ct.tr[[pred.set]]

access.fit <- c(covars.ct$Access[-ct.x.indices[[curr.set]]],
                covars.sm$Access[-sm.x.indices[[curr.set]]])
mean.access.fit <- mean(access.fit)
sd.access.fit <- sd(access.fit)
access.fit <- (access.fit - mean(access.fit)) / sd(access.fit)
access.pred <- c(covars.ct$Access[ct.x.indices[[curr.set]]],
                 covars.tr$Access[tr.x.indices[[curr.set]]])
access.pred <- (access.pred - mean.access.fit) / sd.access.fit 

tree.fit <- c(covars.ct$Tree[-ct.x.indices[[curr.set]]],
              covars.sm$Tree[-sm.x.indices[[curr.set]]])
mean.tree.fit <- mean(tree.fit)
sd.tree.fit <- sd(tree.fit)
tree.fit <- (tree.fit - mean(tree.fit)) / sd(tree.fit)
tree.pred <- c(covars.ct$Tree[ct.x.indices[[curr.set]]],
               covars.tr$Tree[tr.x.indices[[curr.set]]])
tree.pred <- (tree.pred - mean.tree.fit) / sd.tree.fit 

canopy.fit <- c(covars.ct$Canopy[-ct.x.indices[[curr.set]]],
                covars.sm$Canopy[-sm.x.indices[[curr.set]]])
mean.canopy.fit <- mean(canopy.fit)
sd.canopy.fit <- sd(canopy.fit)
canopy.fit <- (canopy.fit - mean(canopy.fit)) / sd(canopy.fit)
canopy.pred <- c(covars.ct$Canopy[ct.x.indices[[curr.set]]],
                 covars.tr$Canopy[tr.x.indices[[curr.set]]])
canopy.pred <- (canopy.pred - mean.canopy.fit) / sd.canopy.fit 

n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.ct <- length(ct.x.indices[[curr.set]])
J.tr <- length(tr.x.indices[[curr.set]])
access.ct.pred <- access.pred[1:J.ct]
access.tr.pred <- access.pred[1:J.tr + J.ct]
tree.ct.pred <- tree.pred[1:J.ct]
tree.tr.pred <- tree.pred[1:J.tr + J.ct]
canopy.ct.pred <- canopy.pred[1:J.ct]
canopy.tr.pred <- canopy.pred[1:J.tr + J.ct]

ct.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                       samples.pred = samples.pred, 
                       access.pred = access.ct.pred, 
                       tree.pred = tree.ct.pred, 
                       canopy.pred = canopy.ct.pred,
                       my.iter = my.iter, 
                       n.iter = n.iter, 
                       I = I, type = 'ct')
elpd.ct.vals[, 2, 5] <- apply(log(ct.elpd.vals.1), 1, sum)

tr.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                       samples.pred = samples.pred, 
                       access.pred = access.tr.pred, 
                       tree.pred = tree.tr.pred, 
                       canopy.pred = canopy.tr.pred,
                       my.iter = my.iter, 
                       n.iter = n.iter, 
                       I = I, type = 'tr')
elpd.tr.vals[, 2, 5] <- apply(log(tr.elpd.vals.1), 1, sum)

##### CT + SM 
samples.fit <- samples.sm[[curr.set]]
samples.pred <- samples.ct.sm[[pred.set]]

access.fit <- c(covars.ct$Access[-ct.x.indices[[curr.set]]],
                covars.sm$Access[-sm.x.indices[[curr.set]]])
mean.access.fit <- mean(access.fit)
sd.access.fit <- sd(access.fit)
access.fit <- (access.fit - mean(access.fit)) / sd(access.fit)
access.pred <- c(covars.ct$Access[ct.x.indices[[curr.set]]],
                 covars.sm$Access[sm.x.indices[[curr.set]]])
access.pred <- (access.pred - mean.access.fit) / sd.access.fit 

tree.fit <- c(covars.ct$Tree[-ct.x.indices[[curr.set]]],
              covars.sm$Tree[-sm.x.indices[[curr.set]]])
mean.tree.fit <- mean(tree.fit)
sd.tree.fit <- sd(tree.fit)
tree.fit <- (tree.fit - mean(tree.fit)) / sd(tree.fit)
tree.pred <- c(covars.ct$Tree[ct.x.indices[[curr.set]]],
               covars.sm$Tree[sm.x.indices[[curr.set]]])
tree.pred <- (tree.pred - mean.tree.fit) / sd.tree.fit 

canopy.fit <- c(covars.ct$Canopy[-ct.x.indices[[curr.set]]],
                covars.sm$Canopy[-sm.x.indices[[curr.set]]])
mean.canopy.fit <- mean(canopy.fit)
sd.canopy.fit <- sd(canopy.fit)
canopy.fit <- (canopy.fit - mean(canopy.fit)) / sd(canopy.fit)
canopy.pred <- c(covars.ct$Canopy[ct.x.indices[[curr.set]]],
                 covars.sm$Canopy[sm.x.indices[[curr.set]]])
canopy.pred <- (canopy.pred - mean.canopy.fit) / sd.canopy.fit 

n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.ct <- length(ct.x.indices[[curr.set]])
J.sm <- length(sm.x.indices[[curr.set]])
access.ct.pred <- access.pred[1:J.ct]
access.sm.pred <- access.pred[1:J.sm + J.ct]
tree.ct.pred <- tree.pred[1:J.ct]
tree.sm.pred <- tree.pred[1:J.sm + J.ct]
canopy.ct.pred <- canopy.pred[1:J.ct]
canopy.sm.pred <- canopy.pred[1:J.sm + J.ct]

ct.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                       samples.pred = samples.pred, 
                       access.pred = access.ct.pred, 
                       tree.pred = tree.ct.pred, 
                       canopy.pred = canopy.ct.pred,
                       my.iter = my.iter, 
                       n.iter = n.iter, 
                       I = I, type = 'ct')
elpd.ct.vals[, 3, 5] <- apply(log(ct.elpd.vals.1), 1, sum)

sm.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                       samples.pred = samples.pred, 
                       access.pred = access.sm.pred, 
                       tree.pred = tree.sm.pred, 
                       canopy.pred = canopy.sm.pred,
                       my.iter = my.iter, 
                       n.iter = n.iter, 
                       I = I, type = 'sm')
elpd.sm.vals[, 2, 5] <- apply(log(sm.elpd.vals.1), 1, sum)

#### TR + SM 
samples.fit <- samples.sm[[curr.set]]
samples.pred <- samples.tr.sm[[pred.set]]

access.fit <- c(covars.ct$Access[-ct.x.indices[[curr.set]]],
                covars.sm$Access[-sm.x.indices[[curr.set]]])
mean.access.fit <- mean(access.fit)
sd.access.fit <- sd(access.fit)
access.fit <- (access.fit - mean(access.fit)) / sd(access.fit)
access.pred <- c(covars.tr$Access[tr.x.indices[[curr.set]]],
                 covars.sm$Access[sm.x.indices[[curr.set]]])
access.pred <- (access.pred - mean.access.fit) / sd.access.fit 

tree.fit <- c(covars.ct$Tree[-ct.x.indices[[curr.set]]],
              covars.sm$Tree[-sm.x.indices[[curr.set]]])
mean.tree.fit <- mean(tree.fit)
sd.tree.fit <- sd(tree.fit)
tree.fit <- (tree.fit - mean(tree.fit)) / sd(tree.fit)
tree.pred <- c(covars.tr$Tree[tr.x.indices[[curr.set]]],
               covars.sm$Tree[sm.x.indices[[curr.set]]])
tree.pred <- (tree.pred - mean.tree.fit) / sd.tree.fit 

canopy.fit <- c(covars.ct$Canopy[-ct.x.indices[[curr.set]]],
                covars.sm$Canopy[-sm.x.indices[[curr.set]]])
mean.canopy.fit <- mean(canopy.fit)
sd.canopy.fit <- sd(canopy.fit)
canopy.fit <- (canopy.fit - mean(canopy.fit)) / sd(canopy.fit)
canopy.pred <- c(covars.tr$Canopy[tr.x.indices[[curr.set]]],
                 covars.sm$Canopy[sm.x.indices[[curr.set]]])
canopy.pred <- (canopy.pred - mean.canopy.fit) / sd.canopy.fit 

n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.tr <- length(tr.x.indices[[curr.set]])
J.sm <- length(sm.x.indices[[curr.set]])
access.tr.pred <- access.pred[1:J.tr]
access.sm.pred <- access.pred[1:J.sm + J.tr]
tree.tr.pred <- tree.pred[1:J.tr]
tree.sm.pred <- tree.pred[1:J.sm + J.tr]
canopy.tr.pred <- canopy.pred[1:J.tr]
canopy.sm.pred <- canopy.pred[1:J.sm + J.tr]

tr.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                       samples.pred = samples.pred, 
                       access.pred = access.tr.pred, 
                       tree.pred = tree.tr.pred, 
                       canopy.pred = canopy.tr.pred,
                       my.iter = my.iter, 
                       n.iter = n.iter, 
                       I = I, type = 'tr')
elpd.tr.vals[, 3, 5] <- apply(log(tr.elpd.vals.1), 1, sum)

sm.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                       samples.pred = samples.pred, 
                       access.pred = access.sm.pred, 
                       tree.pred = tree.sm.pred, 
                       canopy.pred = canopy.sm.pred,
                       my.iter = my.iter, 
                       n.iter = n.iter, 
                       I = I, type = 'sm')
elpd.sm.vals[, 3, 5] <- apply(log(sm.elpd.vals.1), 1, sum)

#### CT + TR + SM 
samples.fit <- samples.sm[[curr.set]]
samples.pred <- samples.ct.tr.sm[[pred.set]]

access.fit <- c(covars.ct$Access[-ct.x.indices[[curr.set]]],
                covars.sm$Access[-sm.x.indices[[curr.set]]])
mean.access.fit <- mean(access.fit)
sd.access.fit <- sd(access.fit)
access.fit <- (access.fit - mean(access.fit)) / sd(access.fit)
access.pred <- c(covars.ct$Access[ct.x.indices[[curr.set]]],
                 covars.tr$Access[tr.x.indices[[curr.set]]],
                 covars.sm$Access[sm.x.indices[[curr.set]]])
access.pred <- (access.pred - mean.access.fit) / sd.access.fit 

tree.fit <- c(covars.ct$Tree[-ct.x.indices[[curr.set]]],
              covars.sm$Tree[-sm.x.indices[[curr.set]]])
mean.tree.fit <- mean(tree.fit)
sd.tree.fit <- sd(tree.fit)
tree.fit <- (tree.fit - mean(tree.fit)) / sd(tree.fit)
tree.pred <- c(covars.ct$Tree[ct.x.indices[[curr.set]]],
               covars.tr$Tree[tr.x.indices[[curr.set]]],
               covars.sm$Tree[sm.x.indices[[curr.set]]])
tree.pred <- (tree.pred - mean.tree.fit) / sd.tree.fit 

canopy.fit <- c(covars.ct$Canopy[-ct.x.indices[[curr.set]]],
                covars.sm$Canopy[-sm.x.indices[[curr.set]]])
mean.canopy.fit <- mean(canopy.fit)
sd.canopy.fit <- sd(canopy.fit)
canopy.fit <- (canopy.fit - mean(canopy.fit)) / sd(canopy.fit)
canopy.pred <- c(covars.ct$Canopy[ct.x.indices[[curr.set]]],
                 covars.tr$Canopy[tr.x.indices[[curr.set]]],
                 covars.sm$Canopy[sm.x.indices[[curr.set]]])
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
elpd.ct.vals[, 4, 5] <- apply(log(ct.elpd.vals.1), 1, sum)

tr.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                       samples.pred = samples.pred, 
                       access.pred = access.tr.pred, 
                       tree.pred = tree.tr.pred, 
                       canopy.pred = canopy.tr.pred,
                       my.iter = my.iter, 
                       n.iter = n.iter, 
                       I = I, type = 'tr')
elpd.tr.vals[, 4, 5] <- apply(log(tr.elpd.vals.1), 1, sum)

sm.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                       samples.pred = samples.pred, 
                       access.pred = access.sm.pred, 
                       tree.pred = tree.sm.pred, 
                       canopy.pred = canopy.sm.pred,
                       my.iter = my.iter, 
                       n.iter = n.iter, 
                       I = I, type = 'sm')
elpd.sm.vals[, 4, 5] <- apply(log(sm.elpd.vals.1), 1, sum)


#############
# TR-SM as predicting model -------------------------------------------------
##### CT 
samples.fit <- samples.sm[[curr.set]]
samples.pred <- samples.ct[[pred.set]]

access.fit <- c(covars.tr$Access[-tr.x.indices[[curr.set]]],
                covars.sm$Access[-sm.x.indices[[curr.set]]])
mean.access.fit <- mean(access.fit)
sd.access.fit <- sd(access.fit)
access.fit <- (access.fit - mean(access.fit)) / sd(access.fit)
access.pred <- c(covars.ct$Access[ct.x.indices[[curr.set]]])
access.pred <- (access.pred - mean.access.fit) / sd.access.fit 

tree.fit <- c(covars.tr$Tree[-tr.x.indices[[curr.set]]],
              covars.sm$Tree[-sm.x.indices[[curr.set]]])
mean.tree.fit <- mean(tree.fit)
sd.tree.fit <- sd(tree.fit)
tree.fit <- (tree.fit - mean(tree.fit)) / sd(tree.fit)
tree.pred <- c(covars.ct$Tree[ct.x.indices[[curr.set]]])
tree.pred <- (tree.pred - mean.tree.fit) / sd.tree.fit 

canopy.fit <- c(covars.tr$Canopy[-tr.x.indices[[curr.set]]],
                covars.sm$Canopy[-sm.x.indices[[curr.set]]])
mean.canopy.fit <- mean(canopy.fit)
sd.canopy.fit <- sd(canopy.fit)
canopy.fit <- (canopy.fit - mean(canopy.fit)) / sd(canopy.fit)
canopy.pred <- c(covars.ct$Canopy[ct.x.indices[[curr.set]]])
canopy.pred <- (canopy.pred - mean.canopy.fit) / sd.canopy.fit 

n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.ct <- length(ct.x.indices[[curr.set]])
elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                    samples.pred = samples.pred, 
                    access.pred = access.pred, 
                    tree.pred = tree.pred, 
                    canopy.pred = canopy.pred,
                    my.iter = my.iter, 
                    n.iter = n.iter, 
                    I = I, type = 'ct')
elpd.ct.vals[, 1, 6] <- apply(log(elpd.vals.1), 1, sum)

##### TR 
samples.fit <- samples.sm[[curr.set]]
samples.pred <- samples.tr[[pred.set]]

access.fit <- c(covars.tr$Access[-tr.x.indices[[curr.set]]],
                covars.sm$Access[-sm.x.indices[[curr.set]]])
mean.access.fit <- mean(access.fit)
sd.access.fit <- sd(access.fit)
access.fit <- (access.fit - mean(access.fit)) / sd(access.fit)
access.pred <- c(covars.tr$Access[tr.x.indices[[curr.set]]])
access.pred <- (access.pred - mean.access.fit) / sd.access.fit 

tree.fit <- c(covars.tr$Tree[-tr.x.indices[[curr.set]]],
              covars.sm$Tree[-sm.x.indices[[curr.set]]])
mean.tree.fit <- mean(tree.fit)
sd.tree.fit <- sd(tree.fit)
tree.fit <- (tree.fit - mean(tree.fit)) / sd(tree.fit)
tree.pred <- c(covars.tr$Tree[tr.x.indices[[curr.set]]])
tree.pred <- (tree.pred - mean.tree.fit) / sd.tree.fit 

canopy.fit <- c(covars.tr$Canopy[-tr.x.indices[[curr.set]]],
                covars.sm$Canopy[-sm.x.indices[[curr.set]]])
mean.canopy.fit <- mean(canopy.fit)
sd.canopy.fit <- sd(canopy.fit)
canopy.fit <- (canopy.fit - mean(canopy.fit)) / sd(canopy.fit)
canopy.pred <- c(covars.tr$Canopy[tr.x.indices[[curr.set]]])
canopy.pred <- (canopy.pred - mean.canopy.fit) / sd.canopy.fit 

n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.tr <- length(tr.x.indices[[curr.set]])
elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                    samples.pred = samples.pred, 
                    access.pred = access.pred, 
                    tree.pred = tree.pred, 
                    canopy.pred = canopy.pred,
                    my.iter = my.iter, 
                    n.iter = n.iter, 
                    I = I, type = 'tr')
elpd.tr.vals[, 1, 6] <- apply(log(elpd.vals.1), 1, sum)

##### SM 
samples.fit <- samples.sm[[curr.set]]
samples.pred <- samples.sm[[pred.set]]

access.fit <- c(covars.tr$Access[-tr.x.indices[[curr.set]]],
                covars.sm$Access[-sm.x.indices[[curr.set]]])
mean.access.fit <- mean(access.fit)
sd.access.fit <- sd(access.fit)
access.fit <- (access.fit - mean(access.fit)) / sd(access.fit)
access.pred <- c(covars.sm$Access[sm.x.indices[[curr.set]]])
access.pred <- (access.pred - mean.access.fit) / sd.access.fit 

tree.fit <- c(covars.tr$Tree[-tr.x.indices[[curr.set]]],
              covars.sm$Tree[-sm.x.indices[[curr.set]]])
mean.tree.fit <- mean(tree.fit)
sd.tree.fit <- sd(tree.fit)
tree.fit <- (tree.fit - mean(tree.fit)) / sd(tree.fit)
tree.pred <- c(covars.sm$Tree[sm.x.indices[[curr.set]]])
tree.pred <- (tree.pred - mean.tree.fit) / sd.tree.fit 

canopy.fit <- c(covars.tr$Canopy[-tr.x.indices[[curr.set]]],
                covars.sm$Canopy[-sm.x.indices[[curr.set]]])
mean.canopy.fit <- mean(canopy.fit)
sd.canopy.fit <- sd(canopy.fit)
canopy.fit <- (canopy.fit - mean(canopy.fit)) / sd(canopy.fit)
canopy.pred <- c(covars.sm$Canopy[sm.x.indices[[curr.set]]])
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
elpd.sm.vals[, 1, 6] <- apply(log(elpd.vals.1), 1, sum)

##### CT + TR 
samples.fit <- samples.sm[[curr.set]]
samples.pred <- samples.ct.tr[[pred.set]]

access.fit <- c(covars.tr$Access[-tr.x.indices[[curr.set]]],
                covars.sm$Access[-sm.x.indices[[curr.set]]])
mean.access.fit <- mean(access.fit)
sd.access.fit <- sd(access.fit)
access.fit <- (access.fit - mean(access.fit)) / sd(access.fit)
access.pred <- c(covars.ct$Access[ct.x.indices[[curr.set]]],
                 covars.tr$Access[tr.x.indices[[curr.set]]])
access.pred <- (access.pred - mean.access.fit) / sd.access.fit 

tree.fit <- c(covars.tr$Tree[-tr.x.indices[[curr.set]]],
              covars.sm$Tree[-sm.x.indices[[curr.set]]])
mean.tree.fit <- mean(tree.fit)
sd.tree.fit <- sd(tree.fit)
tree.fit <- (tree.fit - mean(tree.fit)) / sd(tree.fit)
tree.pred <- c(covars.ct$Tree[ct.x.indices[[curr.set]]],
               covars.tr$Tree[tr.x.indices[[curr.set]]])
tree.pred <- (tree.pred - mean.tree.fit) / sd.tree.fit 

canopy.fit <- c(covars.tr$Canopy[-tr.x.indices[[curr.set]]],
                covars.sm$Canopy[-sm.x.indices[[curr.set]]])
mean.canopy.fit <- mean(canopy.fit)
sd.canopy.fit <- sd(canopy.fit)
canopy.fit <- (canopy.fit - mean(canopy.fit)) / sd(canopy.fit)
canopy.pred <- c(covars.ct$Canopy[ct.x.indices[[curr.set]]],
                 covars.tr$Canopy[tr.x.indices[[curr.set]]])
canopy.pred <- (canopy.pred - mean.canopy.fit) / sd.canopy.fit 

n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.ct <- length(ct.x.indices[[curr.set]])
J.tr <- length(tr.x.indices[[curr.set]])
access.ct.pred <- access.pred[1:J.ct]
access.tr.pred <- access.pred[1:J.tr + J.ct]
tree.ct.pred <- tree.pred[1:J.ct]
tree.tr.pred <- tree.pred[1:J.tr + J.ct]
canopy.ct.pred <- canopy.pred[1:J.ct]
canopy.tr.pred <- canopy.pred[1:J.tr + J.ct]

ct.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                       samples.pred = samples.pred, 
                       access.pred = access.ct.pred, 
                       tree.pred = tree.ct.pred, 
                       canopy.pred = canopy.ct.pred,
                       my.iter = my.iter, 
                       n.iter = n.iter, 
                       I = I, type = 'ct')
elpd.ct.vals[, 2, 6] <- apply(log(ct.elpd.vals.1), 1, sum)

tr.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                       samples.pred = samples.pred, 
                       access.pred = access.tr.pred, 
                       tree.pred = tree.tr.pred, 
                       canopy.pred = canopy.tr.pred,
                       my.iter = my.iter, 
                       n.iter = n.iter, 
                       I = I, type = 'tr')
elpd.tr.vals[, 2, 6] <- apply(log(tr.elpd.vals.1), 1, sum)

##### CT + SM 
samples.fit <- samples.sm[[curr.set]]
samples.pred <- samples.ct.sm[[pred.set]]

access.fit <- c(covars.tr$Access[-tr.x.indices[[curr.set]]],
                covars.sm$Access[-sm.x.indices[[curr.set]]])
mean.access.fit <- mean(access.fit)
sd.access.fit <- sd(access.fit)
access.fit <- (access.fit - mean(access.fit)) / sd(access.fit)
access.pred <- c(covars.ct$Access[ct.x.indices[[curr.set]]],
                 covars.sm$Access[sm.x.indices[[curr.set]]])
access.pred <- (access.pred - mean.access.fit) / sd.access.fit 

tree.fit <- c(covars.tr$Tree[-tr.x.indices[[curr.set]]],
              covars.sm$Tree[-sm.x.indices[[curr.set]]])
mean.tree.fit <- mean(tree.fit)
sd.tree.fit <- sd(tree.fit)
tree.fit <- (tree.fit - mean(tree.fit)) / sd(tree.fit)
tree.pred <- c(covars.ct$Tree[ct.x.indices[[curr.set]]],
               covars.sm$Tree[sm.x.indices[[curr.set]]])
tree.pred <- (tree.pred - mean.tree.fit) / sd.tree.fit 

canopy.fit <- c(covars.tr$Canopy[-tr.x.indices[[curr.set]]],
                covars.sm$Canopy[-sm.x.indices[[curr.set]]])
mean.canopy.fit <- mean(canopy.fit)
sd.canopy.fit <- sd(canopy.fit)
canopy.fit <- (canopy.fit - mean(canopy.fit)) / sd(canopy.fit)
canopy.pred <- c(covars.ct$Canopy[ct.x.indices[[curr.set]]],
                 covars.sm$Canopy[sm.x.indices[[curr.set]]])
canopy.pred <- (canopy.pred - mean.canopy.fit) / sd.canopy.fit 

n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.ct <- length(ct.x.indices[[curr.set]])
J.sm <- length(sm.x.indices[[curr.set]])
access.ct.pred <- access.pred[1:J.ct]
access.sm.pred <- access.pred[1:J.sm + J.ct]
tree.ct.pred <- tree.pred[1:J.ct]
tree.sm.pred <- tree.pred[1:J.sm + J.ct]
canopy.ct.pred <- canopy.pred[1:J.ct]
canopy.sm.pred <- canopy.pred[1:J.sm + J.ct]

ct.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                       samples.pred = samples.pred, 
                       access.pred = access.ct.pred, 
                       tree.pred = tree.ct.pred, 
                       canopy.pred = canopy.ct.pred,
                       my.iter = my.iter, 
                       n.iter = n.iter, 
                       I = I, type = 'ct')
elpd.ct.vals[, 3, 6] <- apply(log(ct.elpd.vals.1), 1, sum)

sm.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                       samples.pred = samples.pred, 
                       access.pred = access.sm.pred, 
                       tree.pred = tree.sm.pred, 
                       canopy.pred = canopy.sm.pred,
                       my.iter = my.iter, 
                       n.iter = n.iter, 
                       I = I, type = 'sm')
elpd.sm.vals[, 2, 6] <- apply(log(sm.elpd.vals.1), 1, sum)

#### TR + SM 
samples.fit <- samples.sm[[curr.set]]
samples.pred <- samples.tr.sm[[pred.set]]

access.fit <- c(covars.tr$Access[-tr.x.indices[[curr.set]]],
                covars.sm$Access[-sm.x.indices[[curr.set]]])
mean.access.fit <- mean(access.fit)
sd.access.fit <- sd(access.fit)
access.fit <- (access.fit - mean(access.fit)) / sd(access.fit)
access.pred <- c(covars.tr$Access[tr.x.indices[[curr.set]]],
                 covars.sm$Access[sm.x.indices[[curr.set]]])
access.pred <- (access.pred - mean.access.fit) / sd.access.fit 

tree.fit <- c(covars.tr$Tree[-tr.x.indices[[curr.set]]],
              covars.sm$Tree[-sm.x.indices[[curr.set]]])
mean.tree.fit <- mean(tree.fit)
sd.tree.fit <- sd(tree.fit)
tree.fit <- (tree.fit - mean(tree.fit)) / sd(tree.fit)
tree.pred <- c(covars.tr$Tree[tr.x.indices[[curr.set]]],
               covars.sm$Tree[sm.x.indices[[curr.set]]])
tree.pred <- (tree.pred - mean.tree.fit) / sd.tree.fit 

canopy.fit <- c(covars.tr$Canopy[-tr.x.indices[[curr.set]]],
                covars.sm$Canopy[-sm.x.indices[[curr.set]]])
mean.canopy.fit <- mean(canopy.fit)
sd.canopy.fit <- sd(canopy.fit)
canopy.fit <- (canopy.fit - mean(canopy.fit)) / sd(canopy.fit)
canopy.pred <- c(covars.tr$Canopy[tr.x.indices[[curr.set]]],
                 covars.sm$Canopy[sm.x.indices[[curr.set]]])
canopy.pred <- (canopy.pred - mean.canopy.fit) / sd.canopy.fit 

n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.tr <- length(tr.x.indices[[curr.set]])
J.sm <- length(sm.x.indices[[curr.set]])
access.tr.pred <- access.pred[1:J.tr]
access.sm.pred <- access.pred[1:J.sm + J.tr]
tree.tr.pred <- tree.pred[1:J.tr]
tree.sm.pred <- tree.pred[1:J.sm + J.tr]
canopy.tr.pred <- canopy.pred[1:J.tr]
canopy.sm.pred <- canopy.pred[1:J.sm + J.tr]

tr.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                       samples.pred = samples.pred, 
                       access.pred = access.tr.pred, 
                       tree.pred = tree.tr.pred, 
                       canopy.pred = canopy.tr.pred,
                       my.iter = my.iter, 
                       n.iter = n.iter, 
                       I = I, type = 'tr')
elpd.tr.vals[, 3, 6] <- apply(log(tr.elpd.vals.1), 1, sum)

sm.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                       samples.pred = samples.pred, 
                       access.pred = access.sm.pred, 
                       tree.pred = tree.sm.pred, 
                       canopy.pred = canopy.sm.pred,
                       my.iter = my.iter, 
                       n.iter = n.iter, 
                       I = I, type = 'sm')
elpd.sm.vals[, 3, 6] <- apply(log(sm.elpd.vals.1), 1, sum)

#### CT + TR + SM 
samples.fit <- samples.sm[[curr.set]]
samples.pred <- samples.ct.tr.sm[[pred.set]]

access.fit <- c(covars.tr$Access[-tr.x.indices[[curr.set]]],
                covars.sm$Access[-sm.x.indices[[curr.set]]])
mean.access.fit <- mean(access.fit)
sd.access.fit <- sd(access.fit)
access.fit <- (access.fit - mean(access.fit)) / sd(access.fit)
access.pred <- c(covars.ct$Access[ct.x.indices[[curr.set]]],
                 covars.tr$Access[tr.x.indices[[curr.set]]],
                 covars.sm$Access[sm.x.indices[[curr.set]]])
access.pred <- (access.pred - mean.access.fit) / sd.access.fit 

tree.fit <- c(covars.tr$Tree[-tr.x.indices[[curr.set]]],
              covars.sm$Tree[-sm.x.indices[[curr.set]]])
mean.tree.fit <- mean(tree.fit)
sd.tree.fit <- sd(tree.fit)
tree.fit <- (tree.fit - mean(tree.fit)) / sd(tree.fit)
tree.pred <- c(covars.ct$Tree[ct.x.indices[[curr.set]]],
               covars.tr$Tree[tr.x.indices[[curr.set]]],
               covars.sm$Tree[sm.x.indices[[curr.set]]])
tree.pred <- (tree.pred - mean.tree.fit) / sd.tree.fit 

canopy.fit <- c(covars.tr$Canopy[-tr.x.indices[[curr.set]]],
                covars.sm$Canopy[-sm.x.indices[[curr.set]]])
mean.canopy.fit <- mean(canopy.fit)
sd.canopy.fit <- sd(canopy.fit)
canopy.fit <- (canopy.fit - mean(canopy.fit)) / sd(canopy.fit)
canopy.pred <- c(covars.ct$Canopy[ct.x.indices[[curr.set]]],
                 covars.tr$Canopy[tr.x.indices[[curr.set]]],
                 covars.sm$Canopy[sm.x.indices[[curr.set]]])
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
elpd.ct.vals[, 4, 6] <- apply(log(ct.elpd.vals.1), 1, sum)

tr.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                       samples.pred = samples.pred, 
                       access.pred = access.tr.pred, 
                       tree.pred = tree.tr.pred, 
                       canopy.pred = canopy.tr.pred,
                       my.iter = my.iter, 
                       n.iter = n.iter, 
                       I = I, type = 'tr')
elpd.tr.vals[, 4, 6] <- apply(log(tr.elpd.vals.1), 1, sum)

sm.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                       samples.pred = samples.pred, 
                       access.pred = access.sm.pred, 
                       tree.pred = tree.sm.pred, 
                       canopy.pred = canopy.sm.pred,
                       my.iter = my.iter, 
                       n.iter = n.iter, 
                       I = I, type = 'sm')
elpd.sm.vals[, 4, 6] <- apply(log(sm.elpd.vals.1), 1, sum)


#############
# CT-TR-SM as predicting model -------------------------------------------------
##### CT 
samples.fit <- samples.sm[[curr.set]]
samples.pred <- samples.ct[[pred.set]]

access.fit <- c(covars.ct$Access[-ct.x.indices[[curr.set]]],
                covars.tr$Access[-tr.x.indices[[curr.set]]],
                covars.sm$Access[-sm.x.indices[[curr.set]]])
mean.access.fit <- mean(access.fit)
sd.access.fit <- sd(access.fit)
access.fit <- (access.fit - mean(access.fit)) / sd(access.fit)
access.pred <- c(covars.ct$Access[ct.x.indices[[curr.set]]])
access.pred <- (access.pred - mean.access.fit) / sd.access.fit 

tree.fit <- c(covars.ct$Tree[-ct.x.indices[[curr.set]]],
              covars.tr$Tree[-tr.x.indices[[curr.set]]],
              covars.sm$Tree[-sm.x.indices[[curr.set]]])
mean.tree.fit <- mean(tree.fit)
sd.tree.fit <- sd(tree.fit)
tree.fit <- (tree.fit - mean(tree.fit)) / sd(tree.fit)
tree.pred <- c(covars.ct$Tree[ct.x.indices[[curr.set]]])
tree.pred <- (tree.pred - mean.tree.fit) / sd.tree.fit 

canopy.fit <- c(covars.ct$Canopy[-ct.x.indices[[curr.set]]],
                covars.tr$Canopy[-tr.x.indices[[curr.set]]],
                covars.sm$Canopy[-sm.x.indices[[curr.set]]])
mean.canopy.fit <- mean(canopy.fit)
sd.canopy.fit <- sd(canopy.fit)
canopy.fit <- (canopy.fit - mean(canopy.fit)) / sd(canopy.fit)
canopy.pred <- c(covars.ct$Canopy[ct.x.indices[[curr.set]]])
canopy.pred <- (canopy.pred - mean.canopy.fit) / sd.canopy.fit 

n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.ct <- length(ct.x.indices[[curr.set]])
elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                    samples.pred = samples.pred, 
                    access.pred = access.pred, 
                    tree.pred = tree.pred, 
                    canopy.pred = canopy.pred,
                    my.iter = my.iter, 
                    n.iter = n.iter, 
                    I = I, type = 'ct')
elpd.ct.vals[, 1, 7] <- apply(log(elpd.vals.1), 1, sum)

##### TR 
samples.fit <- samples.sm[[curr.set]]
samples.pred <- samples.tr[[pred.set]]

access.fit <- c(covars.ct$Access[-ct.x.indices[[curr.set]]],
                covars.tr$Access[-tr.x.indices[[curr.set]]],
                covars.sm$Access[-sm.x.indices[[curr.set]]])
mean.access.fit <- mean(access.fit)
sd.access.fit <- sd(access.fit)
access.fit <- (access.fit - mean(access.fit)) / sd(access.fit)
access.pred <- c(covars.tr$Access[tr.x.indices[[curr.set]]])
access.pred <- (access.pred - mean.access.fit) / sd.access.fit 

tree.fit <- c(covars.ct$Tree[-ct.x.indices[[curr.set]]],
              covars.tr$Tree[-tr.x.indices[[curr.set]]],
              covars.sm$Tree[-sm.x.indices[[curr.set]]])
mean.tree.fit <- mean(tree.fit)
sd.tree.fit <- sd(tree.fit)
tree.fit <- (tree.fit - mean(tree.fit)) / sd(tree.fit)
tree.pred <- c(covars.tr$Tree[tr.x.indices[[curr.set]]])
tree.pred <- (tree.pred - mean.tree.fit) / sd.tree.fit 

canopy.fit <- c(covars.ct$Canopy[-ct.x.indices[[curr.set]]],
                covars.tr$Canopy[-tr.x.indices[[curr.set]]],
                covars.sm$Canopy[-sm.x.indices[[curr.set]]])
mean.canopy.fit <- mean(canopy.fit)
sd.canopy.fit <- sd(canopy.fit)
canopy.fit <- (canopy.fit - mean(canopy.fit)) / sd(canopy.fit)
canopy.pred <- c(covars.tr$Canopy[tr.x.indices[[curr.set]]])
canopy.pred <- (canopy.pred - mean.canopy.fit) / sd.canopy.fit 

n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.tr <- length(tr.x.indices[[curr.set]])
elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                    samples.pred = samples.pred, 
                    access.pred = access.pred, 
                    tree.pred = tree.pred, 
                    canopy.pred = canopy.pred,
                    my.iter = my.iter, 
                    n.iter = n.iter, 
                    I = I, type = 'tr')
elpd.tr.vals[, 1, 7] <- apply(log(elpd.vals.1), 1, sum)

##### SM 
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

canopy.fit <- c(covars.ct$Canopy[-ct.x.indices[[curr.set]]],
                covars.tr$Canopy[-tr.x.indices[[curr.set]]],
                covars.sm$Canopy[-sm.x.indices[[curr.set]]])
mean.canopy.fit <- mean(canopy.fit)
sd.canopy.fit <- sd(canopy.fit)
canopy.fit <- (canopy.fit - mean(canopy.fit)) / sd(canopy.fit)
canopy.pred <- c(covars.sm$Canopy[sm.x.indices[[curr.set]]])
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
elpd.sm.vals[, 1, 7] <- apply(log(elpd.vals.1), 1, sum)

##### CT + TR 
samples.fit <- samples.sm[[curr.set]]
samples.pred <- samples.ct.tr[[pred.set]]

access.fit <- c(covars.ct$Access[-ct.x.indices[[curr.set]]],
                covars.tr$Access[-tr.x.indices[[curr.set]]],
                covars.sm$Access[-sm.x.indices[[curr.set]]])
mean.access.fit <- mean(access.fit)
sd.access.fit <- sd(access.fit)
access.fit <- (access.fit - mean(access.fit)) / sd(access.fit)
access.pred <- c(covars.ct$Access[ct.x.indices[[curr.set]]],
                 covars.tr$Access[tr.x.indices[[curr.set]]])
access.pred <- (access.pred - mean.access.fit) / sd.access.fit 

tree.fit <- c(covars.ct$Tree[-ct.x.indices[[curr.set]]],
              covars.tr$Tree[-tr.x.indices[[curr.set]]],
              covars.sm$Tree[-sm.x.indices[[curr.set]]])
mean.tree.fit <- mean(tree.fit)
sd.tree.fit <- sd(tree.fit)
tree.fit <- (tree.fit - mean(tree.fit)) / sd(tree.fit)
tree.pred <- c(covars.ct$Tree[ct.x.indices[[curr.set]]],
               covars.tr$Tree[tr.x.indices[[curr.set]]])
tree.pred <- (tree.pred - mean.tree.fit) / sd.tree.fit 

canopy.fit <- c(covars.ct$Canopy[-ct.x.indices[[curr.set]]],
                covars.tr$Canopy[-tr.x.indices[[curr.set]]],
                covars.sm$Canopy[-sm.x.indices[[curr.set]]])
mean.canopy.fit <- mean(canopy.fit)
sd.canopy.fit <- sd(canopy.fit)
canopy.fit <- (canopy.fit - mean(canopy.fit)) / sd(canopy.fit)
canopy.pred <- c(covars.ct$Canopy[ct.x.indices[[curr.set]]],
                 covars.tr$Canopy[tr.x.indices[[curr.set]]])
canopy.pred <- (canopy.pred - mean.canopy.fit) / sd.canopy.fit 

n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.ct <- length(ct.x.indices[[curr.set]])
J.tr <- length(tr.x.indices[[curr.set]])
access.ct.pred <- access.pred[1:J.ct]
access.tr.pred <- access.pred[1:J.tr + J.ct]
tree.ct.pred <- tree.pred[1:J.ct]
tree.tr.pred <- tree.pred[1:J.tr + J.ct]
canopy.ct.pred <- canopy.pred[1:J.ct]
canopy.tr.pred <- canopy.pred[1:J.tr + J.ct]

ct.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                       samples.pred = samples.pred, 
                       access.pred = access.ct.pred, 
                       tree.pred = tree.ct.pred, 
                       canopy.pred = canopy.ct.pred,
                       my.iter = my.iter, 
                       n.iter = n.iter, 
                       I = I, type = 'ct')
elpd.ct.vals[, 2, 7] <- apply(log(ct.elpd.vals.1), 1, sum)

tr.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                       samples.pred = samples.pred, 
                       access.pred = access.tr.pred, 
                       tree.pred = tree.tr.pred, 
                       canopy.pred = canopy.tr.pred,
                       my.iter = my.iter, 
                       n.iter = n.iter, 
                       I = I, type = 'tr')
elpd.tr.vals[, 2, 7] <- apply(log(tr.elpd.vals.1), 1, sum)

##### CT + SM 
samples.fit <- samples.sm[[curr.set]]
samples.pred <- samples.ct.sm[[pred.set]]

access.fit <- c(covars.ct$Access[-ct.x.indices[[curr.set]]],
                covars.tr$Access[-tr.x.indices[[curr.set]]],
                covars.sm$Access[-sm.x.indices[[curr.set]]])
mean.access.fit <- mean(access.fit)
sd.access.fit <- sd(access.fit)
access.fit <- (access.fit - mean(access.fit)) / sd(access.fit)
access.pred <- c(covars.ct$Access[ct.x.indices[[curr.set]]],
                 covars.sm$Access[sm.x.indices[[curr.set]]])
access.pred <- (access.pred - mean.access.fit) / sd.access.fit 

tree.fit <- c(covars.ct$Tree[-ct.x.indices[[curr.set]]],
              covars.tr$Tree[-tr.x.indices[[curr.set]]],
              covars.sm$Tree[-sm.x.indices[[curr.set]]])
mean.tree.fit <- mean(tree.fit)
sd.tree.fit <- sd(tree.fit)
tree.fit <- (tree.fit - mean(tree.fit)) / sd(tree.fit)
tree.pred <- c(covars.ct$Tree[ct.x.indices[[curr.set]]],
               covars.sm$Tree[sm.x.indices[[curr.set]]])
tree.pred <- (tree.pred - mean.tree.fit) / sd.tree.fit 

canopy.fit <- c(covars.ct$Canopy[-ct.x.indices[[curr.set]]],
                covars.tr$Canopy[-tr.x.indices[[curr.set]]],
                covars.sm$Canopy[-sm.x.indices[[curr.set]]])
mean.canopy.fit <- mean(canopy.fit)
sd.canopy.fit <- sd(canopy.fit)
canopy.fit <- (canopy.fit - mean(canopy.fit)) / sd(canopy.fit)
canopy.pred <- c(covars.ct$Canopy[ct.x.indices[[curr.set]]],
                 covars.sm$Canopy[sm.x.indices[[curr.set]]])
canopy.pred <- (canopy.pred - mean.canopy.fit) / sd.canopy.fit 

n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.ct <- length(ct.x.indices[[curr.set]])
J.sm <- length(sm.x.indices[[curr.set]])
access.ct.pred <- access.pred[1:J.ct]
access.sm.pred <- access.pred[1:J.sm + J.ct]
tree.ct.pred <- tree.pred[1:J.ct]
tree.sm.pred <- tree.pred[1:J.sm + J.ct]
canopy.ct.pred <- canopy.pred[1:J.ct]
canopy.sm.pred <- canopy.pred[1:J.sm + J.ct]

ct.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                       samples.pred = samples.pred, 
                       access.pred = access.ct.pred, 
                       tree.pred = tree.ct.pred, 
                       canopy.pred = canopy.ct.pred,
                       my.iter = my.iter, 
                       n.iter = n.iter, 
                       I = I, type = 'ct')
elpd.ct.vals[, 3, 7] <- apply(log(ct.elpd.vals.1), 1, sum)

sm.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                       samples.pred = samples.pred, 
                       access.pred = access.sm.pred, 
                       tree.pred = tree.sm.pred, 
                       canopy.pred = canopy.sm.pred,
                       my.iter = my.iter, 
                       n.iter = n.iter, 
                       I = I, type = 'sm')
elpd.sm.vals[, 2, 7] <- apply(log(sm.elpd.vals.1), 1, sum)

#### TR + SM 
samples.fit <- samples.sm[[curr.set]]
samples.pred <- samples.tr.sm[[pred.set]]

access.fit <- c(covars.ct$Access[-ct.x.indices[[curr.set]]],
                covars.tr$Access[-tr.x.indices[[curr.set]]],
                covars.sm$Access[-sm.x.indices[[curr.set]]])
mean.access.fit <- mean(access.fit)
sd.access.fit <- sd(access.fit)
access.fit <- (access.fit - mean(access.fit)) / sd(access.fit)
access.pred <- c(covars.tr$Access[tr.x.indices[[curr.set]]],
                 covars.sm$Access[sm.x.indices[[curr.set]]])
access.pred <- (access.pred - mean.access.fit) / sd.access.fit 

tree.fit <- c(covars.ct$Tree[-ct.x.indices[[curr.set]]],
              covars.tr$Tree[-tr.x.indices[[curr.set]]],
              covars.sm$Tree[-sm.x.indices[[curr.set]]])
mean.tree.fit <- mean(tree.fit)
sd.tree.fit <- sd(tree.fit)
tree.fit <- (tree.fit - mean(tree.fit)) / sd(tree.fit)
tree.pred <- c(covars.tr$Tree[tr.x.indices[[curr.set]]],
               covars.sm$Tree[sm.x.indices[[curr.set]]])
tree.pred <- (tree.pred - mean.tree.fit) / sd.tree.fit 

canopy.fit <- c(covars.ct$Canopy[-ct.x.indices[[curr.set]]],
                covars.tr$Canopy[-tr.x.indices[[curr.set]]],
                covars.sm$Canopy[-sm.x.indices[[curr.set]]])
mean.canopy.fit <- mean(canopy.fit)
sd.canopy.fit <- sd(canopy.fit)
canopy.fit <- (canopy.fit - mean(canopy.fit)) / sd(canopy.fit)
canopy.pred <- c(covars.tr$Canopy[tr.x.indices[[curr.set]]],
                 covars.sm$Canopy[sm.x.indices[[curr.set]]])
canopy.pred <- (canopy.pred - mean.canopy.fit) / sd.canopy.fit 

n.iter <- 500
my.iter <- (nrow(samples.fit) - n.iter + 1):nrow(samples.fit)
J.tr <- length(tr.x.indices[[curr.set]])
J.sm <- length(sm.x.indices[[curr.set]])
access.tr.pred <- access.pred[1:J.tr]
access.sm.pred <- access.pred[1:J.sm + J.tr]
tree.tr.pred <- tree.pred[1:J.tr]
tree.sm.pred <- tree.pred[1:J.sm + J.tr]
canopy.tr.pred <- canopy.pred[1:J.tr]
canopy.sm.pred <- canopy.pred[1:J.sm + J.tr]

tr.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                       samples.pred = samples.pred, 
                       access.pred = access.tr.pred, 
                       tree.pred = tree.tr.pred, 
                       canopy.pred = canopy.tr.pred,
                       my.iter = my.iter, 
                       n.iter = n.iter, 
                       I = I, type = 'tr')
elpd.tr.vals[, 3, 7] <- apply(log(tr.elpd.vals.1), 1, sum)

sm.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                       samples.pred = samples.pred, 
                       access.pred = access.sm.pred, 
                       tree.pred = tree.sm.pred, 
                       canopy.pred = canopy.sm.pred,
                       my.iter = my.iter, 
                       n.iter = n.iter, 
                       I = I, type = 'sm')
elpd.sm.vals[, 3, 7] <- apply(log(sm.elpd.vals.1), 1, sum)

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

canopy.fit <- c(covars.ct$Canopy[-ct.x.indices[[curr.set]]],
                covars.tr$Canopy[-tr.x.indices[[curr.set]]],
                covars.sm$Canopy[-sm.x.indices[[curr.set]]])
mean.canopy.fit <- mean(canopy.fit)
sd.canopy.fit <- sd(canopy.fit)
canopy.fit <- (canopy.fit - mean(canopy.fit)) / sd(canopy.fit)
canopy.pred <- c(covars.ct$Canopy[ct.x.indices[[curr.set]]],
                 covars.tr$Canopy[tr.x.indices[[curr.set]]],
                 covars.sm$Canopy[sm.x.indices[[curr.set]]])
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
elpd.ct.vals[, 4, 7] <- apply(log(ct.elpd.vals.1), 1, sum)

tr.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                       samples.pred = samples.pred, 
                       access.pred = access.tr.pred, 
                       tree.pred = tree.tr.pred, 
                       canopy.pred = canopy.tr.pred,
                       my.iter = my.iter, 
                       n.iter = n.iter, 
                       I = I, type = 'tr')
elpd.tr.vals[, 4, 7] <- apply(log(tr.elpd.vals.1), 1, sum)

sm.elpd.vals.1 <- elpd(samples.fit = samples.fit, 
                       samples.pred = samples.pred, 
                       access.pred = access.sm.pred, 
                       tree.pred = tree.sm.pred, 
                       canopy.pred = canopy.sm.pred,
                       my.iter = my.iter, 
                       n.iter = n.iter, 
                       I = I, type = 'sm')
elpd.sm.vals[, 4, 7] <- apply(log(sm.elpd.vals.1), 1, sum)

## Save all!
save(elpd.ct.vals, elpd.tr.vals, 
     elpd.sm.vals, file = paste('cross-val-results-', 
                                 curr.set, '.R', sep = ''))
