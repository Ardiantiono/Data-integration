### R code for data input for Cross validation test (create two datasets; one half = predict, one half = holdout)
### This script is to summarize cross-validation results
### Script from Ardiantiono et al. 2024. Improved occupancy and cost-effectiveness of species monitoring programs through data integration
### Jags adaptation from Nimble script by Jeffrey W. Doser (doserjef@msu.edu)
#==================================================================================================================

rm(list = ls())
library(coda)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(ggpubr)

# Set working directory (where the all data required to run the models are stored)
# Please amend the filepath to ensure compatibility with your folder structure
setwd("D:\\Data analysis\\Chapter 2_Integrated\\Ardiantiono et al. 2024. Data integration\\3. Accuracy test")
getwd()

# Read in results ---------------------------------------------------------
load("cross-val-results-1.R")
elpd.ct.vals.1 <- elpd.ct.vals
elpd.tr.vals.1 <- elpd.tr.vals
elpd.sm.vals.1 <- elpd.sm.vals
load("cross-val-results-2.R")
elpd.ct.vals.2 <- elpd.ct.vals
elpd.tr.vals.2 <- elpd.tr.vals
elpd.sm.vals.2 <- elpd.sm.vals

# Take average of two-fold cross-validation. 
elpd.ct.vals.icom <- (elpd.ct.vals.1 + elpd.ct.vals.2) / 2
elpd.tr.vals.icom <- (elpd.tr.vals.1 + elpd.tr.vals.2) / 2
elpd.sm.vals.icom <- (elpd.sm.vals.1 + elpd.sm.vals.2) / 2

# Species by model results
ct.sp.model <- apply(elpd.ct.vals.icom, c(1, 3), mean, na.rm = TRUE)
tr.sp.model <- apply(elpd.tr.vals.icom, c(1, 3), mean, na.rm = TRUE)
sm.sp.model <- apply(elpd.sm.vals.icom, c(1, 3), mean, na.rm = TRUE)

ct.sp.model_t <- t(ct.sp.model)
tr.sp.model_t <- t(tr.sp.model)
sm.sp.model_t <- t(sm.sp.model)

# Order of maximum ELPD values by species
# Use function rank instead
ct.orders <- apply(ct.sp.model, 1, function(x) rank(-x))
tr.orders <- apply(tr.sp.model, 1, function(x) rank(-x))
sm.orders <- apply(sm.sp.model, 1, function(x) rank(-x))

# Compute average rank for a given model and data set
# CT
apply(ct.orders, 1, mean, na.rm = TRUE)

# TR
apply(tr.orders, 1, mean, na.rm = TRUE)

# SM
apply(sm.orders, 1, mean, na.rm = TRUE)

# All
all.sp.model <- ct.sp.model + tr.sp.model + sm.sp.model
all.sp.model_t <- t(all.sp.model)
all.orders <- apply(all.sp.model, 1, function(x) rank(-x))
apply(all.orders, 1, mean, na.rm = TRUE)


# All species model results
ct.model <- apply(ct.sp.model, 2, sum, na.rm = TRUE)
tr.model <- apply(tr.sp.model, 2, sum, na.rm = TRUE)
sm.model <- apply(sm.sp.model, 2, sum, na.rm = TRUE)

# All species all data sets
ct.model + tr.model + sm.model
all.model <- ct.model + tr.model + sm.model

#merge dataset
cross.val.df <- as.data.frame(cbind(ct.sp.model_t, ct.orders,
                tr.sp.model_t, tr.orders,
                sm.sp.model_t, sm.orders,
                all.sp.model_t, all.orders))
write.csv(cross.val.df, file = "1. cross.validation result.csv", row.names = FALSE)
