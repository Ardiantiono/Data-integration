### R code for data input for Cross validation test (create two datasets; one half = predict, one half = holdout)
### This script is to run occupancy model with half data
### Script from Ardiantiono et al. 2024. Improved occupancy and cost-effectiveness of species monitoring programs through data integration
### Jags adaptation from Nimble script by Jeffrey W. Doser (doserjef@msu.edu)
#==================================================================================================================


rm(list = ls())

# Set working directory (where the all data required to run the models are stored)
# Please amend the filepath to ensure compatibility with your folder structure
setwd("D:\\Data analysis\\Chapter 2_Integrated\\Ardiantiono et al. 2024. Data integration\\3. Accuracy test")
getwd()

# First, let's read our data preparation file
# Please amend the filepath to ensure compatibility with your folder structure
# RUN SEPARATELY as we want to run each half of the dataset

# 1 holdout = current set 1 
source("D:\\Data analysis\\Chapter 2_Integrated\\Ardiantiono et al. 2024. Data integration\\3. Accuracy test\\0. Input file integrated_Global_CrossVal_1_holdout.R") 

# 2 holdout = current set 2 
source("D:\\Data analysis\\Chapter 2_Integrated\\Ardiantiono et al. 2024. Data integration\\3. Accuracy test\\0. Input file integrated_Global_CrossVal_2_holdout.R") 

##---------------------------------------------------------------------------
#Add the model script

sink("Icom_CT_SWTS_SMART.txt")
cat("

model{

  # Priors ------------------------------------------------------------------
  
  # Occurrence ------------------------
  int.beta.mean ~ dunif(0, 1) # overall  
  beta.0.mean <- logit(int.beta.mean)
  beta.1.mean ~ dnorm(0, 0.368) # linear access
  beta.2.mean ~ dnorm(0, 0.368) # quadratic access
  beta.3.mean ~ dnorm(0, 0.368) # linear tree height
  beta.4.mean ~ dnorm(0, 0.368) # quadratic tree height
  beta.5.mean ~ dnorm(0, 0.368) # linear tree variance
  
  #-----------------------------------
  # Detection -------------------------
  
  # Overall CT detection--------------------
  # Intercept = Slope Camera Type (effect parameterization)
  # First we specify one hyperparameter for each camera type
    for(h in 1:n.type.ct){
       alpha.0.mean[h] ~ dnorm(0, 0.368)
       tau.alpha.0[h] ~ dgamma(0.1, 0.1) #add tau in the loop
    }
       
  alpha.1.mean ~ dnorm(0, 0.368) # TRI
  alpha.2.mean ~ dnorm(0, 0.368) # Trap days
       
       
  # Overall SWTS detection--------------------
  # Intercept = Slope Habitat Type (effect parameterization)
  # First we specify one hyperparameter for each habitat type
    for(h in 1:n.hab.tr){
       gamma.1.0.mean[h] ~ dnorm(0, 0.368)
       tau.gamma.1.0[h] ~ dgamma(0.1, 0.1) #add tau in the loop
    }
      
  gamma.1.1.mean ~ dnorm(0, 0.368) # TRI
  gamma.1.2.mean ~ dnorm(0, 0.368) # Transect length


  # overall SMART detection--------------------
  # Intercept = Slope Habitat Type (effect parameterization)
  # First we specify one hyperparameter for each habitat type
    for(h in 1:n.hab.sm){
       gamma.2.0.mean[h] ~ dnorm(0, 0.368)
       tau.gamma.2.0[h] ~ dgamma(0.1, 0.1) #add tau in the loop
    }

  gamma.2.1.mean ~ dnorm(0, 0.368) # TRI
  gamma.2.2.mean ~ dnorm(0, 0.368) # WP
  
  
  # Hyperprior for half-Cauchy scale parameter for process model ---------
  xi.tau <- pow(xi.sd, -2)
  xi.sd ~ dunif(0, 10) 
     
  #-----------------------------------
  # Precision Parameters --------------
  tau.beta.0 ~ dgamma(0.1, 0.1)
  tau.beta.1 ~ dgamma(0.1, 0.1)
  tau.beta.2 ~ dgamma(0.1, 0.1)
  tau.beta.3 ~ dgamma(0.1, 0.1)
  tau.beta.4 ~ dgamma(0.1, 0.1)
  tau.beta.5 ~ dgamma(0.1, 0.1)
  
  
  tau.alpha.1 ~ dgamma(0.1, 0.1)
  tau.alpha.2 ~ dgamma(0.1, 0.1)

  tau.gamma.1.1 ~ dgamma(0.1, 0.1)
  tau.gamma.1.2 ~ dgamma(0.1, 0.1)

  tau.gamma.2.1 ~ dgamma(0.1, 0.1)
  tau.gamma.2.2 ~ dgamma(0.1, 0.1)
  
  
  #-----------------------------------
  # Species-specific coefficients -----------------------------------------
  for (i in 1:I) {
    beta.0[i] ~ dnorm(beta.0.mean, tau.beta.0)
    beta.1[i] ~ dnorm(beta.1.mean, tau.beta.1)
    beta.2[i] ~ dnorm(beta.2.mean, tau.beta.2)
    beta.3[i] ~ dnorm(beta.3.mean, tau.beta.3)
    beta.4[i] ~ dnorm(beta.4.mean, tau.beta.4)
    beta.5[i] ~ dnorm(beta.5.mean, tau.beta.5)
    
    
    # Year-related Random effect
    for(year in 1:n.year){
        eps[i,year] ~ dnorm(0, eps.tau[i])
        }

      eps.tau[i] ~ dgamma(0.5, 0.5)
      xi[i] ~ dnorm(0, xi.tau)
      sigma.cauchy[i] <- abs(xi[i]) / sqrt(eps.tau[i])
  
  } #i
  
  #CT detection
  for (i in 1:I) {
    alpha.1[i] ~ dnorm(alpha.1.mean, tau.alpha.1)
    alpha.2[i] ~ dnorm(alpha.2.mean, tau.alpha.2)
    
    # For the sums-to-zero approach we need to specify a temporary parameter drawn from the hyperparameters
    # Then we subtract the mean of all categories from this
    for(h in 1:n.type.ct){
        alpha.0[i,h] ~ dnorm(alpha.0.mean[h], tau.alpha.0[h])
      }
  } # i
  
  #SWTS detection
  for (i in 1:I) {
    gamma.1.1[i] ~ dnorm(gamma.1.1.mean, tau.gamma.1.1)
    gamma.1.2[i] ~ dnorm(gamma.1.2.mean, tau.gamma.1.2)   
    
    # For the sums-to-zero approach we need to specify a temporary parameter drawn from the hyperparameters
    # Then we subtract the mean of all categories from this
    for(h in 1:n.hab.tr){
        gamma.1.0[i,h] ~ dnorm(gamma.1.0.mean[h], tau.gamma.1.0[h])
      }
  } # i
  
  #SMART detection
  for (i in 1:I) {
    gamma.2.1[i] ~ dnorm(gamma.2.1.mean, tau.gamma.2.1)
    gamma.2.2[i] ~ dnorm(gamma.2.2.mean, tau.gamma.2.2)
    
    # For the sums-to-zero approach we need to specify a temporary parameter drawn from the hyperparameters
    # Then we subtract the mean of all categories from this
    for(h in 1:n.hab.sm){
        gamma.2.0[i,h] ~ dnorm(gamma.2.0.mean[h], tau.gamma.2.0[h])
      }
  } #i
  
  
  #-----------------------------------
  # Process Models and Likelihoods ----------------------------------------
  # CT
  for (j in 1:J.ct) {
    for (i in 1:I) {
      logit(psi.ct[i, j]) <- beta.0[i] + 
        beta.1[i] * ACCESS.ct[j] + beta.2[i] * pow(ACCESS.ct[j], 2) +
        beta.3[i] * TREE.ct[j] + beta.4[i] * pow(TREE.ct[j], 2) +
        beta.5[i] * CANOPY.ct[j] + 
        xi[i]*eps[i,YEAR.ct[j]] #we added year random effect
        
      z.ct[i, j] ~ dbern(psi.ct[i, j])
    } # i
  } # j


  # SWTS
  for (j in 1:J.tr) {
    for (i in 1:I) {
       logit(psi.tr[i, j]) <- beta.0[i] + 
                             beta.1[i] * ACCESS.tr[j] + beta.2[i] * pow(ACCESS.tr[j], 2) +
                             beta.3[i] * TREE.tr[j] + beta.4[i] * pow(TREE.tr[j], 2) +
                             beta.5[i] * CANOPY.tr[j] + 
                             xi[i]*eps[i,YEAR.tr[j]] #we added year random effect
      
      z.tr[i, j] ~ dbern(psi.tr[i, j])
    } # i
  } # j


  # SMART
  for (j in 1:J.sm) {
    for (i in 1:I) {
      logit(psi.sm[i, j]) <- beta.0[i] + 
        beta.1[i] * ACCESS.sm[j] + beta.2[i] * pow(ACCESS.sm[j], 2) +
        beta.3[i] * TREE.sm[j] + beta.4[i] * pow(TREE.sm[j], 2) +
        beta.5[i] * CANOPY.sm[j] + 
        xi[i]*eps[i,YEAR.sm[j]] #we added year random effect
      
      z.sm[i, j] ~ dbern(psi.sm[i, j])
    } # i
  } # j
  
  
#-------------------------------------------------------------------------------  
  # Camera Trap Leuser Data ----------------------------------------------------
  # Data are stacked in a single vector as opposed to a multi-dimensional 
  # array to improve computational performance. 
  for (i in 1:n.vals.ct) {
    logit(p.ct[i]) <- alpha.0[sp.indx.ct[i], type.ct[i]] +
      alpha.1[sp.indx.ct[i]] * TRI.ct[i] +
      alpha.2[sp.indx.ct[i]] * trap.ct[i] 
    
    y[i] ~ dbern(p.ct[i] * z.ct[sp.indx.ct[i], site.ct[i]]) #*y=det/non
    
    # Replicate of data for Bayesian p-value
    y.rep[i] ~ dbern(p.ct[i] * z.ct[sp.indx.ct[i], site.ct[i]]) #expected
    E.ct[i] <- p.ct[i] * z.ct[sp.indx.ct[i], site.ct[i]]
    TRI.ct[i] ~ dnorm(0, 1)
  } # i
  
  
  # Transect Leuser Data ----------------------------------------------------
  # Data are stacked in a single vector as opposed to a multi-dimensional 
  # array to improve computational performance. 
  # Create two loops (1st replicate and subsequent replicates)
  
  for(k1 in 1:n.vals.tr) {                                                                                 ## Nick comment: looping over all sites where segment == 1

      TRI.tr[k1] ~ dnorm(0, 1)                                                                       ## Nick comment: What is TRI.tr? It looks like a covariate drawn from a normal distribution - what does this account for?
      logit(p.tr[k1]) <- gamma.1.0[sp.indx.tr[k1], hab.tr[k1]] + 
                         gamma.1.1[sp.indx.tr[k1]] * TRI.tr[k1] +
                         gamma.1.2[sp.indx.tr[k1]] * leng.tr[k1] 
           
          dtr[k1] ~ dbern(p.tr[k1] * z.tr[sp.indx.tr[k1], site.tr[k1]])
          
          # Replicate of data for Bayesian p-value
          dtr.rep[k1] ~ dbern(p.tr[k1] * z.tr[sp.indx.tr[k1], site.tr[k1]])
          E.tr[k1] <- p.tr[k1] * z.tr[sp.indx.tr[k1], site.tr[k1]]
  } #k1

     
  # SMART Leuser Data ----------------------------------------------------
  # Data are stacked in a single vector as opposed to a multi-dimensional 
  # array to improve computational performance. 
  for (i in 1:n.vals.sm) {
    logit(p.sm[i]) <- gamma.2.0[sp.indx.sm[i], hab.sm[i]] +
      gamma.2.1[sp.indx.sm[i]] * TRI.sm[i] +
      gamma.2.2[sp.indx.sm[i]] * wp.sm[i] 
    
    dsm[i] ~ dbern(p.sm[i] * z.sm[sp.indx.sm[i], site.sm[i]]) #*y=det/non
    
    # Replicate of data for Bayesian p-value
    dsm.rep[i] ~ dbern(p.sm[i] * z.sm[sp.indx.sm[i], site.sm[i]]) #expesmed
    E.sm[i] <- p.sm[i] * z.sm[sp.indx.sm[i], site.sm[i]]
    TRI.sm[i] ~ dnorm(0, 1)
  } # i
  
  
#----------------------------------------------------------------------------  
    # Compute Bayesian P-value GoF Assessment -------------------------------
    for (i in 1:I) {
      y.grouped[i] <- sum(y[low.bp.y[i]:high.bp.y[i]])
      y.rep.grouped[i] <- sum(y.rep[low.bp.y[i]:high.bp.y[i]])
      E.grouped.y[i] <- sum(E.ct[low.bp.y[i]:high.bp.y[i]])
      x2.y[i] <- pow((y.grouped[i] - E.grouped.y[i]), 2) / (E.grouped.y[i] + e) #equal to chi2.actual
      x2.y.rep[i] <- pow((y.rep.grouped[i] - E.grouped.y[i]), 2) / (E.grouped.y[i] + e) #equal to chi2.sim
    
      dtr.grouped[i] <- sum(dtr[low.bp.dtr[i]:high.bp.dtr[i]])
      dtr.rep.grouped[i] <- sum(dtr.rep[low.bp.dtr[i]:high.bp.dtr[i]])
      E.grouped.dtr[i] <- sum(E.tr[low.bp.dtr[i]:high.bp.dtr[i]])
      x2.dtr[i] <- pow((dtr.grouped[i] - E.grouped.dtr[i]), 2) / (E.grouped.dtr[i] + e) #equal to chi2.actual
      x2.dtr.rep[i] <- pow((dtr.rep.grouped[i] - E.grouped.dtr[i]), 2) / (E.grouped.dtr[i] + e) #equal to chi2.sim

      dsm.grouped[i] <- sum(dsm[low.bp.dsm[i]:high.bp.dsm[i]])
      dsm.rep.grouped[i] <- sum(dsm.rep[low.bp.dsm[i]:high.bp.dsm[i]])
      E.grouped.dsm[i] <- sum(E.sm[low.bp.dsm[i]:high.bp.dsm[i]])
      x2.dsm[i] <- pow((dsm.grouped[i] - E.grouped.dsm[i]), 2) / (E.grouped.dsm[i] + e) #equal to chi2.actual
      x2.dsm.rep[i] <- pow((dsm.rep.grouped[i] - E.grouped.dsm[i]), 2) / (E.grouped.dsm[i] + e) #equal to chi2.sim
    } #i 
  
    # Add up Chi-square discrepancies
    fit.actual.ct <- sum(x2.y[1:I]) #change the name from Nimble Doser et ak (initally chi.2.y)
    fit.sim.ct <- sum(x2.y.rep[1:I])

    fit.actual.tr <- sum(x2.dtr[1:I]) #change the name from Nimble Doser et ak (initally chi.2.y)
    fit.sim.tr <- sum(x2.dtr.rep[1:I])
    
    fit.actual.sm <- sum(x2.dsm[1:I]) #change the name from Nimble Doser et ak (initally chi.2.y)
    fit.sim.sm <- sum(x2.dsm.rep[1:I])
    
    # Calculate overall chi-squared discrepency measure
    #==================================================
    c.hat.ct <- fit.actual.ct/fit.sim.ct
    c.hat.tr <- fit.actual.tr/fit.sim.tr
    c.hat.sm <- fit.actual.sm/fit.sim.sm
    
    
    bpv.ct <- step(fit.sim.ct - fit.actual.ct) 
    bpv.tr <- step(fit.sim.tr - fit.actual.tr)
    bpv.sm <- step(fit.sim.sm - fit.actual.sm)
      
  
    # Derived quantities
    # Number of occupied sites per SPECIES (perhaps per data source)
    #=========================
    for(i in 1:I) {
        Nocc.fs.ct[i] <- sum(z.ct[i,]) 
        Nocc.fs.tr[i] <- sum(z.tr[i,])
        Nocc.fs.sm[i] <- sum(z.sm[i,])
    }
    
    
    # Average occupancy for each species
    #=========================================
    for(i in 1:I) {
        mean.psi.ct[i] <- mean(psi.ct[i,])
        mean.psi.tr[i] <- mean(psi.tr[i,])
        mean.psi.sm[i] <- mean(psi.sm[i,])
    }  
        
}

",fill=TRUE)
sink()


# 1d. Model preparation and initialization
#==================================================================================================================
# Load library to call JAGS from R
library(jagsUI) #need to use JAGs 3.0.1 for R > 4.2.0
library(plyr)
library(reshape2)
library(boot)  # for logit functions
library(verification)  # to calc AUC statistic (used in functions)


# Specify the fundamentals of the MCMC chains
ni <- 200000; nb <- 100000; nt <- 100; nc <- 3 # Final iterations
ni <- 1000; nb <- 0; nt <- 1; nc <- 3 # Exploratory iterations (nt=thinning rate, nc= chain)


# For reproduceability
set.seed(2000)


# Run univariate models across all covariates 
{
  m.fit <- setNames(data.frame(matrix(ncol = 8, nrow = 1)), c("Model", "DIC", "BPV", "Chat", "lppd", "pD", "WAIC", "CPO"))
  ### Analyze the data
  ###=================
  # Specify the parameters to be monitored  #what's the difference between int and beta.0 (logit of int)
  # Focus only psi and z
  params <- c('int.beta.mean', 'beta.0.mean', 'beta.1.mean', 
              'beta.2.mean', 'beta.3.mean', 'beta.4.mean',
              'beta.5.mean', 
              'beta.0', 'beta.1', 'beta.2', 'beta.3', 
              'beta.4', 'beta.5', 
              'z.ct', 'z.tr', 'z.sm',
              'bpv.ct', 'bpv.tr', 'bpv.sm')
  
  # Load the data
  data <- list(I = I, 
               J.ct = J.ct.fit, J.tr = J.tr.fit, J.sm = J.sm.fit,
               n.vals.ct = n.vals.ct.fit, n.vals.tr = n.vals.tr.fit, n.vals.sm = n.vals.sm.fit,
               sp.indx.ct = y.fit.df$Species, sp.indx.tr = dtr.fit.df$Species,  
               sp.indx.sm = dsm.fit.df$Species, 
               site.ct = y.fit.df$SiteID, site.tr = dtr.fit.df$SiteID, site.sm = dsm.fit.df$SiteID,
               low.bp.y = low.bp.y.fit, low.bp.dtr = low.bp.dtr.fit, low.bp.dsm = low.bp.dsm.fit,  
               high.bp.y = high.bp.y.fit, high.bp.dtr = high.bp.dtr.fit, high.bp.dsm = high.bp.dsm.fit, 
               e = 0.0001,
               ACCESS.ct = access.ct.fit, ACCESS.tr = access.tr.fit, ACCESS.sm = access.sm.fit,
               TREE.ct = tree.ct.fit, TREE.tr = tree.tr.fit, TREE.sm = tree.sm.fit, 
               CANOPY.ct = canopy.ct.fit, CANOPY.tr = canopy.tr.fit, CANOPY.sm = canopy.sm.fit,
               YEAR.ct = year.ct.fit, YEAR.tr = year.tr.fit, YEAR.sm = year.sm.fit, 
               n.year = n.year, 
               TRI.ct = y.fit.df$tri, TRI.tr = dtr.fit.df$tri, TRI.sm = dsm.fit.df$tri,
               hab.tr = dtr.fit.df$hab, n.hab.tr = n.hab.tr, 
               hab.sm = dsm.fit.df$hab, n.hab.sm = n.hab.sm,
               trap.ct = y.fit.df$trap, type.ct = y.fit.df$type, n.type.ct = n.type.ct,
               leng.tr = dtr.fit.df$length, 
               wp.sm = dsm.fit.df$wp,
               y = y.fit.df$Count, dtr = dtr.fit.df$Count, dsm = dsm.fit.df$Count)
  
  #Specify initial values
  z.ct.init <- array(1, dim=c(n.sp.ct, J.ct.fit))
  z.tr.init <- array(1, dim=c(n.sp.tr, J.tr.fit))
  z.sm.init <- array(1, dim=c(n.sp.sm, J.sm.fit))
  
  inits <- function() list(z.ct = z.ct.init, z.tr = z.tr.init, z.sm = z.sm.init,  
                           int.beta.mean = runif(1, 0.3, 0.9), 
                           beta.1.mean = rnorm(1), beta.2.mean = rnorm(1),
                           beta.3.mean = rnorm(1), beta.4.mean = rnorm(1),
                           beta.5.mean = rnorm(1), 
                           alpha.0.mean = rnorm(n.type.ct), alpha.1.mean = rnorm(1), 
                           alpha.2.mean = rnorm(1), 
                           gamma.1.0.mean = rnorm(n.hab.tr), gamma.1.1.mean = rnorm(1), 
                           gamma.1.2.mean = rnorm(1), 
                           gamma.2.0.mean = rnorm(n.hab.sm), gamma.2.1.mean = rnorm(1),
                           gamma.2.2.mean = rnorm(1), 
                           tau.beta.0 = runif(1, 0.1, 2), tau.beta.1 = runif(1, 0.1, 2),
                           tau.beta.2 = runif(1, 0.1, 2), tau.beta.3 = runif(1, 0.1, 2),
                           tau.beta.4 = runif(1, 0.1, 2), tau.beta.5 = runif(1, 0.1, 2), 
                           tau.alpha.0 = runif(n.type.ct, 0.1, 2), tau.alpha.1 = runif(1, 0.1, 2),
                           tau.alpha.2 = runif(1, 0.1, 2), 
                           tau.gamma.1.0 = runif(n.hab.tr, 0.1, 2), tau.gamma.1.1 = runif(1, 0.1, 2),
                           tau.gamma.1.2 = runif(1, 0.1, 2),
                           tau.gamma.2.0 = runif(n.hab.sm, 0.1, 2), tau.gamma.2.1 = runif(1, 0.1, 2), 
                           tau.gamma.2.2 = runif(1, 0.1, 2), 
                           xi.sd=runif(1,0,10), eps.tau=runif(n.sp.tr))
  
  #Run the model and extract the results
  print(paste('All data integrated Global', 'Time: ', Sys.time()))
  
  out <- jags(data, inits, params, "Icom_CT_SWTS_SMART.txt", n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nt, parallel=TRUE)
  save(out, file=paste('icom-CT-TR-SM-results-', ni, '-iterations-', curr.set, '-holdout-', 
                       '.R', sep = ''))
  jags.summ <- out$summary
  write.csv(jags.summ, file=paste('icom-CT-TR-SM-results-', ni, '-iterations-', curr.set, '-holdout-', 
                                  '.csv', sep=''))
  
}








