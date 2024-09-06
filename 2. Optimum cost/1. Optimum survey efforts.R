### R code for for calculating optimum number of grid cells
### With specified occupancy precision and detection range
### Adapted from Mackenzie & Royle 2005 (https://besjournals.onlinelibrary.wiley.com/doi/10.1111/j.1365-2664.2005.01098.x)
### Script from Ardiantiono et al. 2024. Improved occupancy and cost-effectiveness of species monitoring programs through data integration
#==================================================================================================================

rm(list = ls())

# Set working directory (where the all data required to run the models are stored)
setwd("D:\\Data analysis\\Chapter 2_Integrated\\Ardiantiono et al. 2024. Data integration\\2. Optimum cost\\")
getwd()

#####################
sc.list <- c("SP1", "SP2", "SP3")

# Occupancy precision estimates (we use standard error 0.05-0.50; 
# SE 0.05 = 10% deviation from occupancy estimates; SE 0.50 = 100% deviation)
prec <- seq(0.05, 0.50, by = 0.01)

# Occupancy (psi) estimates 
# Use mean estimates from model output "Output occupancy landscape estimates.csv"
PSI <- c(0.326, #SP1
         0.320, #SP2
         0.744) #SP3
         
        
# Detection estimates (separate between mean, minimum 2.5% BCI, and maximum 97.5% BCI)
# Use mean estimates from model output "Output detection landscape estimates.csv"
pmean <- c(0.150, #SP1
           0.099, #SP2
           0.078) #SP3

pmin <- c(0.023,
          0.003,
          0.044)

pmax <- c(0.412,
          0.287,
          0.125)

# Replicates (adjust based on dataset replicates)
# Patrol = 4, Transect = 12, Camera = 18 (sum replicates for integrated models)
SO <- c(34, 34, 34)

#################################
## Calculate optimum survey effort

# Create empty dataset
beta_cat <- data.frame(matrix(ncol = 139, nrow = 0)) #46 estimates x 3 species

for(k in 1:length(pmean)) {
  
  #Create empty vector
  vector1 <- c()
  vector2 <- c()
  vector3 <- c()
  
  for(j in 1:length(prec)) {
  ### MEAN *****************************
  occ.survey.optimiser <- function(psi, p, K, SE) {
    p1 <- 1- (1-pmean[k])^K #remember to adjust based on precision (SE)
    (psi/(SE^2))*((1-psi)+((1-p1)/(p1 - (K*(p*(1-p)^(K-1))))))
  }
  opt.mean <- occ.survey.optimiser(psi=PSI[k], p=pmean[k], K=SO[k], SE=prec[j])
  
  ### MIN *****************************
  occ.survey.optimiser <- function(psi, p, K, SE) {
    p1 <- 1- (1-pmin[k])^K #remember to adjust based on p-value
    (psi/(SE^2))*((1-psi)+((1-p1)/(p1 - (K*(p*(1-p)^(K-1))))))
  }
  opt.min <- occ.survey.optimiser(psi=PSI[k], p=pmin[k], K=SO[k], SE=prec[j]) 
  
  ### MAX *****************************
  occ.survey.optimiser <- function(psi, p, K, SE) {
    p1 <- 1- (1-pmax[k])^K #remember to adjust based on p-value
    (psi/(SE^2))*((1-psi)+((1-p1)/(p1 - (K*(p*(1-p)^(K-1))))))
  }
  opt.max <- occ.survey.optimiser(psi=PSI[k], p=pmax[k], K=SO[k], SE=prec[j])
     
  vector1 <- c(vector1, opt.mean)
  vector2 <- c(vector2, opt.max)
  vector3 <- c(vector3, opt.min)
  
  vector.all <- c(vector1, vector2, vector3)
  }
  
# Create new dataset
  beta.est <- as.data.frame(t(vector.all)) 
  beta_cat <- rbind(beta_cat, beta.est)
}


#### Now add species name and save!
beta.table <- cbind(sc.list, beta_cat)
# This table inform optimum number to survey based on specified occupancy estimates + precision
# And detection estimates and 95% BCI range
# columns 2-46: Mean optimum numbers (from SE 0.05 [high precision] - 0.50 [low precision]); 
# Columns 47-92: Minimum numbers of cells; 93-138: Maximum number of cells

###=================
write.csv(beta.table, file=paste("Output optimum site numbers", ".csv", sep=""))


#==================================================================================================================
# Calculate the optimum cost
# Here we use operational cost needed for one team to survey one grid cell
# Calculation refer to WebTable 2 in Supplementary Material

# For integrated model the cost = Sum total cost (patrol + transect + camera) / Sum replicates (patrol + transect + camera)
Cost = (3123	+ 1725 +	3090) / (12 + 8 + 30)

beta.cost <- beta.table[,2:138] * Cost
# This table inform the operational cost needed to survey optimum number of grid cells to achieve certain precision

###=================
write.csv(beta.cost, file=paste("Output optimum costs", ".csv", sep=""))


