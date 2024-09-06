### R code for for extracting landscape-wide occupancy and detection probability
### Combining three datasets: patrol, transect, and camera (using example datasets)
### Use this output to visualise occupancy and detection probability with 95% BCIs
### Script from Ardiantiono et al. 2024. Improved occupancy and cost-effectiveness of species monitoring programs through data integration
#==================================================================================================================

rm(list = ls())

# Set working directory (where the all data required to run the models are stored)
# Please amend the filepath to ensure compatibility with your folder structure
setwd("D:\\Data analysis\\Chapter 2_Integrated\\Ardiantiono et al. 2024. Data integration\\1. Integrated occupancy model\\")
getwd()


#Call the ICOM model output
jags.out <- load("D:\\Data analysis\\Chapter 2_Integrated\\Ardiantiono et al. 2024. Data integration\\1. Integrated occupancy model\\Output Icom Example.RData")
#will be read as "out" instead of jags.out

# Load in  the covariate (for calculating site proportion)
covars.sm <- read.table("Covariates_patrol_example.csv", header=TRUE,sep=",",na.strings=c("NA"))
covars.tr <- read.table("Covariates_transect_example.csv", header=TRUE,sep=",",na.strings=c("NA"))
covars.ct <- read.table("Covariates_camera_example.csv", header=TRUE,sep=",",na.strings=c("NA"))

#==================================================================================================================
#==================================================================================================================
# To predict occupancy across landscape 
# In grid cell level (for creating prediction maps)

# Load in  the covariate for all landscape grid cells
covars.pred <- read.table("All subgrid covariates.csv", header=TRUE,sep=",",na.strings=c("NA"))

nGrid <- as.numeric(length(unique(covars.pred$Subgrid_ID)))
#we have 500 grid cells that we will predict their occupancy

Site_ID <- c(covars.pred$Subgrid_ID)
predict.table <- cbind(Site_ID,
                       out$mean$Occu.pred[,1], out$sd$Occu.pred[,1], out$q2.5$Occu.pred[,1], out$q97.5$Occu.pred[,1], #SP1
                       out$mean$Occu.pred[,2], out$sd$Occu.pred[,2], out$q2.5$Occu.pred[,2], out$q97.5$Occu.pred[,2], #SP2
                       out$mean$Occu.pred[,3], out$sd$Occu.pred[,3], out$q2.5$Occu.pred[,3], out$q97.5$Occu.pred[,3]) #SP3

predict.df <- as.data.frame(predict.table)
names(predict.df)
names(predict.df) <- c('Site_ID',  
                       'SP1', 'SP1_SD', 'SP1_q2.5', 'SP1_q97.5', 
                       'SP2', 'SP2_SD', 'SP2_q2.5', 'SP2_q97.5', 
                       'SP3', 'SP3_SD', 'SP3_q2.5', 'SP3_q97.5')

# Save the outputs --------------
write.csv(predict.df, file=paste("Output occupancy grid cell estimates", ".csv", sep=""))


#==================================================================================================================
#==================================================================================================================
# To predict occupancy across landscape 

# Create template dataset
sc.list <- as.data.frame(c("SP1", "SP2", "SP3"))
I <- length(sc.list)
J <- nGrid #number of grid cells across landscape that we want to predict
Site_Data <- 1:nGrid #same length as J

for(j in 1:J) {
  SM.TR.CT_SP1 <- c(base::mean(out$sims.list$Occu.pred[,Site_Data[j],1]), sd(out$sims.list$Occu.pred[,Site_Data[j],1]), 
                    quantile(out$sims.list$Occu.pred[,Site_Data[j],1],                                  
                             probs = c(0.025, 0.125, 0.25, 0.5, 0.75, 0.875, 0.975, 1)))
  SM.TR.CT_SP2 <- c(base::mean(out$sims.list$Occu.pred[,Site_Data[j],2]), sd(out$sims.list$Occu.pred[,Site_Data[j],2]), 
                    quantile(out$sims.list$Occu.pred[,Site_Data[j],2],                                  
                             probs = c(0.025, 0.125, 0.25, 0.5, 0.75, 0.875, 0.975, 1)))
  SM.TR.CT_SP3 <- c(base::mean(out$sims.list$Occu.pred[,Site_Data[j],3]), sd(out$sims.list$Occu.pred[,Site_Data[j],3]), 
                    quantile(out$sims.list$Occu.pred[,Site_Data[j],3],                                  
                             probs = c(0.025, 0.125, 0.25, 0.5, 0.75, 0.875, 0.975, 1)))
}

###############
beta_dog <- rbind(SM.TR.CT_SP1, SM.TR.CT_SP2, SM.TR.CT_SP3) 
beta_dog <- cbind(sc.list, beta_dog)

#Change the column names
names(beta_dog) <- c('Species', 'Mean','SD','2.5%', '12.5%',
                       '25%', '50%', "75%", '87.5%', '97.5%', '100%')

# Save the outputs --------------
write.csv(beta_dog, file=paste("Output occupancy landscape estimates", ".csv", sep=""))



#==================================================================================================================
#==================================================================================================================
# To predict detection across landscape (using weighted average)

# Number of site
J.ct <- dim(covars.ct)[1] #CT
J.tr <- dim(covars.tr)[1] #SWTS
J.sm <- dim(covars.sm)[1] #SMART

# Calculate the proportion of each intercept category
n.type.ct <- 2 #camera data
n.type.hab <- 4 #patrol and transect data

ct.prop <- tabulate(covars.ct$CamType)/J.ct #(N per camera type/total grids)
habtr.prop <- tabulate(covars.tr$Hab)/J.tr #(N per habitat type/total grids)
habsm.prop <- tabulate(covars.sm$Hab)/J.sm #(N per habitat type/total grids)


#####################################################
# Create template for df 
beta_cat <- data.frame(matrix(ncol = 11, nrow = 0))

sc.list <- c("SP1", "SP2", "SP3")
I <- length(sc.list)


for(i in 1:I) {
  alpha.weight <- out$sims.list$alpha.0 #Camera
  gamma.1.weight <- out$sims.list$gamma.1.0 #Transect
  gamma.2.weight <- out$sims.list$gamma.2.0 #Patrol
  
  #Camera
  for(k in 1:n.type.ct) {
    alpha.weight[,,k] <- out$sims.list$alpha.0[,,k] * ct.prop[k]
  }
  alpha.weight <- apply(alpha.weight,c(1,2),base::sum) #sum in third dimension (p across CT type: weighted average)
  
  #Transect
  for(k in 1:n.type.hab) {
    gamma.1.weight[,,k] <- out$sims.list$gamma.1.0[,,k] * habtr.prop[k]
  }
  gamma.1.weight <- apply(gamma.1.weight,c(1,2),base::sum)
  
  #Patrol
  for(k in 1:n.type.hab) {
    gamma.2.weight[,,k] <- out$sims.list$gamma.2.0[,,k] * habsm.prop[k]
  }
  gamma.2.weight <- apply(gamma.2.weight,c(1,2),base::sum) #sum in third dimension (p across habitat: weighted average)
  
  #combine all method estimates
  det.est <- rbind(alpha.weight, gamma.1.weight, gamma.2.weight) 
  
  backtransformed.det <- plogis(det.est[,i]) #a backtransformed version of b0 (species intercept)
  p.ctw.bt <- c(base::mean(backtransformed.det), sd(backtransformed.det), quantile(backtransformed.det, 
                                                                                   probs = c(0.025, 0.125, 0.25, 0.5, 0.75, 0.875, 0.975, 1)))
  
  # Create new dataset
  sp  <- sc.list[i]
  beta.est <- as.data.frame(t(p.ctw.bt)) 
  beta.est.sp <- cbind(sp, beta.est) 
  beta_cat <- rbind(beta_cat, beta.est.sp)
}

# Change the column names
names(beta_cat) <- c('Species', 'Mean','SD','2.5%', '12.5%',
                     '25%', '50%', "75%", '87.5%', '97.5%', '100%')

# Save the output------------------------
write.csv(beta_cat, file=paste("Output detection landscape estimates", ".csv", sep=""))
