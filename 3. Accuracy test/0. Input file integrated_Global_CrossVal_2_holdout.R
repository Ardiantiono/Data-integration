### R code for data input for Cross validation test (create two datasets; one half = predict, one half = holdout)
### This script is to prepare the dataset (species detection matrix and covariate)
### Script from Ardiantiono et al. 2024. Improved occupancy and cost-effectiveness of species monitoring programs through data integration
### Jags adaptation from Nimble script by Jeffrey W. Doser (doserjef@msu.edu)
#==================================================================================================================

rm(list = ls())

# Set working directory (where the all data required to run the models are stored)
# Please amend the filepath to ensure compatibility with your folder structure
setwd("D:\\Data analysis\\Chapter 2_Integrated\\Ardiantiono et al. 2024. Data integration\\3. Accuracy test")
getwd()

library(coda)
library(tidyverse)
library(jagsUI)

# Curr.set will determine which half of the dataset that we will run
# 1 = predict dataset, 2 = holdout dataset/model fit
curr.set <- 2 #THIS TIME WE RUN THE SECOND HALF


#==================================================================================================================
#==================================================================================================================
# Load the data ------------------------

# Load the patrol data
smart.dat <- read.table("Species_patrol_example.csv", header=TRUE, sep=",", na.strings=TRUE,
                        colClasses = c("factor", "numeric", "numeric", "numeric", "numeric", 
                                       "factor")) #specify class perhaps cause NAs

names(smart.dat)
head(smart.dat)
str(smart.dat)

# Convert dataframe into a 3-dimensional array
# Site (row) * Temporal Replicate (column) * Species (3rd dimension)
# Dimensions renamed for consistency and ease of recognition
smart.m <- simplify2array(by(smart.dat[2:5], smart.dat$Species, as.matrix, row.names=smart.dat$Subgrid_ID))
occasions <- rep(1:4,1) #create vector        
dimnames(smart.m)[[2]] <- occasions #name dim 2
smart.m #see the result

#==================================================================================================================
#Load the transect data 
tr.dat <- read.table("Species_transect_example.csv", header=TRUE, sep=",", na.strings=TRUE,
                     colClasses = c("factor", 
                                    "numeric", "numeric", "numeric", 
                                    "numeric", "numeric", "numeric", 
                                    "numeric", "numeric", "numeric", 
                                    "numeric", "numeric", "numeric", 
                                    "factor")) #specify class perhaps cause NAs

names(tr.dat)
head(tr.dat)
str(tr.dat)


# Convert dataframe into a 3-dimensional array
# Site (row) * Temporal Replicate (column) * Species (3rd dimension)
# Dimensions renamed for consistency and ease of recognition
tr.dat.m <- simplify2array(by(tr.dat[2:13], tr.dat$Species, as.matrix, row.names=tr.dat$Subgrid_ID))
occasions <- rep(1:12,1) #create vector        
dimnames(tr.dat.m)[[2]] <- occasions #name dim 2 (*actually simplify the matrix?)
tr.dat.m #see the result

#==================================================================================================================
# Load in the camera data. 
ct.dat <- read.table("Species_camera_example.csv", header=TRUE, sep=",", na.strings=TRUE,
                     colClasses = c("factor", 
                                    "numeric", "numeric", "numeric", "numeric", "numeric", "numeric",
                                    "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", 
                                    "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", 
                                    "factor")) #specify class  cause NAs
names(ct.dat)
head(ct.dat)
str(ct.dat)

# Convert dataframe into a 3-dimensional array
# Site (row) * Temporal Replicate (column) * Species (3rd dimension)
# Dimensions renamed for consistency and ease of recognition
ct.dat.m <- simplify2array(by(ct.dat[2:19], ct.dat$Species, as.matrix, row.names=ct.dat$Subgrid_ID))
occasions <- rep(1:18,1) #create vector        
dimnames(ct.dat.m)[[2]] <- occasions #name dim 2
ct.dat.m #see the result

#==================================================================================================================
#==================================================================================================================
# Prepare detection history data ---------------------
# Convert all count data to detection/non-detection data
dsm <- ifelse(smart.m > 0, 1, smart.m) #Patrol
dtr <- ifelse(tr.dat.m > 0, 1, tr.dat.m) #Transect
y <- ifelse(ct.dat.m > 0, 1, ct.dat.m) #Camera

#==================================================================================================================
# Define number of species
I <- dim(ct.dat.m)[3] #use camera data, all datasets have same number of species

# Number of sites
J.ct <- dim(ct.dat.m)[1] #CT
J.tr <- dim(tr.dat.m)[1] #SWTS
J.sm <- dim(smart.m)[1] #SMART

# Repeat visits at CT
K.ct <- dim(ct.dat.m)[2]
K.tr <- dim(tr.dat.m)[2]
K.sm <- dim(smart.m)[2]


###########################################################################
# Create Cross Validation Data Sets ---------------------------------------
# We split our dataset into half each

# SKIP THIS STEPS IF YOU HAVE RUN IT ALREADY
# INSTEAD call "load("0. cross-val-indices.R.R")"
#set.seed(2021)
#chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))
#ct.nums <- sample(1:J.ct, replace = FALSE)
#n.x.val <- 2
#ct.x.indices <- chunk2(ct.nums, n.x.val)

#tr.nums <- sample(1:J.tr, replace = FALSE)
#tr.x.indices <- chunk2(tr.nums, n.x.val)

#sm.nums <- sample(1:J.sm, replace = FALSE)
#sm.x.indices <- chunk2(sm.nums, n.x.val)

# Save cross validation sets to maintain consistency across models
# USE THIS AS REFERENCE AND CALL THIS WHEN RUNNING THE MODEL
# save(ct.x.indices, tr.x.indices, sm.x.indices, file = '0. cross-val-indices.R')

#Call the indices that we just created
load("0. cross-val-indices.R")


y.fit <- y[-ct.x.indices[[curr.set]], ,  ] #this is for holdout data (later model fit to see if predict model is accurate)
y.pred <- y[ct.x.indices[[curr.set]], ,  ] #this is for predict model (run the model then use it to predict in holdout locations)
dtr.fit <- dtr[-tr.x.indices[[curr.set]], , ]
dtr.pred <- dtr[tr.x.indices[[curr.set]], , ]
dsm.fit <- dsm[-sm.x.indices[[curr.set]], , ]
dsm.pred <- dsm[sm.x.indices[[curr.set]], , ]

# Site ID in number 
site.ct.fit <- as.numeric(rownames(y.fit))
site.ct.pred <- as.numeric(rownames(y.pred))
site.tr.fit <- as.numeric(rownames(dtr.fit))
site.tr.pred <- as.numeric(rownames(dtr.pred))
site.sm.fit <- as.numeric(rownames(dsm.fit))
site.sm.pred <- as.numeric(rownames(dsm.pred))

#length of sites
J.ct.fit <- length(site.ct.fit)
J.ct.pred <- length(site.ct.pred)
J.tr.fit <- length(site.tr.fit)
J.tr.pred <- length(site.tr.pred)
J.sm.fit <- length(site.sm.fit)
J.sm.pred <- length(site.sm.pred)

#==================================================================================================================
#==================================================================================================================
# Load in  the covariate------------------------------
covars.sm <- read.table("Covariates_patrol_example.csv", header=TRUE,sep=",",na.strings=c("NA"))
covars.tr <- read.table("Covariates_transect_example.csv", header=TRUE,sep=",",na.strings=c("NA"))
covars.ct <- read.table("Covariates_camera_example.csv", header=TRUE,sep=",",na.strings=c("NA"))

# We use: Accessibility, Tree Height, Tree Height Variability as occupancy covariates
# We use : Habitat type (except camera), TRI, and survey efforts as detection covariates
# We use: Year as temporal random factor

#-----------------------------------------------------
# OCCUPANCY COVARIATES

# Accessibility
access.ct.fit <- c(covars.ct$Access[-ct.x.indices[[curr.set]]])
access.ct.copy <- access.ct.fit
access.ct.fit <- (access.ct.fit - mean(access.ct.fit)) / sd(access.ct.fit)
access.ct.pred <- c(covars.ct$Access[ct.x.indices[[curr.set]]])
access.ct.pred <- (access.ct.pred - mean(access.ct.copy)) / sd(access.ct.copy)

access.tr.fit <- c(covars.tr$Access[-tr.x.indices[[curr.set]]])
access.tr.copy <- access.tr.fit
access.tr.fit <- (access.tr.fit - mean(access.tr.fit)) / sd(access.tr.fit)
access.tr.pred <- c(covars.tr$Access[tr.x.indices[[curr.set]]])
access.tr.pred <- (access.tr.pred - mean(access.tr.copy)) / sd(access.tr.copy)

access.sm.fit <- c(covars.sm$Access[-sm.x.indices[[curr.set]]])
access.sm.copy <- access.sm.fit
access.sm.fit <- (access.sm.fit - mean(access.sm.fit)) / sd(access.sm.fit)
access.sm.pred <- c(covars.sm$Access[sm.x.indices[[curr.set]]])
access.sm.pred <- (access.sm.pred - mean(access.sm.copy)) / sd(access.sm.copy)


# Tree height
tree.ct.fit <- c(covars.ct$Tree[-ct.x.indices[[curr.set]]])
tree.ct.copy <- tree.ct.fit
tree.ct.fit <- (tree.ct.fit - mean(tree.ct.fit)) / sd(tree.ct.fit)
tree.ct.pred <- c(covars.ct$Tree[ct.x.indices[[curr.set]]])
tree.ct.pred <- (tree.ct.pred - mean(tree.ct.copy)) / sd(tree.ct.copy)

tree.tr.fit <- c(covars.tr$Tree[-tr.x.indices[[curr.set]]])
tree.tr.copy <- tree.tr.fit
tree.tr.fit <- (tree.tr.fit - mean(tree.tr.fit)) / sd(tree.tr.fit)
tree.tr.pred <- c(covars.tr$Tree[tr.x.indices[[curr.set]]])
tree.tr.pred <- (tree.tr.pred - mean(tree.tr.copy)) / sd(tree.tr.copy)

tree.sm.fit <- c(covars.sm$Tree[-sm.x.indices[[curr.set]]])
tree.sm.copy <- tree.sm.fit
tree.sm.fit <- (tree.sm.fit - mean(tree.sm.fit)) / sd(tree.sm.fit)
tree.sm.pred <- c(covars.sm$Tree[sm.x.indices[[curr.set]]])
tree.sm.pred <- (tree.sm.pred - mean(tree.sm.copy)) / sd(tree.sm.copy)

# Tree height variance
canopy.ct.fit <- c(covars.ct$Tree_var[-ct.x.indices[[curr.set]]])
canopy.ct.copy <- canopy.ct.fit
canopy.ct.fit <- (canopy.ct.fit - mean(canopy.ct.fit)) / sd(canopy.ct.fit)
canopy.ct.pred <- c(covars.ct$Tree_var[ct.x.indices[[curr.set]]])
canopy.ct.pred <- (canopy.ct.pred - mean(canopy.ct.copy)) / sd(canopy.ct.copy)

canopy.tr.fit <- c(covars.tr$Tree_var[-tr.x.indices[[curr.set]]])
canopy.tr.copy <- canopy.tr.fit
canopy.tr.fit <- (canopy.tr.fit - mean(canopy.tr.fit)) / sd(canopy.tr.fit)
canopy.tr.pred <- c(covars.tr$Tree_var[tr.x.indices[[curr.set]]])
canopy.tr.pred <- (canopy.tr.pred - mean(canopy.tr.copy)) / sd(canopy.tr.copy)

canopy.sm.fit <- c(covars.sm$Tree_var[-sm.x.indices[[curr.set]]])
canopy.sm.copy <- canopy.sm.fit
canopy.sm.fit <- (canopy.sm.fit - mean(canopy.sm.fit)) / sd(canopy.sm.fit)
canopy.sm.pred <- c(covars.sm$Tree_var[sm.x.indices[[curr.set]]])
canopy.sm.pred <- (canopy.sm.pred - mean(canopy.sm.copy)) / sd(canopy.sm.copy)


# Year (specify random factor): we use subset year (by most recent data)
year.ct.fit <- c(covars.ct$Year[-ct.x.indices[[curr.set]]]) 
unique(year.ct.fit) #make sure have both year 2 & 3
n.year.ct.fit <- as.numeric(length(unique(year.ct.fit)))
year.ct.pred <- c(covars.ct$Year[ct.x.indices[[curr.set]]]) 
unique(year.ct.pred) #make sure have both year 2 & 3
n.year.ct.pred <- as.numeric(length(unique(year.ct.pred)))
#year.ct <- c(covars.ct$Year_CT) #USE THIS if running CT single model only (year 2-3 converted into year 1-2)

year.tr.fit <- c(covars.tr$Year[-tr.x.indices[[curr.set]]]) 
unique(year.tr.fit) #make sure have both year 1 & 2
n.year.tr.fit <- as.numeric(length(unique(year.tr.fit)))
year.tr.pred <- c(covars.tr$Year[tr.x.indices[[curr.set]]]) 
unique(year.tr.pred) #make sure have both year 1 & 2
n.year.tr.pred <- as.numeric(length(unique(year.tr.pred)))

year.sm.fit <- c(covars.sm$Year[-sm.x.indices[[curr.set]]]) 
unique(year.sm.fit) #make sure have both year 1, 2 & 3
n.year.sm.fit <- as.numeric(length(unique(year.sm.fit)))
year.sm.pred <- c(covars.sm$Year[sm.x.indices[[curr.set]]]) 
unique(year.sm.pred) #make sure have both year 1, 2 & 3
n.year.sm.pred <- as.numeric(length(unique(year.sm.pred)))

n.year <- as.numeric(length(unique(covars.sm$Year))) #the only dataset with three year period! For integrated model

#-----------------------------------------------------
# DETECTION COVARIATES

# TRI---------------------------
# TRI CT
tri.ct.fit <- c(covars.ct$Tri[-ct.x.indices[[curr.set]]])
tri.ct.copy <- tri.ct.fit
tri.ct.fit <- (tri.ct.fit - mean(tri.ct.fit)) / sd(tri.ct.fit)
tri.ct.fit <- cbind.data.frame((covars.ct$Site_ID[-ct.x.indices[[curr.set]]]), tri.ct.fit) 
names(tri.ct.fit) <- c('Site', 'tri') #we will integrate this into dataset
tri.ct.fit <- tri.ct.fit %>%
  mutate(Site = as.numeric(site.ct.fit),
         tri = as.numeric(tri))

tri.ct.pred <- c(covars.ct$Tri[ct.x.indices[[curr.set]]])
tri.ct.pred <- (tri.ct.pred - mean(tri.ct.copy)) / sd(tri.ct.copy)
tri.ct.pred <- cbind.data.frame((covars.ct$Site_ID[ct.x.indices[[curr.set]]]), tri.ct.pred) 
names(tri.ct.pred) <- c('Site', 'tri') #we will integrate this into dataset
tri.ct.pred <- tri.ct.pred %>%
  mutate(Site = as.numeric(site.ct.pred),
         tri = as.numeric(tri))

tri.tr.fit <- c(covars.tr$Tri[-tr.x.indices[[curr.set]]])
tri.tr.copy <- tri.tr.fit
tri.tr.fit <- (tri.tr.fit - mean(tri.tr.fit)) / sd(tri.tr.fit)
tri.tr.fit <- cbind.data.frame((covars.tr$Site_ID[-tr.x.indices[[curr.set]]]), tri.tr.fit) 
names(tri.tr.fit) <- c('Site', 'tri') #we will integrate this into dataset
tri.tr.fit <- tri.tr.fit %>%
  mutate(Site = as.numeric(site.tr.fit),
         tri = as.numeric(tri))
tri.tr.pred <- c(covars.tr$Tri[tr.x.indices[[curr.set]]])
tri.tr.pred <- (tri.tr.pred - mean(tri.tr.copy)) / sd(tri.tr.copy)
tri.tr.pred <- cbind.data.frame((covars.tr$Site_ID[tr.x.indices[[curr.set]]]), tri.tr.pred) 
names(tri.tr.pred) <- c('Site', 'tri') #we will integrate this into dataset
tri.tr.pred <- tri.tr.pred %>%
  mutate(Site = as.numeric(site.tr.pred),
         tri = as.numeric(tri))

tri.sm.fit <- c(covars.sm$Tri[-sm.x.indices[[curr.set]]])
tri.sm.copy <- tri.sm.fit
tri.sm.fit <- (tri.sm.fit - mean(tri.sm.fit)) / sd(tri.sm.fit)
tri.sm.fit <- cbind.data.frame((covars.sm$Site_ID[-sm.x.indices[[curr.set]]]), tri.sm.fit) 
names(tri.sm.fit) <- c('Site', 'tri') #we will integrate this into dataset
tri.sm.fit <- tri.sm.fit %>%
  mutate(Site = as.numeric(site.sm.fit),
         tri = as.numeric(tri))
tri.sm.pred <- c(covars.sm$Tri[sm.x.indices[[curr.set]]])
tri.sm.pred <- (tri.sm.pred - mean(tri.sm.copy)) / sd(tri.sm.copy)
tri.sm.pred <- cbind.data.frame((covars.sm$Site_ID[sm.x.indices[[curr.set]]]), tri.sm.pred) 
names(tri.sm.pred) <- c('Site', 'tri') #we will integrate this into dataset
tri.sm.pred <- tri.sm.pred %>%
  mutate(Site = as.numeric(site.sm.pred),
         tri = as.numeric(tri))

# HABITAT (categorical factor)-------------------------------

# Habitat CT (skip because only one category)

# Habitat SWTS
hab.tr.fit <- c(covars.tr$Hab[-tr.x.indices[[curr.set]]])
unique(hab.tr.fit) #make sure 4 covariates
hab.tr.fit <- cbind.data.frame((covars.tr$Site_ID[-tr.x.indices[[curr.set]]]), hab.tr.fit)
names(hab.tr.fit) <- c('Site', 'hab')
hab.tr.fit <- hab.tr.fit %>%
  mutate(Site = as.numeric(site.tr.fit),
         hab = as.numeric(hab))

hab.tr.pred <- c(covars.tr$Hab[tr.x.indices[[curr.set]]])
unique(hab.tr.pred) #make sure 4 covariates
hab.tr.pred <- cbind.data.frame((covars.tr$Site_ID[tr.x.indices[[curr.set]]]), hab.tr.pred)
names(hab.tr.pred) <- c('Site', 'hab')
hab.tr.pred <- hab.tr.pred %>%
  mutate(Site = as.numeric(site.tr.pred),
         hab = as.numeric(hab))

n.hab.tr <- as.numeric(length(unique(covars.tr$Hab)))

# Habitat SMART
hab.sm.fit <- c(covars.sm$Hab[-sm.x.indices[[curr.set]]])
unique(hab.sm.fit) #make sure 4 covariates
hab.sm.fit <- cbind.data.frame((covars.sm$Site_ID[-sm.x.indices[[curr.set]]]), hab.sm.fit)
names(hab.sm.fit) <- c('Site', 'hab')
hab.sm.fit <- hab.sm.fit %>%
  mutate(Site = as.numeric(site.sm.fit),
         hab = as.numeric(hab))

hab.sm.pred <- c(covars.sm$Hab[sm.x.indices[[curr.set]]])
unique(hab.sm.pred) #make sure 4 covariates
hab.sm.pred <- cbind.data.frame((covars.sm$Site_ID[sm.x.indices[[curr.set]]]), hab.sm.pred)
names(hab.sm.pred) <- c('Site', 'hab')
hab.sm.pred <- hab.sm.pred %>%
  mutate(Site = as.numeric(site.sm.pred),
         hab = as.numeric(hab))

n.hab.sm <- as.numeric(length(unique(covars.sm$Hab)))

#---------------------
#Trap days CT
trap.ct.fit <- c(covars.ct$CTN_90[-ct.x.indices[[curr.set]]])
trap.ct.copy <- trap.ct.fit
trap.ct.fit <- (trap.ct.fit - mean(trap.ct.fit)) / sd(trap.ct.fit)
trap.ct.fit <- cbind.data.frame((covars.ct$Site_ID[-ct.x.indices[[curr.set]]]), trap.ct.fit)
names(trap.ct.fit) <- c('Site', 'trap')
trap.ct.fit <- trap.ct.fit %>%
  mutate(Site = as.numeric(site.ct.fit),
         trap = as.numeric(trap))

trap.ct.pred <- c(covars.ct$CTN_90[ct.x.indices[[curr.set]]])
trap.ct.pred <- (trap.ct.pred - mean(trap.ct.copy)) / sd(trap.ct.copy)
trap.ct.pred <- cbind.data.frame((covars.ct$Site_ID[ct.x.indices[[curr.set]]]), trap.ct.pred)
names(trap.ct.pred) <- c('Site', 'trap')
trap.ct.pred <- trap.ct.pred %>%
  mutate(Site = as.numeric(site.ct.pred),
         trap = as.numeric(trap))


#CT Type (Categorical factor as intercept)
type.ct.fit <- c(covars.ct$CamType[-ct.x.indices[[curr.set]]])
unique(type.ct.fit) #make sure 2 types
type.ct.fit <- cbind.data.frame((covars.ct$Site_ID[-ct.x.indices[[curr.set]]]), type.ct.fit)
names(type.ct.fit) <- c('Site', 'type')
type.ct.fit <- type.ct.fit %>%
  mutate(Site = as.numeric(site.ct.fit),
         type = as.numeric(type))

type.ct.pred <- c(covars.ct$CamType[ct.x.indices[[curr.set]]])
unique(type.ct.pred) #make sure 2 types
type.ct.pred <- cbind.data.frame((covars.ct$Site_ID[ct.x.indices[[curr.set]]]), type.ct.pred)
names(type.ct.pred) <- c('Site', 'type')
type.ct.pred <- type.ct.pred %>%
  mutate(Site = as.numeric(site.ct.pred),
         type = as.numeric(type))

n.type.ct <- as.numeric(length(unique(covars.ct$CamType)))

#-------------------------
# Transect length SWTS
leng.tr.fit <- c(covars.tr$Transect_length[-tr.x.indices[[curr.set]]])
leng.tr.copy <- leng.tr.fit
leng.tr.fit <- (leng.tr.fit - mean(leng.tr.fit)) / sd(leng.tr.fit)
leng.tr.fit <- cbind.data.frame((covars.tr$Site_ID[-tr.x.indices[[curr.set]]]), leng.tr.fit) 
names(leng.tr.fit) <- c('Site', 'length') #we will integrate this into dataset
leng.tr.fit <- leng.tr.fit %>%
  mutate(Site = as.numeric(site.tr.fit),
         length = as.numeric(length))

leng.tr.pred <- c(covars.tr$Transect_length[tr.x.indices[[curr.set]]])
leng.tr.pred <- (leng.tr.pred - mean(leng.tr.copy)) / sd(leng.tr.copy)
leng.tr.pred <- cbind.data.frame((covars.tr$Site_ID[tr.x.indices[[curr.set]]]), leng.tr.pred) 
names(leng.tr.pred) <- c('Site', 'length') #we will integrate this into dataset
leng.tr.pred <- leng.tr.pred %>%
  mutate(Site = as.numeric(site.tr.pred),
         length = as.numeric(length))

#---------------------------
# Number of waypoint SMART
wp.sm.fit <- c(covars.sm$Num_WP[-sm.x.indices[[curr.set]]])
wp.sm.copy <- wp.sm.fit
wp.sm.fit <- (wp.sm.fit - mean(wp.sm.fit)) / sd(wp.sm.fit)
wp.sm.fit <- cbind.data.frame((covars.sm$Site_ID[-sm.x.indices[[curr.set]]]), wp.sm.fit) 
names(wp.sm.fit) <- c('Site', 'wp') #we will integrate this into dataset
wp.sm.fit <- wp.sm.fit %>%
  mutate(Site = as.numeric(site.sm.fit),
         wp = as.numeric(wp))

wp.sm.pred <- c(covars.sm$Num_WP[sm.x.indices[[curr.set]]])
wp.sm.pred <- (wp.sm.pred - mean(wp.sm.copy)) / sd(wp.sm.copy)
wp.sm.pred <- cbind.data.frame((covars.sm$Site_ID[sm.x.indices[[curr.set]]]), wp.sm.pred) 
names(wp.sm.pred) <- c('Site', 'wp') #we will integrate this into dataset
wp.sm.pred <- wp.sm.pred %>%
  mutate(Site = as.numeric(site.sm.pred),
         wp = as.numeric(wp))

##############################################################################
# Get All Data in Long Format ---------------------------------------------
# Modeling the data in this long format drastically speeds up run times (*Oooooo)
# CT Leuser ----------------------------------
y.fit.df <- as.data.frame.table(y.fit) #1620 observation, half of dataset
names(y.fit.df) <- c('Site', 'Visit', 'Species', 'Count')
y.fit.df$SiteID <- as.numeric(y.fit.df$Site) #create ordered site ID, keep unique site numbers below
y.fit.df <- y.fit.df %>%
  mutate(Site = as.numeric(as.character(Site)),
         Visit = as.numeric(Visit),
         Species = as.numeric(Species)) #change all data into numeric
y.fit.df <- y.fit.df %>% #remove NA values
  filter(!is.na(Count))
n.vals.ct.fit <- nrow(y.fit.df) #1608

y.pred.df <- as.data.frame.table(y.pred) 
names(y.pred.df) <- c('Site', 'Visit', 'Species', 'Count')
y.pred.df$SiteID <- as.numeric(y.pred.df$Site)
y.pred.df <- y.pred.df %>%
  mutate(Site = as.numeric(as.character(Site)), #make as.character to keep the value the same
         Visit = as.numeric(Visit),
         Species = as.numeric(Species)) #change all data into numeric
y.pred.df <- y.pred.df %>% #remove NA values
  filter(!is.na(Count))
n.vals.ct.pred <- nrow(y.pred.df) 

dtr.fit.df <- as.data.frame.table(dtr.fit) #6228 half of dataset
names(dtr.fit.df) <- c('Site', 'Visit', 'Species', 'Count')
dtr.fit.df$SiteID <- as.numeric(dtr.fit.df$Site)
dtr.fit.df <- dtr.fit.df %>%
  mutate(Site = as.numeric(as.character(Site)),
         Visit = as.numeric(Visit),
         Species = as.numeric(Species)) #change all data into numeric
dtr.fit.df <- dtr.fit.df %>% #remove NA values
  filter(!is.na(Count))
n.vals.tr.fit <- nrow(dtr.fit.df) #3222

dtr.pred.df <- as.data.frame.table(dtr.pred)
names(dtr.pred.df) <- c('Site', 'Visit', 'Species', 'Count')
dtr.pred.df$SiteID <- as.numeric(dtr.pred.df$Site)
dtr.pred.df <- dtr.pred.df %>%
  mutate(Site = as.numeric(as.character(Site)), #make as.character to keep the value the same
         Visit = as.numeric(Visit),
         Species = as.numeric(Species)) #change all data into numeric
dtr.pred.df <- dtr.pred.df %>% #remove NA values
  filter(!is.na(Count))
n.vals.tr.pred <- nrow(dtr.pred.df)

dsm.fit.df <- as.data.frame.table(dsm.fit) #1872 half of dataset
names(dsm.fit.df) <- c('Site', 'Visit', 'Species', 'Count')
dsm.fit.df$SiteID <- as.numeric(dsm.fit.df$Site)
dsm.fit.df <- dsm.fit.df %>%
  mutate(Site = as.numeric(as.character(Site)),
         Visit = as.numeric(Visit),
         Species = as.numeric(Species)) #change all data into numeric
dsm.fit.df <- dsm.fit.df %>% #remove NA values
  filter(!is.na(Count))
n.vals.sm.fit <- nrow(dsm.fit.df) #879

dsm.pred.df <- as.data.frame.table(dsm.pred) 
names(dsm.pred.df) <- c('Site', 'Visit', 'Species', 'Count')
dsm.pred.df$SiteID <- as.numeric(dsm.pred.df$Site)
dsm.pred.df <- dsm.pred.df %>%
  mutate(Site = as.numeric(as.character(Site)), #make as.character to keep the value the same
         Visit = as.numeric(Visit),
         Species = as.numeric(Species)) #change all data into numeric
dsm.pred.df <- dsm.pred.df %>% #remove NA values
  filter(!is.na(Count))
n.vals.sm.pred <- nrow(dsm.pred.df)


#include species parameters---------------------------
n.sp.ct <- length(unique(y.fit.df$Species)) #num species #use fit as pred also have same species number
n.sp.tr <- length(unique(dtr.fit.df$Species)) #num species
n.sp.sm <- length(unique(dsm.fit.df$Species)) #num species

#-----------------------------------------------
# join detection covariates into full data
y.fit.df <- left_join(y.fit.df, tri.ct.fit, by = c('Site'))
y.fit.df <- left_join(y.fit.df, trap.ct.fit, by = c('Site'))
y.fit.df <- left_join(y.fit.df, type.ct.fit, by = c('Site'))

y.pred.df <- left_join(y.pred.df, tri.ct.pred, by = c('Site'))
y.pred.df <- left_join(y.pred.df, trap.ct.pred, by = c('Site'))
y.pred.df <- left_join(y.pred.df, type.ct.pred, by = c('Site'))

dtr.fit.df<- left_join(dtr.fit.df, tri.tr.fit, by = c('Site'))
dtr.fit.df<- left_join(dtr.fit.df, leng.tr.fit, by = c('Site'))
dtr.fit.df<- left_join(dtr.fit.df, hab.tr.fit, by = c('Site'))

dtr.pred.df<- left_join(dtr.pred.df, tri.tr.pred, by = c('Site'))
dtr.pred.df<- left_join(dtr.pred.df, leng.tr.pred, by = c('Site'))
dtr.pred.df<- left_join(dtr.pred.df, hab.tr.pred, by = c('Site'))

dsm.fit.df <- left_join(dsm.fit.df, tri.sm.fit, by = c('Site'))
dsm.fit.df <- left_join(dsm.fit.df, wp.sm.fit, by = c('Site'))
dsm.fit.df <- left_join(dsm.fit.df, hab.sm.fit, by = c('Site'))

dsm.pred.df <- left_join(dsm.pred.df, tri.sm.pred, by = c('Site'))
dsm.pred.df <- left_join(dsm.pred.df, wp.sm.pred, by = c('Site'))
dsm.pred.df <- left_join(dsm.pred.df, hab.sm.pred, by = c('Site'))

##############################################################################
# Get indices necessary for Bayesian p-value
# low.bp: lower index for each species and year (we only have one year)
# high.bp: higher index for each species and year
y.fit.df <- y.fit.df %>% arrange(Species)
low.bp.y.fit <- c(matrix(NA, nrow = I, ncol = 1)) #bcs we only have one year
high.bp.y.fit <- c(matrix(NA, nrow = I, ncol = 1))

y.pred.df <- y.pred.df %>% arrange(Species)
low.bp.y.pred <- c(matrix(NA, nrow = I, ncol = 1)) #bcs we only have one year
high.bp.y.pred <- c(matrix(NA, nrow = I, ncol = 1))

for (i in 1:I) {
  tmp <- which(y.fit.df$Species == i)
  low.bp.y.fit[i] <- tmp[1]
  high.bp.y.fit[i] <- tmp[length(tmp)]
} # i

for (i in 1:I) {
  tmp <- which(y.pred.df$Species == i)
  low.bp.y.pred[i] <- tmp[1]
  high.bp.y.pred[i] <- tmp[length(tmp)]
} # i


#--------------------------
dtr.fit.df <- dtr.fit.df %>% arrange(Species)
low.bp.dtr.fit <- c(matrix(NA, nrow = I, ncol = 1)) 
high.bp.dtr.fit <- c(matrix(NA, nrow = I, ncol = 1))

dtr.pred.df <- dtr.pred.df %>% arrange(Species)
low.bp.dtr.pred <- c(matrix(NA, nrow = I, ncol = 1)) 
high.bp.dtr.pred <- c(matrix(NA, nrow = I, ncol = 1))

for (i in 1:I) {
  tmp <- which(dtr.fit.df$Species == i)
  low.bp.dtr.fit[i] <- tmp[1]
  high.bp.dtr.fit[i] <- tmp[length(tmp)]
} # i

for (i in 1:I) {
  tmp <- which(dtr.pred.df$Species == i)
  low.bp.dtr.pred[i] <- tmp[1]
  high.bp.dtr.pred[i] <- tmp[length(tmp)]
} # i

#--------------------------
dsm.fit.df <- dsm.fit.df %>% arrange(Species)
low.bp.dsm.fit <- c(matrix(NA, nrow = I, ncol = 1)) 
high.bp.dsm.fit <- c(matrix(NA, nrow = I, ncol = 1))

dsm.pred.df <- dsm.pred.df %>% arrange(Species)
low.bp.dsm.pred <- c(matrix(NA, nrow = I, ncol = 1)) 
high.bp.dsm.pred <- c(matrix(NA, nrow = I, ncol = 1))

for (i in 1:I) {
  tmp <- which(dsm.fit.df$Species == i)
  low.bp.dsm.fit[i] <- tmp[1]
  high.bp.dsm.fit[i] <- tmp[length(tmp)]
} # i

for (i in 1:I) {
  tmp <- which(dsm.pred.df$Species == i)
  low.bp.dsm.pred[i] <- tmp[1]
  high.bp.dsm.pred[i] <- tmp[length(tmp)]
} # i

##############################################################################


