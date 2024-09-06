### R code for data input for integrated occupancy model (using example datasets)
### The while script will be run in "2. Integrated occupancy_model script.R" file
### Script from Ardiantiono et al. 2024. Improved occupancy and cost-effectiveness of species monitoring programs through data integration
#==================================================================================================================

rm(list = ls())

# Set working directory (where the all data required to run the models are stored)
# Please amend the filepath to ensure compatibility with your folder structure
setwd("D:\\Data analysis\\Chapter 2_Integrated\\Ardiantiono et al. 2024. Data integration\\1. Integrated occupancy model\\")
getwd()

library(coda)
library(tidyverse)

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
access.c <- c(covars.ct$Access)
access.real.ct <- access.c
access.ct <- c(scale(access.c))

access.t <- c(covars.tr$Access)
access.real.tr <- access.t
access.tr <- c(scale(access.t))

access.s <- c(covars.sm$Access)
access.real.sm <- access.s
access.sm <- c(scale(access.s))


# Tree height
tree.c <- c(covars.ct$Tree)
tree.real.ct <- tree.c
tree.ct <- c(scale(tree.c))

tree.t <- c(covars.tr$Tree)
tree.real.tr <- tree.t
tree.tr <- c(scale(tree.t))

tree.s <- c(covars.sm$Tree)
tree.real.sm <- tree.s
tree.sm <- c(scale(tree.s))

# Tree height variance
canopy.c <- c(covars.ct$Tree_var)
canopy.real.ct <- canopy.c
canopy.ct <- c(scale(canopy.c))

canopy.t <- c(covars.tr$Tree_var)
canopy.real.tr <- canopy.t
canopy.tr <- c(scale(canopy.t))

canopy.s <- c(covars.sm$Tree_var)
canopy.real.sm <- canopy.s
canopy.sm <- c(scale(canopy.s))

#----------------------------
# Year (temporal random factor): we use subset year (by most recent data)
year.c <- c(covars.ct$Year) 
n.year.ct <- as.numeric(length(unique(covars.ct$Year))) 

year.tr <- c(covars.tr$Year) 
n.year.tr <- as.numeric(length(unique(covars.tr$Year)))

year.sm <- c(covars.sm$Year) 
n.year.sm <- as.numeric(length(unique(covars.sm$Year)))

n.year <- n.year.sm #the only dataset with three year period.

#------------------------
# DETECTION COVARIATES

# TRI---------------------------
tri.c <- c(covars.ct$Tri)
tri.real.ct <- tri.c
tri.ct <- c(scale(tri.c))
tri.ct <- cbind.data.frame(covars.ct$Site_ID, tri.ct) 
names(tri.ct) <- c('Site', 'tri') #we will integrate this into dataset
tri.ct <- tri.ct %>%
  mutate(Site = as.numeric(factor(Site)),
         tri = as.numeric(tri))

tri.t <- c(covars.tr$Tri)
tri.real.tr <- tri.t
tri.tr <- c(scale(tri.t))
tri.tr <- cbind.data.frame(covars.tr$Site_ID, tri.tr) 
names(tri.tr) <- c('Site', 'tri') #we will integrate this into dataset
tri.tr <- tri.tr %>%
  mutate(Site = as.numeric(factor(Site)),
         tri = as.numeric(tri))

tri.s <- c(covars.sm$Tri)
tri.real.sm <- tri.s
tri.sm <- c(scale(tri.s))
tri.sm <- cbind.data.frame(covars.sm$Site_ID, tri.sm) 
names(tri.sm) <- c('Site', 'tri') #we will integrate this into dataset
tri.sm <- tri.sm %>%
  mutate(Site = as.numeric(factor(Site)),
         tri = as.numeric(tri))


# Habitat (categorical factor)-------------------------------
# Skip camera dataset because only one habitat type (forest)

hab.t <- c(covars.tr$Hab)
hab.real.tr <- hab.t
hab.tr <- cbind.data.frame(covars.tr$Site_ID, hab.t)
names(hab.tr) <- c('Site', 'hab')
hab.tr <- hab.tr %>%
  mutate(Site = as.numeric(factor(Site)),
         hab = as.numeric(hab))

n.hab.tr <- as.numeric(length(unique(covars.tr$Hab)))

# Habitat SMART
hab.s <- c(covars.sm$Hab)
hab.real.sm <- hab.s
hab.sm <- cbind.data.frame(covars.sm$Site_ID, hab.s)
names(hab.sm) <- c('Site', 'hab')
hab.sm <- hab.sm %>%
  mutate(Site = as.numeric(factor(Site)),
         hab = as.numeric(hab))

n.hab.sm <- as.numeric(length(unique(covars.sm$Hab)))

# Survey efforts ---------------------

# Trap days camera
trap.c <- c(covars.ct$CTN_90)
trap.real.ct <- trap.c
trap.ct <- c(scale(trap.c))
trap.ct <- cbind.data.frame(covars.ct$Site_ID, trap.ct)
names(trap.ct) <- c('Site', 'trap')
trap.ct <- trap.ct %>%
  mutate(Site = as.numeric(factor(Site)),
         trap = as.numeric(trap))

# Camera brand type (Categorical factor as intercept)
type.c <- c(covars.ct$CamType)
type.real.ct <- type.c
type.ct <- cbind.data.frame(covars.ct$Site_ID, type.c)
names(type.ct) <- c('Site', 'type')
type.ct <- type.ct %>%
  mutate(Site = as.numeric(factor(Site)),
         type = as.numeric(type))

n.type.ct <- as.numeric(length(unique(covars.ct$CamType)))

#-------------------------
# Transect length 
leng.t <- c(covars.tr$Transect_length)
leng.real.tr <- leng.t
leng.tr <- c(scale(leng.t))
leng.tr <- cbind.data.frame(covars.tr$Site_ID, leng.tr) 
names(leng.tr) <- c('Site', 'length') #we will integrate this into dataset
leng.tr <- leng.tr %>%
  mutate(Site = as.numeric(factor(Site)),
         length = as.numeric(length))

#---------------------------
# Number of waypoint SMART
wp.s <- c(covars.sm$Num_WP)
wp.real.sm <- wp.s
wp.sm <- c(scale(wp.s))
wp.sm <- cbind.data.frame(covars.sm$Site_ID, wp.sm) 
names(wp.sm) <- c('Site', 'wp') #we will integrate this into dataset
wp.sm <- wp.sm %>%
  mutate(Site = as.numeric(factor(Site)),
         wp = as.numeric(wp))



#==================================================================================================================
#==================================================================================================================
# Prepare detection history data ---------------------
# Convert all count data to detection/non-detection data
dsm <- ifelse(smart.m > 0, 1, smart.m) #Patrol
dtr <- ifelse(tr.dat.m > 0, 1, tr.dat.m) #Transect
y <- ifelse(ct.dat.m > 0, 1, ct.dat.m) #Camera

#==================================================================================================================
# Get all data in long format 
# Modeling the data in this long format drastically speeds up run times

# Patrol
dsm.df <- as.data.frame.table(dsm) #7896
names(dsm.df) <- c('Site', 'Visit', 'Species', 'Count')
dsm.df <- dsm.df %>%
  mutate(Site = as.numeric(Site),
         Visit = as.numeric(Visit),
         Species = as.numeric(Species)) #change all data into numeric

# Transect
dtr.df <- as.data.frame.table(dtr)
names(dtr.df) <- c('Site', 'Visit', 'Species', 'Count')
dtr.df <- dtr.df %>%
  mutate(Site = as.numeric(Site),
         Visit = as.numeric(Visit),
         Species = as.numeric(Species)) #change all data into numeric

# Camera
y.df <- as.data.frame.table(y)
names(y.df) <- c('Site', 'Visit', 'Species', 'Count')
y.df <- y.df %>%
  mutate(Site = as.numeric(Site),
         Visit = as.numeric(Visit),
         Species = as.numeric(Species)) #change all data into numeric

#==================================================================================================================
# Remove rows without data
dsm.df <- dsm.df %>%
  filter(!is.na(Count))
n.vals.sm <- nrow(dsm.df) 

dtr.df <- dtr.df %>%
  filter(!is.na(Count))
n.vals.tr <- nrow(dtr.df) 

y.df <- y.df %>%
  filter(!is.na(Count))
n.vals.ct <- nrow(y.df)

#include species parameters---------------------------
n.sp.sm <- length(unique(dsm.df$Species)) #num species
n.sp.tr <- length(unique(dtr.df$Species)) #num species
n.sp.ct <- length(unique(y.df$Species)) #num species

#==================================================================================================================
# join detection covariates ---------------------
dsm.df <- left_join(dsm.df, tri.sm, by = c('Site'))
dsm.df <- left_join(dsm.df, wp.sm, by = c('Site'))
dsm.df <- left_join(dsm.df, hab.sm, by = c('Site'))

dtr.df<- left_join(dtr.df, tri.tr, by = c('Site'))
dtr.df<- left_join(dtr.df, leng.tr, by = c('Site'))
dtr.df<- left_join(dtr.df, hab.tr, by = c('Site'))

y.df <- left_join(y.df, tri.ct, by = c('Site'))
y.df <- left_join(y.df, trap.ct, by = c('Site'))
y.df <- left_join(y.df, type.ct, by = c('Site'))


#==================================================================================================================
#==================================================================================================================
# Get indices for Bayesian p-value ---------------------------
# low.bp: lower index for each species and year 
# high.bp: higher index for each species and year

dsm.df <- dsm.df %>% arrange(Species)
low.bp.dsm <- c(matrix(NA, nrow = I, ncol = 1)) #bcs we only have one year
high.bp.dsm <- c(matrix(NA, nrow = I, ncol = 1))

for (i in 1:I) {
  tmp <- which(dsm.df$Species == i)
  low.bp.dsm[i] <- tmp[1]
  high.bp.dsm[i] <- tmp[length(tmp)]
} # i

#--------------------------
dtr.df <- dtr.df %>% arrange(Species, Site, Visit) #make sure visit 1:n in each sites
low.bp.dtr <- c(matrix(NA, nrow = I, ncol = 1)) #bcs we only have one year
high.bp.dtr <- c(matrix(NA, nrow = I, ncol = 1))

for (i in 1:I) {
  tmp <- which(dtr.df$Species == i)
  low.bp.dtr[i] <- tmp[1]
  high.bp.dtr[i] <- tmp[length(tmp)]
} # i

#--------------------------
y.df <- y.df %>% arrange(Species)
low.bp.y <- c(matrix(NA, nrow = I, ncol = 1)) #because we only have one year
high.bp.y <- c(matrix(NA, nrow = I, ncol = 1))

for (i in 1:I) {
  tmp <- which(y.df$Species == i)
  low.bp.y[i] <- tmp[1]
  high.bp.y[i] <- tmp[length(tmp)]
} # i

#==================================================================================================================
#==================================================================================================================
# Prepare layer for prediction -----------------------

# Load in  the covariate.
covars.pred <- read.table("All subgrid covariates.csv", header=TRUE,sep=",",na.strings=c("NA"))

nGrid <- as.numeric(length(unique(covars.pred$Subgrid_ID)))
#we have 500 grid cells that we will predict their occupancy

Site_ID <- c(covars.pred$Subgrid_ID)


# Standardize the occupancy covariates

# Accessibility
access.leuser <- c(scale(covars.pred$Access_pred))

# Tree height
tree.leuser <- c(scale(covars.pred$Tree_pred))

# Tree height variance
canopy.leuser <- c(scale(covars.pred$TreeVar_pred))

#==================================================================================================================
#==================================================================================================================
######### end of the code