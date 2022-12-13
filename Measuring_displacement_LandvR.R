
#Quantifying landmark displacement 

####################################
#                                  #
#   PENNA and COOMBS LMs           #
#                                  #
####################################

if(!require(devtools)) install.packages("devtools")
library(devtools)
install_github("TGuillerme/landvR")
library(landvR)

############# LANDVR ###################
#Check the values
AP_LMs[1,,2] == EC_LMs[1,,2]
# FALSE FALSE FALSE
## This is due to rounding, in fact they have the same 9 digits - round them 
round(AP_LMs[1,,2], digits = 9) == round(EC_LMs[1,,2], digits = 9)
# TRUE TRUE TRUE

#I’ve updated landvR to version 0.3 where the coordinates.difference function now have a tolerance optional argument.
#You can use the following to get the 0 difference results:

differences_between_lms <- coordinates.difference(coordinates = AP_LMs[,,1],
                                                  reference = EC_LMs[,,1],
                                                  type = "spherical",
                                                  rounding = 9)

#Remove errornous missing landmarks (these should be zero because they are static)
differences_between_lms[[1]][1:66, 1:3] <- c(0.000000, 0.000000, 0.000000)


#Ellen's colour function 
colfunc <- colorRampPalette(c("red", "yellow", "white"))
colfunc(10)
plot(rep(1,10),col=colfunc(10),pch=19,cex=3)

get.col.spectrum <- landvR::procrustes.var.plot(AP_LMs[,,1], EC_LMs[,,1], col.val = differences_between_lms[[1]][,1], col = colfunc)


test=differences_between_lms[[1]][,1] #this is a test for specimen 1 to look at the differences between lms 
test

##### LOOKING AT AN AVERAGE SPECIMEN ######

library(landvR)
N=26 #number of total landmarks 
nfixed = 26 #number of fixed landmarks 
specs = 3 #number of specimens 
k = 3 #number of dimensions in the matrix
colfunc <- colorRampPalette(c("red", "yellow", "white")) #create colour function for visualising landmarks
colfunc(10) #choose number of increments for colour scale 

#make the array for analyses: 
all_combined = array(dim=c(N,3,specs)) #3 is the columns of data we need (radii, azimuth, polar)
#manual_data is the fully landmarked skull (reference data)
#mirrored_data is the half-landmarked skull that has been mirrored 
#calculate the coordinates.differences between these data sets (i.e. how much the landmarks move between the manually placed landmarks and the mirrored landmarks) 

i=1
for (i in 1:specs)
{
  all_differences <- coordinates.difference(coordinates = AP_LMs[,,i],
                                            reference = EC_LMs[,,i],
                                            type = "spherical",
                                            rounding = 9)
  
  all_combined[,,i]=all_differences[[1]]
  
  i=i+1
}


#landmarks 1:66 (nfixed) in this example are fixed and therefore have the value of zero 
all_combined[1:nfixed, 1:k, 1:specs] <- c(0.000000, 0.000000, 0.000000)

#save output if desired using: write.csv(all_combined, file = 'all_combined.csv')

radii=all_combined[,1,] #radius per landmark for each specimen (second column of whole dataset with just the radii [,1,])

radii_mean=apply(radii, c(1), mean) #c(1) look at the first column - the radii 
#radii_mean is a mean radius value per landmark 

#save radii and radii_mean as .csv files for further analyses 

#example
#looking at the average radii compared to specimen 21 (or choose an average specimen if preferred)
get.col.spectrum <- landvR::procrustes.var.plot(AP_LMs[,,1], EC_LMs[,,1], col.val = radii_mean, col = colfunc)

#Checking for significant differences bwtween the two sides 
#Code from Agnese Lanzetti - see https://github.com/AgneseLan for details 

#Matrix with just the radius column from landVR output
radii <- all_combined[,1,]

#Create array to be filled with radii only
radii_mirror_array <- array(dim = c(nrow(radii), 1, ncol(radii))) #1 is the columns of data we need  - see output of previous line

#Loop linear distances calculations
for (b in 1:ncol(radii)){
  radii_mirror_array[,,b] <- radii[,b]
  
}
#Set dimnames as specimens
dimnames(radii_mirror_array)[[3]] <- dimnames(AP_LMs)[[3]]


#Create array to be filled with radii only
radii_mirror_real <- array(dim = c(nrow(radii), 1, ncol(radii))) #1 is the columns of data we need  - see output of previous line

#Overall variance in the dataset - max specimen - min specimen
var_range_radii <- rowSums(apply(radii_mirror_real, c(1,2), var))

landmarks_fix <- AP_LMs[c(1:26),,]

#Test difference between sides in entire dataset - is the variance on the right higher than on the left+midline?
t_test_LMs_diff <- t.test(var_range_radii[landmarks_fix], var_range_radii[-landmarks_fix])

#Summary with p-value
t_test_LMs_diff

t.test(COORDS[,,c(1:3)], COORDS[,,c(4:6)])
