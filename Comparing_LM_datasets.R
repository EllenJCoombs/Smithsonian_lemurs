
#Testing microscribed specimens against checkpoint specimens - how different is landmark placement? 
#One data set is read in in morphologika format (microscribed by Anna P), the second data set is in .pts Checkpoint format (landmarked by Ellen C). 
#Both data sets use the SAME landmarks added to the SAME specimens in the SAME order. 

#load your libraries
library(Rvcg)
library(Morpho)
library(rgl)
library(geomorph)
library(paleomorph)
library(geomorph)


#Read in AP's landmarks (landmark data set number 1)
landmarks <- read.morphologika('AP_lemur_test.txt') #this is how to read in morphologika data, see below for .pts from checkpoint

arraylm <- landmarks #first data set

#Check LM numbers match your desciptions and are placed where you would expect
spheres3d(arraylm[,,2]) #number denotes specimen number
text3d(arraylm[c(1:52),,3], text = 1:52) #This numbers the landmarks in the order they were placed

#Read in Ellen's LMs (landmark data set number 2)

###=========== LOADING DATA SET 1:LE LANDMARKED SKULL 
#Read in landmarks manually placed on the whole skull 

ntaxa<-3 ## number of specimens in the test data set - NB you also put this in the code below (x,y,z) 

#data set .pts from Checkpoint - read in 
ptslist<-dir(pattern='.pts',recursive=T)
ptsarray<-array(dim=c(52,3,3)) #dim=c(number of landmarks and semilandmarks, number of dimensions, number of specimens)
for(i in 1:length(ptslist))
{
  ptsarray[,,i]<-as.matrix(read.table(file=ptslist[i],skip=2,header=F,sep="",row.names=1))
}

#Pull the names from the .pts filesso they match the coords
#[3] stays the same 
dimnames(ptsarray)[3]<-list(
  substr(dir("./pts",pattern=".pts"),1,(nchar(dir("./pts",pattern=".pts"))-4)))


array2 <- ptsarray #this is your second data set

#Have a look at the two data sets 
spheres3d(arraylm[,,2], col = 'green')
spheres3d(array2[,,2], col = 'red')


#Bind the two data sets
#Top set of specimens are data set 1, bottom set of specimens are data set 2 (e.g. here, specimens 1-3 and 4-6, accordingly) 
merged_dataset=abind::abind(arraylm, array2, along = 3) #column 3 is specimens - this binds the coords data from the two data sets in alphabetical order         


#Sort out missing data 
#Deal with missing LMs (estimated from surrounding specimens)
merged_dataset[which(merged_dataset==9999)] <- NA
merged_dataset <- estimate.missing(merged_dataset,method="TPS")


#Procrustes the two data sets together (important!)

merged_PROC=geomorph::gpagen(merged_dataset) #Remove non-shape aspects 
COORDS=merged_PROC$coords #Subset out the coords 

#PCA - Principle components analysis
pca_res <- geomorph::gm.prcomp(COORDS)

#See how they look in the morphospace 
#Basic plot with specimen numbers on - specimen numbers relate to the .csv fil
#Here, we want 1 to sit near 4, 2 to sit near 5, and 3 to sit near 6 because of how we merged the two data sets (see above) - the closer the two data sets sit, the 
#more similar the landmarking
plot(pca_res, axis1 = 1, axis2 = 2) #to plot 
#to look at specimen position 
text(pca_res$x[,1],pca_res$x[,2], pos = 4)

#Look at one of the first specimens (data set 1) against one of the later specimens (data set 2) to compare 
spheres3d(COORDS[,,c(3)], radius = 0.004, col = 'red')
spheres3d(COORDS[,,c(6)], radius = 0.004, col = 'green')

#Do the same with numbering of landmarks to check everything is in the right place
text3d(COORDS[,,c(3)], text = 1:52)
text3d(COORDS[,,c(6)], text = 1:52)
texts3d(COORDS[,,c(8)], text = 1:52)


#Run a t-test between the two data sets - if the landmarks are placed correctly, this shouldn't be significant 
t.test(, postTreat, paired = TRUE)

#If you are getting a sig difference in the t.test result (the landmarks are not positioned as closely as you'd like) then run the following to look for erroneous landmarks 

#HOT DOTS to check any error - if you like
#Check distances between data sets 
plotRefToTarget(COORDS[,,1], COORDS[,,4], method = "vector", radius = 0.004)

# 1. Create colour palette (rainbow spectrum colours)
col.var <- colorRampPalette( c("violet","blue","cyan","yellow","red" ) )


# 2. Calculate variances for each individual landmark
variances <- rowSums(apply(COORDS, c(1,2), var))

# 3. Log variance values to give more continuous values
log.variances <- log( variances ) 

# 4. Create colour gradient of variance values based on colour palette created earlier
varianceCol <- col.var( 30 )[ as.numeric( cut( log.variances, breaks = 30 ) ) ]
# varianceCol <- col.gp2( 30 )[ as.numeric( cut( variances, breaks = 30 ) ) ]

# 5. Check output by plotting points as per-point variance
plot( log.variances, col = varianceCol, cex = 2, pch = 19 )

# 6. Create 3D plot of per-landmark variance, based on single reference specimen
plot3d( array2[,,3], col = varianceCol, type = "s",
        radius = 0.8, aspect = T, main = "mean",axes = F, main = F, fov = 0 )


