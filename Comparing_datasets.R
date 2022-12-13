
#Testing microscribed specimens against checkpoint specimens

library(Rvcg)
library(Morpho)
library(rgl)
library(geomorph)
library(paleomorph)
library(geomorph)

landmarks <- read.morphologika('AP_lemur_test.txt')

arraylm <- landmarks



#Check LM numbers match Anna's desciptions 
spheres3d(arraylm[,,2])
text3d(arraylm[c(1:52),,3], text = 1:52)


#Read in Ellen's LMs 

###=========== LOADING DATA SET 1: WHOLE LANDMARKED SKULL 
#Read in landmarks manually placed on the whole skull 

ntaxa<-3 ## number of specimens (extant only) - NB can also put this in the code below (x,y,z) 
#data set .pts from Checkpoint

ptslist<-dir(pattern='.pts',recursive=T)
ptsarray<-array(dim=c(52,3,3)) #dim=c(number of landmarks and semilandmarks, number of dimensions, number of specimens)
for(i in 1:length(ptslist))
{
  ptsarray[,,i]<-as.matrix(read.table(file=ptslist[i],skip=2,header=F,sep="",row.names=1))
}


#Pull the names from the .pts files 
#[3] stays the same 
dimnames(ptsarray)[3]<-list(
  substr(dir("./pts",pattern=".pts"),1,(nchar(dir("./pts",pattern=".pts"))-4)))


array2 <- ptsarray

#Have a look at the two data sets 
spheres3d(arraylm[,,2], col = 'green')
spheres3d(array2[,,2], col = 'red')


#Bind the two data sets
merged_dataset=abind::abind(arraylm, array2, along = 3) #column 3 is specimens - this binds the coords data in alphabetical order         


#Sort out missing data 
#Deal with missing LMs (estimated from surrounding specimens)
merged_dataset[which(merged_dataset==9999)] <- NA
merged_dataset <- estimate.missing(merged_dataset,method="TPS")


#Procrustes the data 
#Procrustes the data 

merged_PROC=geomorph::gpagen(merged_dataset) #Remove non-shape aspects 
COORDS=merged_PROC$coords #Subset out the coords 


#PCA
pca_res <- geomorph::gm.prcomp(COORDS)

#Basic plot with specimen numbers on - specimen numbers relate to the .csv file 
plot(pca_res, axis1 = 1, axis2 = 2) #to plot 
#to look at specimen position 
text(pca_res$x[,1],pca_res$x[,2], pos = 4)



spheres3d(COORDS[,,c(3)], radius = 0.004, col = 'red')
spheres3d(COORDS[,,c(6)], radius = 0.004, col = 'green')


text3d(COORDS[,,c(3)], text = 1:52)
text3d(COORDS[,,c(6)], text = 1:52)
texts3d(COORDS[,,c(8)], text = 1:52)


#HOT DOTS to check any error 
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




