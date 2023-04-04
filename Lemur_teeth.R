#Reading in PTS files, running Procrustes alignment, a PCA, and plotting a morphospace 

#load your libraries
library(Rvcg)
library(Morpho)
library(rgl)
library(geomorph)
library(paleomorph)
library(geomorph)


###=========== LOADING DATA SET 1:THE LANDMARKED SKULL/DATA
#Read in landmarks manually placed on the whole skull 

ntaxa<-26 ## number of specimens in the test data set - NB you also put this in the code below (x,y,z) 

#data set .pts from Checkpoint - read in 
ptslist<-dir(pattern='.pts',recursive=T)
ptsarray<-array(dim=c(70,3,26)) #dim=c(number of landmarks and semilandmarks, number of dimensions, number of specimens)
for(i in 1:length(ptslist))
{
  ptsarray[,,i]<-as.matrix(read.table(file=ptslist[i],skip=2,header=F,sep="",row.names=1))
}

#Pull the names from the .pts files so they match the coords
#[3] stays the same 
#!! Make sure to navigate to the folder above with the .pts in it - in the code below, it is looking for a folder called .pts
dimnames(ptsarray)[3]<-list(
  substr(dir("./pts",pattern=".pts"),1,(nchar(dir("./pts",pattern=".pts"))-4)))


#Let's rename your ptsarry to arraylm (just because this is what it's called in my code and we'll keep consistent)
arraylm <- ptsarray

##### MISSING LANDMARKS #########
#Find a specimen with missing data before run the code and have a look at what it looks like with 'spheres3D'
#Then run the below code and have another look at the specimen (look at both the coordinates data and the spheres3D plot - hopefully the missing has gone
arraylm[which(arraylm==9999)] <- NA
arraylm <- estimate.missing(arraylm,method="TPS")

#Have a look at some specimens - CHECK EACH ONE for any mistakes
#Note that the radius is larger because we haven't Procrusted the data yet 
#Have a look especially at the specimens where missing data has been estimated - do they look ok?

spheres3d(arraylm[,,3], radius = 0.1, col = 'pink') #this is looking at specimen number 3, have a play with colour and radius 

#Once you've worked with the missing data, you can then plot your PCA 

#First Procrustes the data (this removes all info on translation, rotation, and size - leaving shape which is what we want)
Y.gpa=gpagen(arraylm) #Procrustes
data=Y.gpa$coords #Subset out the coords 
size=Y.gpa$Csize #centroid size - if we want to look at size later


#This runs a PCA and pulls out PC scores (you can also use this for the mesh later if you want to put nice shapes on the axes)
PCA.ALL <- gm.prcomp(data)
plot(PCA.ALL)

#to look at specimen position 
#These numbers correspond to the numbers of your specimen as they are listed in your .pts folder 
text(PCA.ALL$x[,1],PCA.ALL$x[,2], pos = 4)
#Any thing odd? Take a look at specific specimens (using spheres3d) to check landmarks


#What's the shape of a mean specimen?
meanspec <- findMeanSpec(Y.gpa$coords)
ref <- Y.gpa$coords[, , meanspec]
plot(ref) #This is what a 'mean' species looks like 

#Let's have a look at what extreme specimens look like (shape) on PC21 and PC2 (the two axes)

#############
#Having a look at the min and max of the PCAs
#Have a play with the radius and colours - NOTE, the radius may is much smaller npw because these data have been Procrusted
spheres3d(PCA.ALL$shapes$shapes.comp1$min, radius=.008,color = "grey") #PC1 min (left side of x axis)
spheres3d(PCA.ALL$shapes$shapes.comp1$max, radius=.01,color = "grey") #PC2 max (right side of x axis)
spheres3d(PCA.ALL$shapes$shapes.comp2$min, radius=.008,color = "grey") #PC2 min (bottom of y axis)
spheres3d(PCA.ALL$shapes$shapes.comp2$max, radius=.01,color = "grey") #PC2 max (top of y axis)

#Save the output or take a screen shot, I'd love to see these!

#Let's plot a better morphospace
library(dplyr)
library(ggplot2)

#First of all, let's pull out the PC scores from the data so we can put them on the axes
summary(PCA.ALL)
#You need to pull out 'proportion of variance' for each comp PCA 
#eg 0.897 = 89.7% 
#Put these numbers in the ggplot2 code below


#pca_res <- plotTangentSpace(shape[,,-which(dimnames(shape)[[3]]%in%reconstructed)])
pca_scores <- as_tibble(PCA.ALL$x)
pca_scores <- pca_scores %>% mutate(., full=rownames(PCA.ALL$x)) 


# Plot the data
ggplot(pca_scores, aes(x = Comp1, y = Comp2)) +
  geom_point(size = 3, color = "blue") +
  xlab("PC1 (65.1%)") + ylab("PC2 (15.7%)") +
  ggtitle("Morphospace Plot")

