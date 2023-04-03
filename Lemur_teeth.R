
#load your libraries
library(Rvcg)
library(Morpho)
library(rgl)
library(geomorph)
library(paleomorph)
library(geomorph)


###=========== LOADING DATA SET 1:THE LANDMARKED SKULL 
#Read in landmarks manually placed on the whole skull 

ntaxa<-4 ## number of specimens in the test data set - NB you also put this in the code below (x,y,z) 

#data set .pts from Checkpoint - read in 
ptslist<-dir(pattern='.pts',recursive=T)
ptsarray<-array(dim=c(70,3,4)) #dim=c(number of landmarks and semilandmarks, number of dimensions, number of specimens)
for(i in 1:length(ptslist))
{
  ptsarray[,,i]<-as.matrix(read.table(file=ptslist[i],skip=2,header=F,sep="",row.names=1))
}

#Pull the names from the .pts filesso they match the coords
#[3] stays the same 
#Make sure to navigate to the folder above with the .pts in it - in the code below, it is looking for a folder called .pts
dimnames(ptsarray)[3]<-list(
  substr(dir("./pts",pattern=".pts"),1,(nchar(dir("./pts",pattern=".pts"))-4)))


#Let's rename your ptsarry to arraylm (jus because this is what it's called in my code)
arraylm <- ptsarray

##### MISSING LANDMARKS #########
#Find a specimen with missing data before run the code and have a look at what it looks like with 'spheres3D'
#Then run the below code and have another look at the specimen (look at both the coordinates data and the spheres3D plot - hopefully the missing has gone
arraylm[which(arraylm==9999)] <- NA
arraylm <- estimate.missing(arraylm,method="TPS")

#Have a look at some specimens 
#Note that the radius is larger because we haven't Procrusted the data yet 
#Have a look especially at the specimens wehre missing data has been estimated - do they look ok?

spheres3d(arraylm[,,,3], radius = 0.1, col = 'pink') #this is looking at specimen number 3, have a play with colour and radius 

#Once you've worked with the missing data, you can then plot your PCA 

#First Procrustes the data (this removes all info on translation, rotation, and size - leaving shape which is what we want)
Y.gpa=gpagen(arraylm) #Procrustes
data=Y.gpa$coords #Subset out the coords 
size=Y.gpa$Csize #centroid size


#This runs a PCA and pulls out PC scores (you can also use this for the mesh later if you want to put nice shapes on the axes)
PCA.ALL <- gm.prcomp(data)
plot(PCA.ALL)

#to look at specimen position 
#These numbers correspond to the numbers of your specimen as they are listed in your .pts folder 
text(PCA.ALL$x[,1],PCA.ALL$x[,2], pos = 4)



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


#Let's plot a better morphospace
library(dplyr)
library(ggplot2)


#LETS READ IN SOME ECOLOGICAL DATA TO MAKE THIS MORE INTERESTING 
#Add your own actual data for the specimens - this is just a fake data set 

species_data <- read.csv('lemur_data.csv')

#let's also look at the data to pull out the PC values (to put on the axes)
summary(PCA.ALL)

#pca_res <- plotTangentSpace(shape[,,-which(dimnames(shape)[[3]]%in%reconstructed)])
pca_scores <- as_tibble(PCA.ALL$x)
pca_scores <- pca_scores %>% mutate(., full=rownames(PCA.ALL$x)) 

#comp2<-species.data %>% filter(., Tip_Label %in% pca_scores$Species) %>% arrange(., Tip_Label) %>% filter(., !Filename%in%c('Crocodylus_acutus_AMNH_7857_cranium',"Crocodilus_biporcatus_cranium_MNHN_A_5316"))
#pca_scores <- pca_scores %>% mutate(., CladeA = comp2$Clade2, CladeB = as.factor(comp2$Clade),Extinct=as.factor(comp2$`extinct?`),Diet=as.factor(comp2$Diet2),Habitat=as.factor(comp2$Habitat))

library(hablar)
pca_scores2 <- left_join(species_data, pca_scores, by = c("full"))
pca_scores2 <- pca_scores2 %>% convert(fct(sex)) #this will colour code by sex 


pca_res <- PCA.ALL #the Procrusted, aligned data 

#Set the palette
cbPalette <- c("#481567FF", "#20A387FF") #only need the number of colours that you have for categories you're plotting, so sex = 2

g <- ggplot(pca_scores2, aes(x=Comp1, y=Comp2, color=sex))+
  geom_point()
g <-ggplot(pca_scores2, aes(x=Comp1, y=Comp2, fill=sex))+ 
  
  geom_point(size=6,aes(colour=family, fill= sex))+
  theme_bw()+
  theme(aspect.ratio = 1) +
  scale_alpha(range=c(0.4, 1),guide=FALSE)+
  #ggtitle("Principal Components Analysis") +
  
  xlab(paste0("PC Axis 1 (", signif((pca_res$pc.summary$importance[2,1]*100),3), "78.4 % of Total Variance)")) + 
  ylab(paste0("PC Axis 2 (",signif((pca_res$pc.summary$importance[2,2]*100),3),"8.6 % of Total Variance)")) +
  
  labs(fill = "family") + 
  scale_color_manual(values= cbPalette)

g

