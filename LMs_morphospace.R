
#Read in LM data (pts files), check names match, Procrustes, PCA, morphospace (with ecological data if required)

#read in the libraries
library(geomorph)
library(rgl)
library(ggplot2)
library(hablar)
library(dplyr)

#Read in Ellen's LMs - these should be in a .pts folder

###=========== LOADING DATA SET 1:LE LANDMARKED SKULL 
#Read in landmarks manually placed on the whole skull 

ntaxa<-201 ## number of specimens in the test data set - NB you also put this in the code below (x,y,z) 

#data set .pts from Checkpoint - read in 
ptslist<-dir(pattern='.pts',recursive=T)
ptsarray<-array(dim=c(124,3,201)) #dim=c(number of landmarks, number of dimensions, number of specimens)
for(i in 1:length(ptslist))
{
  ptsarray[,,i]<-as.matrix(read.table(file=ptslist[i],skip=2,header=F,sep="",row.names=1))
}

#Pull the names from the .pts files so they match the coords data
#[3] stays the same 
dimnames(ptsarray)[3]<-list(
  substr(dir("./pts",pattern=".pts"),1,(nchar(dir("./pts",pattern=".pts"))-4)))


arraylm <- ptsarray #this is your data set

#Have a look at the data in 3D 
spheres3d(arraylm[,,2], col = 'black', radius = 2)

#Check it on the mesh (your mesh must be in ASCII .ply format)
dimnames(arraylm)[3] #check your specimen names

#Load ply file
mesh=ply2mesh(file="D:/Ply ASCII/ply ASCII/Mysts/ply/Diorocetus hiatus USNM 16783.ply")
shade3d(mesh,col='white')
spheres3d(arraylm[c(1:66),,21], radius = 2, color = 'black') #here I am just plotting LMs 1:66 on specimen 21 (Diorocetus hiatus USNM 16783)

#Pulling in data to build a nicer morphospace later
#Read in species data (this is your ecology data, geological age etc)
#MAKE SURE YOUR SHAPE DATA HAS THE SAME NAME AS YOUR DATA IN THE SPECIES DATA SPEADSHEET - otherwise the shape data and the eco data won't be combined
#I tend to give everything the same name as it has in the phylogeny (in 'species_data' there is a column called 'phylo_name). 

species_data <- read.csv('species_data.csv',fileEncoding="UTF-8-BOM")

#If needed: 

Full_names_phylo=species_data$phylo_name #pull the names out 
dimnames(arraylm)[[3]]<-Full_names_phylo
View(cbind(dimnames(arraylm)[[3]], species_data$phylo_name)) #check the names match 

#Note, this has not phylogenetically corrected the data, it has just given it the name of the specimen as it is in the phylogeny
#Double check your data are still in the same order using spheres3D, for example, is specimen #1 still plotting as a sperm whale (or whatever your specimen number 1 is)
#If you want to phylo correct before running the Procrustes, check the code here: https://github.com/EllenJCoombs/Matt_odontocetes/blob/master/Phylo_corrected_phylomorphospace.R

#Procrustes the data (not phylo correction, see above)
Y.gpa=gpagen(arraylm) #Procrustes
shape_data=Y.gpa$coords #Subset out the coords 
size=Y.gpa$Csize #Centroid size

plot(shape_data) #raw data

#This runs a PC and pulls out PC scores (you can also use this for the mesh)
PCA.ALL <- gm.prcomp(shape_data)
plot(PCA.ALL) #raw morphospace

#to look at specimen position (numbers)
text(PCA.ALL$x[,1],PCA.ALL$x[,2], pos = 4) #these numbers are the same as dimnames(arraylm)[3]


#Having a look at the min and max of the PCAs

meanspec <- findMeanSpec(Y.gpa$coords)
ref <- Y.gpa$coords[, , meanspec]
plot(ref) #This is what a 'mean' species looks like 

#Pull out the max and min shapes 
spheres3d(PCA.ALL$shapes$shapes.comp1$min, radius=.008,color = "grey")
spheres3d(PCA.ALL$shapes$shapes.comp1$max, radius=.01,color = "grey")
spheres3d(PCA.ALL$shapes$shapes.comp2$min, radius=.008,color = "grey")
spheres3d(PCA.ALL$shapes$shapes.comp2$max, radius=.01,color = "grey")


#Visualise the points for each specimen if wanted 
spheres3d(shape_data[,,4],radius= 0.0002) #note that the radius has changed because these data are Procrusted


#Building a nicer morphospace with ecological data 
pca_scores <- as_tibble(PCA.ALL$x)
pca_scores <- pca_scores %>% mutate(., full=rownames(PCA.ALL$x))

pca_scores2 <- left_join(species_data, pca_scores, by = c("full"))
pca_scores2 <- pca_scores2 %>% convert(fct(family)) #or whichever factor you are keen to look at 


pca_res <- PCA.ALL #the Procrusted, aligned data 

#Set the palette
cbPalette <- c("#481567FF", "#20A387FF", "#DCE319FF")

g <- ggplot(pca_scores2, aes(x=Comp1, y=Comp2, color=family))+
  geom_point()
g <-ggplot(pca_scores2, aes(x=Comp1, y=Comp2, fill=family))+ 
  
  geom_point(size=6,aes(colour=family, fill= family))+
  theme_bw()+
  theme(aspect.ratio = 1) +
  scale_alpha(range=c(0.4, 1),guide=FALSE)+
  #ggtitle("Principal Components Analysis") +
  
  xlab(paste0("PC Axis 1 (", signif((pca_res$pc.summary$importance[2,1]*100),3), "% of Total Variance)")) + 
  ylab(paste0("PC Axis 2 (",signif((pca_res$pc.summary$importance[2,2]*100),3),"% of Total Variance)")) +
  
  labs(fill = "family") + 
  scale_color_manual(values= cbPalette)

g #plot

#To add specimen names 
g2<- g + geom_text_repel(inherit.aes = FALSE, aes(fontface=1,x=Comp1, y=Comp2,label=sub("_", " ", full), size = 1 ))
g2



