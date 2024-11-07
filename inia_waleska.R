#morphometric code for Waleska

#load the packages 
library(geomorph)
library(paleomorph)


setwd("X:xxxxxxx/xxxxxx") 
#NB the plys have been transformed to ASCII (from binary) so that landmarks can be visualised on the mesh 
#See code 'Binary_ASCII_ply.R" if needed for this step

ntaxa<-36 ## number of specimens (extant only) - NB can also put this in the code below (x,y,z) 
#data set .pts from Checkpoint

ptslist<-dir(pattern='.pts',recursive=T)
ptsarray<-array(dim=c(26,3,36)) #dim=c(number of landmarks and semilandmarks, number of dimensions, number of specimens)
for(i in 1:length(ptslist))
{
  ptsarray[,,i]<-as.matrix(read.table(file=ptslist[i],skip=2,header=F,sep="",row.names=1))
}

#Set to the whole project directory for this - the code wants to look for 'ply' folder 
dimnames(ptsarray)[3]<-list(
  substr(dir("./pts",pattern=".pts"),1,(nchar(dir("./pts",pattern=".pts"))-4)))###donne nom de scan a ptsarray et -4 pour retirer derniere lettre de nom
arraylm <- ptsarray
arraylm #run to check how this looks


##### CHANGING MISSING LANDMARKS BEFORE MIRRORING #########
#Otherwise you get weird numbers that aren't 9999 mirroring
#estimate.missing from Geomorph is used for estimating missing landmarks 

arraylm[which(arraylm==9999)] <- NA
arraylm <- estimate.missing(arraylm,method="TPS")


#have a look at how these look without the mesh first 
checkLM(arraylm,path="./ply/", pt.size = 2,suffix=".ply",render="s", begin = 1)

#or can check manually 
spheres3d(arraylm[,,1],radius=5) #check lms


#Check lms on the mesh 


#remove translation, rotation, and size info
#Procrustes the data 
arraylm_proc <- arraylm
Y.gpa=gpagen(arraylm_proc) #Procrustes
data=Y.gpa$coords #Subset out the coords 
size=Y.gpa$Csize #centroid size

plot(data)


#PCA - run the pca
pca.all <- gm.prcomp(data)
plot(pca.all) #simple morphospace

#to look at specimen position 
text(pca.all$x[,1],pca.all$x[,2], pos = 4)


spheres3d(pca.all$shapes$shapes.comp1$min, radius=.008,color = "grey") #PC1 min (left side of x axis)
spheres3d(pca.all$shapes$shapes.comp1$max, radius=.01,color = "grey") #PC2 max (right side of x axis)
spheres3d(pca.all$shapes$shapes.comp2$min, radius=.008,color = "grey") #PC2 min (bottom of y axis)
spheres3d(pca.all$shapes$shapes.comp2$max, radius=.01,color = "grey") #PC2 max (top of y axis)


#can also look at specimens on the mesh 
atarfa=ply2mesh(file="D:/Inia skulls/ply binary/AMNH147503.ply")
shade3d(atarfa,col='white')

spheres3d(arraylm_proc[,,6], col = 'blue', radius = 1.5) #add specimen number 


#build a nicer morphospace w/ ecological data if you have

library(geomorph)
library(rgl)
library(ggplot2)
library(hablar)
library(tidyr)
library(dplyr)
library(ggrepel)

#Read in your ecological data spreadsheet (this should be an excel file saved as a .csv file). 
#Add a column in your data sheet called 'landmarking_side' and then put 'left' or 'right' next to each of the specimens depending on which side you landmarked 
#Change the name of the species column (with the species names in) to 'full', this is how we join the shape data with the species data. 
#Make sure the name matches your .pts names 

pca_res <- PCA.ALL #the Procrusted, aligned data

species_data <- read.csv('EC_cleaned_data.csv', fileEncoding="UTF-8-BOM") #this is your species excel csv with the new coloum you added 

pca_scores <- as_tibble(pca_res$x)   #this is your PCA data as before 
pca_scores <- pca_scores %>% mutate(., full=rownames(pca_res$x))

pca_scores2 <- left_join(species_data, pca_scores, by = c("full")) 
pca_scores2 <- pca_scores2 %>% convert(fct(age)) #pulling out landmarking_side as the factor (could also pull out diet, sex, age, etc)

#Set the palette
cbPalette <- c("#481567FF", "#20A387FF", "green", "pink", "black", "blue", "yellow", "brown", "orange",
               'red', 'grey') #You only need two colours, for 'left' and 'right' 

g <- ggplot(pca_scores2, aes(x=Comp1, y=Comp2, color=age))+
  geom_point()
g <-ggplot(pca_scores2, aes(x=Comp1, y=Comp2, fill=age))+
  geom_point(size=6,aes(colour=genetic_id, fill= genetic_id))+
  theme_bw()+
  theme(aspect.ratio = 1) +
  scale_alpha(range=c(0.4, 1),guide='none')+
  #ggtitle("Principal Components Analysis looking at landmarking side") +
  
  xlab(paste0("PC Axis 1 (", signif((pca_res$pc.summary$importance[2,1]*100),3), "% of Total Variance)")) +
  ylab(paste0("PC Axis 2 (",signif((pca_res$pc.summary$importance[2,2]*100),3),"% of Total Variance)")) +
  
  labs(fill = "genetic_id") +
  scale_color_manual(values= cbPalette)

g #plot

#To add specimen names if you like! Always helpeful to see what's where (this pulls the names from the 'full' column 
g2<- g + geom_text_repel(inherit.aes = FALSE, aes(fontface=1,x=Comp1, y=Comp2,label=sub("_", " ", label), size = 1 ))
g2

