#############################################
#                                           #
#  Code to import, resample, and slide LMs  #
#                                           #
#############################################

#load packages 

rm(list=ls()) #to clear workspace 
library(tidyverse)
library(broom) #needed with tidyverse
library(latticeExtra) #qgraph
library(readr)
library(Morpho)
library(geomorph)
library(Rvcg)
library(paleomorph)
library(EMMLi)
library(qgraph) #load in packages bit 
library(ape)
library(geiger)
library(abind)
library(tibble)
library(dplyr)
library(tidyr)
library(rgl)
library("devtools")
devtools::install_github("rnfelice/SURGE")
library(SURGE)
library(RColorBrewer) # for color palettes
library(magick)


#set path to PTS final LHS - this is where all of the curves are 

#import table defining curves
#different curves table for the odontocetes 
curve_table <- read_csv('new_curves.csv')

#identify the folder where your pts files are
ptsfolder <- "D:/Missy lemurs/Elizabeths_ PCA and lemur stuff/landmarks/pts"


#import the pts file names
ptslist <- dir(ptsfolder, pattern='.pts', recursive=F)

my_curves <- create_curve_info(curve_table, n_fixed = 70)
setwd(ptsfolder)
subsampled.lm <- import_chkpt_data(ptslist, my_curves, subsampl = TRUE, verbose=TRUE)

#capture the output if it's too big to read in the console 
capture.output(import_chkpt_data(ptslist, my_curves, subsampl = TRUE, verbose=TRUE), file = "odonts.txt", append = T)
#if you have any missing points, Checkpoint will set them to x,y,z=9999
#this makes for bad plots in checkLM, so switch them to NA
subsampled.lm[subsampled.lm == 9999] <- NA
subsampled.lm[subsampled.lm == Inf] <- NA
subsampled.lm[subsampled.lm == -Inf] <- NA

arraylm <- estimate.missing(subsampled.lm,method="TPS")

#SET WORKING DIRECTORY TO ASCII PLY
#check to make sure your curves look okay on each specimen
checkLM(subsampled.lm,path="./ply/", pt.size = 2,suffix=".ply",render="s", begin = 5)


newpts <- arraylm

#Create missing list 
misslist<-createMissingList(dim(newpts)[3])
for (j in 1:dim(newpts)[[3]]){
  misslist[[j]]<-which(is.na(newpts[,1,j]))
} 
newpts2<-fixLMtps(newpts)
#check that your pts and ply match
plyfolder <- "D:/Missy lemurs/Elizabeths_ PCA and lemur stuff"
ptslist2<-gsub(pattern="\\.pts$","",ptslist)
plylist <-  dir(plyfolder, pattern='.ply', recursive=F)
plylist2<-gsub(pattern="\\.ply$","",plylist)
setdiff(plylist2,ptslist2) #should be zero if all matching up OK
#######NOTE
#####my plys are awkwardly in a folder called  "C:/Data/whale_ply/ply/ply"
##### I am setting the working directory ONE FOLDER UP from where the meshes are

#mc.cores = 1 for WINDOWS. mc.cores = 3 for MAC
####so for you it should be 

#################################################
#                                               #
#        RUN THE SLIDER 3D_2 CODE               #
#                                               #
#################################################

##################################
#                                #
#    STEP SIZE - 0.2 or 0.002    #
#                                #
##################################


setwd("D:/Missy lemurs/Elizabeths_ PCA and lemur stuff")
slided4.all <- slider3d_2(newpts2, SMvector= my_curves$Sliding.LMs,
                          outlines = my_curves$Curves, sur.path = "./ply", sur.name = NULL, 
                          meshlist = paste("./ply/",dimnames(newpts2)[[3]],".ply",sep=""), ignore = NULL,
                          sur.type = "ply", tol = 1e-10, deselect = FALSE, inc.check = FALSE,
                          recursive = TRUE, iterations = 3, initproc = TRUE,
                          pairedLM = 0, mc.cores = 1, bending=TRUE,
                          fixRepro = FALSE,stepsize=0.02,
                          missingList=misslist)
dimnames(slided4.all[["dataslide"]])[3]<-dimnames(newpts2)[3]
#save(slided4.all,file="~/Google Drive/NHM/crocs/data/slid.crocs.all.apr25.R")
#load("~/Google Drive/NHM/crocs/data/slid.crocs.sept23.R")
}



#re-estimate missing post sliding
slided4.all$dataslide[which(is.na(newpts))]<-NA
#Fix missing landmarks
slid.lms<-fixLMtps(slided4.all$dataslide)
#the landmarks ready for mirroring are in an object called slid.lms$out


slidedlms <- slid.lms$out

checkLM(slidedlms,path="./ply/", pt.size = 2,suffix=".ply",render="s", begin = 5)

save(slidedlms, file = 'Final_dataset_LMs_102.R')


load('slidedlms.R')

#Pull out right hand side to reflect 

RHS <- slidedlms[,,c(5, 15, 16, 17, 20)]
LHS <- slidedlms[,,-c(5, 15, 16, 17, 20)]

#Reflect RHS to left 
arr.refl <-(-1)*(RHS[1:length(1:nrow(RHS)),,])

#Bind arr.refl with LHS 
#Bind the 2 datasets 
binded_data=abind::abind(LHS, arr.refl, along = 3)

#pull out the canines 

spheres3d(binded_data[c(7, 9, 11, 71:90),, 1])
spheres3d(binded_data[c(7, 9, 11),, 1], col = 'red')
spheres3d(binded_data[c(71:90),, 1], col = 'black')


#pull out the canines 
canines <- binded_data[c(7, 9, 11, 71:90),,]

#remove the blank specimen 
canines_no_blank <- canines[,,-c(1)]


#Procrustes the whole skull
Y.gpa=gpagen(canines_no_blank) #Procrustes
canines_no_blank_proc=Y.gpa$coords #Subset out the coords 
size=Y.gpa$Csize #Centroid size

plot(canines_no_blank_proc) #raw data


#RUN PCA ON PROCRUSTED DATA WITH MIRRORED LMS REMOVED
#This runs a PC and pulls out PC scores (you can also use this for the mesh)
PCA.ALL <- gm.prcomp(canines_no_blank_proc)
plot(PCA.ALL) #raw morphospace

#to look at specimen position (numbers)
text(PCA.ALL$x[,1],PCA.ALL$x[,2], pos = 4) #these numbers are the same as dimnames(arraylm)[3]


#Pull out the max and min shapes 
spheres3d(PCA.ALL$shapes$shapes.comp1$min, radius=.02,color = "grey")
spheres3d(PCA.ALL$shapes$shapes.comp1$max, radius=.02,color = "grey")
spheres3d(PCA.ALL$shapes$shapes.comp2$min, radius=.02,color = "grey")
spheres3d(PCA.ALL$shapes$shapes.comp2$max, radius=.02,color = "grey")


#Building a nicer morphospace with ecological data 
#This is ordered the same as the binned data 
species_data <- read.csv('species_teeth_binded_no_blanks.csv', fileEncoding="UTF-8-BOM")

library(hablar)
library(ggrepel)
library(ggConvexHull)

pca_scores <- as_tibble(PCA.ALL$x)
pca_scores <- pca_scores %>% mutate(., full=rownames(PCA.ALL$x))

pca_scores2 <- left_join(species_data, pca_scores, by = c("full"))
pca_scores2 <- pca_scores2 %>% convert(fct(genetic_id)) #or whichever factor you are keen to look at 


pca_res <- PCA.ALL #the Procrusted, aligned data 

#Set the palette
cbPalette <- c("#481567FF", "#20A387FF", "#DCE319FF")

g <- ggplot(pca_scores2, aes(x=Comp1, y=Comp2, color=genetic_id))+
  geom_point()
g <-ggplot(pca_scores2, aes(x=Comp1, y=Comp2, fill=genetic_id))+ 
  
  geom_point(size=8,aes(colour=genetic_id, fill= genetic_id))+
  theme_bw()+
  theme(aspect.ratio = 1) +
  scale_alpha(range=c(0.4, 1),guide='none')+
  #ggtitle("Principal Components Analysis") +
  
  xlab(paste0("PC Axis 1 (", signif((pca_res$pc.summary$importance[2,1]*100),3), "% of Total Variance)")) + 
  ylab(paste0("PC Axis 2 (",signif((pca_res$pc.summary$importance[2,2]*100),3),"% of Total Variance)")) +
  geom_text_repel(aes(label = species_data$number), 
                  box.padding = unit(0.45, "lines"),  max.overlaps = 30) + #uncomment this if you want names on the plot
  labs(fill = "genetic_id") + 
  scale_color_manual(values= cbPalette) +
  #geom_convexhull(alpha = 0.3,aes(fill = genetic_id))
  
  
  g #plot

g + scale_y_reverse()
g

#To add specimen names 
g2<- g + geom_text_repel(inherit.aes = FALSE, aes(fontface=1,x=Comp1, y=Comp2,label=sub("_", " ", full), size = 1 ))
g2


#Pull out M3

#pull out the canines 
atarfa=ply2mesh(file="D:/Missy lemurs/Elizabeths_ PCA and lemur stuff/landmarks/ply/ply ASCII/ply/Propithecus diadema USNM_63347.ply")
shade3d(atarfa,col='white')
spheres3d(binded_data[c(7, 9, 11, 71:90),, 14], radius = 0.7, col = 'red')

M3 <- (binded_data[c(60:69),,])
#remove the blank specimen 
M3_no_blank <- M3[,,-c(1)]

#Procrustes the whole skull
Y.gpa=gpagen(M3_no_blank) #Procrustes
M3_no_blank_proc=Y.gpa$coords #Subset out the coords 
size=Y.gpa$Csize #Centroid size

plot(M3_no_blank_proc) #raw data


#RUN PCA ON PROCRUSTED DATA WITH MIRRORED LMS REMOVED
#This runs a PC and pulls out PC scores (you can also use this for the mesh)
PCA.ALL <- gm.prcomp(M3_no_blank_proc)
plot(PCA.ALL) #raw morphospace

#to look at specimen position (numbers)
text(PCA.ALL$x[,1],PCA.ALL$x[,2], pos = 4) #these numbers are the same as dimnames(arraylm)[3]


#Pull out the max and min shapes 
spheres3d(PCA.ALL$shapes$shapes.comp1$min, radius=.06,color = "grey")
spheres3d(PCA.ALL$shapes$shapes.comp1$max, radius=.06,color = "grey")
spheres3d(PCA.ALL$shapes$shapes.comp2$min, radius=.06,color = "grey")
spheres3d(PCA.ALL$shapes$shapes.comp2$max, radius=.06,color = "grey")


#Building a nicer morphospace with ecological data 
#This is ordered the same as the binned data 
species_data <- read.csv('species_teeth_binded_no_blanks.csv', fileEncoding="UTF-8-BOM")

library(hablar)
library(ggrepel)
library(ggConvexHull)

pca_scores <- as_tibble(PCA.ALL$x)
pca_scores <- pca_scores %>% mutate(., full=rownames(PCA.ALL$x))

pca_scores2 <- left_join(species_data, pca_scores, by = c("full"))
pca_scores2 <- pca_scores2 %>% convert(fct(genetic_id)) #or whichever factor you are keen to look at 


pca_res <- PCA.ALL #the Procrusted, aligned data 

#Set the palette
cbPalette <- c('grey', "#481567FF", "#20A387FF", "#DCE319FF")

g <- ggplot(pca_scores2, aes(x=Comp1, y=Comp2, color=genetic_id))+
  geom_point()
g <-ggplot(pca_scores2, aes(x=Comp1, y=Comp2, fill=genetic_id))+ 
  
  geom_point(size=6,aes(colour=genetic_id, fill= genetic_id))+
  theme_bw()+
  theme(aspect.ratio = 1) +
  scale_alpha(range=c(0.4, 1),guide='none')+
  #ggtitle("Principal Components Analysis") +
  
  xlab(paste0("PC Axis 1 (", signif((pca_res$pc.summary$importance[2,1]*100),3), "% of Total Variance)")) + 
  ylab(paste0("PC Axis 2 (",signif((pca_res$pc.summary$importance[2,2]*100),3),"% of Total Variance)")) +
  geom_text_repel(aes(label = species_data$number), 
                  box.padding = unit(0.45, "lines"),  max.overlaps = 30) + #uncomment this if you want names on the plot
  labs(fill = "genetic_id") + 
  scale_color_manual(values= cbPalette) +
  #geom_convexhull(alpha = 0.3,aes(fill = genetic_id))
  
  
  g #plot

g + scale_y_reverse()
g

#To add specimen names 
g2<- g + geom_text_repel(inherit.aes = FALSE, aes(fontface=1,x=Comp1, y=Comp2,label=sub("_", " ", full), size = 1 ))
g2



#MANOVAs

gdf <- geomorph.data.frame(Y.gpa, genetic_id = species_data$genetic_id, sex = species_data$sex, age = species_data$age, Csize = size)

#Allometry on coords (Y.gpa holds all coords Csize etc)
fit_size <- procD.lm(coords ~ log(Csize), data=gdf, iter=9999, print.progress = FALSE)
summary(fit_size) 


#ecological factors 
fit <- procD.lm(coords ~ genetic_id, iter=9999, data = gdf)
summary(fit)


#Missy's curves 
#Pull out Missy's original curves 
shape_data_lhs_MH <- newarray[c(1:26,32:41, 62:76, 77:91),,]
spheres3d(shape_data_lhs_MH[,,11], radius = 1)#Check 

#Run through above

spheres3d(shape_data_lhs_MH[,,13], radius = 0.005, col = 'blue')
spheres3d(shape_data_lhs_MH[,,3], radius = 0.005, col = 'pink')


spheres3d(shape_data[c(1:101),,13], radius = 0.002, col = 'blue') #LHS
spheres3d(shape_data[c(102:210),,13], radius = 0.002, col = 'pink') #RHS

spheres3d(shape_data_lhs[c(1:26),,13], radius = 0.002, col = 'blue') #LMs
spheres3d(shape_data_lhs[c(27:101),,13], radius = 0.002, col = 'pink') #curves

spheres3d(shape_data_lhs[c(1:26),,13], radius = 0.002, col = 'blue')
spheres3d(shape_data_lhs[c(32:41, 62:76, 77:91),,13], radius = 0.002, col = 'pink')
