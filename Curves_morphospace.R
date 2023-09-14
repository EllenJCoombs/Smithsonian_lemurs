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
ptsfolder <- "D:/Missy lemurs/Diadema/w curves/pts"


#import the pts file names
ptslist <- dir(ptsfolder, pattern='.pts', recursive=F)

my_curves <- create_curve_info(curve_table, n_fixed = 26)
setwd(ptsfolder)
subsampled.lm <- import_chkpt_data(ptslist, my_curves, subsampl = TRUE, verbose=TRUE)

#capture the output if it's too big to read in the console 
capture.output(import_chkpt_data(ptslist, my_curves, subsampl = TRUE, verbose=TRUE), file = "odonts.txt", append = T)
#if you have any missing points, Checkpoint will set them to x,y,z=9999
#this makes for bad plots in checkLM, so switch them to NA
subsampled.lm[subsampled.lm == 9999] <- NA
subsampled.lm[subsampled.lm == Inf] <- NA
subsampled.lm[subsampled.lm == -Inf] <- NA


#SET WORKING DIRECTORY TO ASCII PLY
#check to make sure your curves look okay on each specimen
checkLM(subsampled.lm,path="./ply/", pt.size = 2,suffix=".ply",render="s", begin = 5)


newpts <- subsampled.lm

#Create missing list 
misslist<-createMissingList(dim(newpts)[3])
for (j in 1:dim(newpts)[[3]]){
  misslist[[j]]<-which(is.na(newpts[,1,j]))
} 
newpts2<-fixLMtps(newpts)
#check that your pts and ply match
plyfolder <- "D:/Missy lemurs/Diadema/w curves/plys"
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


setwd("D:/Missy lemurs/Diadema/w curves/ply")
slided4.all <- slider3d_2(newpts2$out, SMvector= my_curves$Sliding.LMs,
                          outlines = my_curves$Curves, sur.path = "./ply", sur.name = NULL, 
                          meshlist = paste("./ply/",dimnames(newpts2$out)[[3]],".ply",sep=""), ignore = NULL,
                          sur.type = "ply", tol = 1e-10, deselect = FALSE, inc.check = FALSE,
                          recursive = TRUE, iterations = 3, initproc = TRUE,
                          pairedLM = 0, mc.cores = 1, bending=TRUE,
                          fixRepro = FALSE,stepsize=0.02,
                          missingList=misslist)
dimnames(slided4.all[["dataslide"]])[3]<-dimnames(newpts2$out)[3]
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


#######################################################
#                                                     #
#      MIRRORING SYMMETRIC SPECIMENS LANDMARKS        #
#                                                     #
#######################################################

#Code for mirroring easy symmetrical specimens 

## REMEMBER THE MIDLINE IS DIFFERENT ##
midline<-as.integer(c(1:9, 27:58))
left.curves<-c(5:7) #62-99
left.lm <- c(10:26)
right.lm <- c(59:75)
left.curve.list<-unlist(my_curves$Curve.in[left.curves])
leftside<-c(left.lm,left.curve.list) # RHS LMs+ RHS curves+all patch points
num.missing<-(length(leftside)+length(right.lm)) # number of LMs to create= total RHS-current LHS LMs
blanks<-c((dim(slidedlms)[1]+1):(dim(slidedlms)[1]+num.missing))
# to fill in blanks from one row past the last current point, for the number of rows needed (num.missing)
rightside<-c(right.lm,blanks)
add_col_or_row = function(x, n = 1, add_col = T, fill = 0)
{
  m1 = matrix(x, ncol = if(add_col) nrow(x) * ncol(x) else nrow(x), byrow = T)
  m2 = matrix(fill, nrow = if(add_col) dim(x)[3] else prod(dim(x)[-1]),
              ncol = if(add_col) nrow(x) * n else n)
  array(t(cbind(m1, m2)),
        c(nrow(x) + ((!add_col) * n), ncol(x) + (add_col * n), dim(x)[3]))
}
specimens2<-add_col_or_row(slidedlms,n=num.missing,add_col=FALSE,fill=NA)
dimnames(specimens2)[3]<-dimnames(slidedlms)[3]
bilats<-cbind(leftside,rightside)
newarray<-mirrorfill(specimens2,l1=midline,l2=bilats)
dimnames(newarray)[3]<-dimnames(slidedlms)[3]


#Remove duplicates
newarray <- newarray[-c(159:175),,]



#Double check the extra landmarks
open3d();
spheres3d(newarray[,,14],radius=1, col = 'green')
spheres3d(newarray[,,13],radius=1, col = 'brown')
spheres3d(newarray[bilats[,1],,13],col='red',radius=1)
spheres3d(newarray[bilats[,2],,13],col='blue',radius=1)

#double check it's mirrored correctly 
spheres3d(newarray[c(1:9, 27:61),,16], radius = 2, col = 'yellow') #midline 
spheres3d(newarray[c(10:26),,16], radius = 2, col = 'blue') #LHS landmarks 
spheres3d(newarray[c(62:101),,16], radius = 2, col = 'pink') #LHS curves 
spheres3d(newarray[c(102:141),,16], radius = 2, col = 'black') #RHS curves 
spheres3d(newarray[c(142:158),,16], radius = 2, col = 'gold') #RHS landmarks#


#Pull out all LMs only (for LHS and RHS) and check - want to include the midline too 
LMs_only <- newarray[c(1:9, 10:26, 142:158),,]
spheres3d(LMs_only[c(1:9),,11], radius = 2, col = 'black')
spheres3d(LMs_only[c(10:26),,11], radius = 2, col = 'red')
spheres3d(LMs_only[c(27:43),,11], radius = 2, col = 'blue')


#Procrustes the whole skull
Y.gpa=gpagen(LMs_only) #Procrustes
LMs_proc=Y.gpa$coords #Subset out the coords 
size=Y.gpa$Csize #Centroid size

plot(LMs_proc) #raw data


#Remove the extra landmarks from the raw data - leave the midline 
#LMs only 
LMs_only_lhs_midline <- LMs_proc[c(1:26),,]
spheres3d(LMs_only_lhs_midline[,,11], radius = 0.009)


#RUN PCA ON PROCRUSTED DATA WITH MIRRORED LMS REMOVED
#This runs a PC and pulls out PC scores (you can also use this for the mesh)
PCA.ALL <- gm.prcomp(LMs_only_lhs_midline)
plot(PCA.ALL) #raw morphospace

#to look at specimen position (numbers)
text(PCA.ALL$x[,1],PCA.ALL$x[,2], pos = 4) #these numbers are the same as dimnames(arraylm)[3]


#Pull out the max and min shapes 
spheres3d(PCA.ALL$shapes$shapes.comp1$min, radius=.007,color = "grey")
spheres3d(PCA.ALL$shapes$shapes.comp1$max, radius=.007,color = "grey")
spheres3d(PCA.ALL$shapes$shapes.comp2$min, radius=.007,color = "grey")
spheres3d(PCA.ALL$shapes$shapes.comp2$max, radius=.007,color = "grey")


#Building a nicer morphospace with ecological data 

species_data <- read.csv('species_data_diadema.csv', fileEncoding="UTF-8-BOM")

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


#########################
#
#     CURVES
#


#Use 'mirrored_curves_raw' (newarray) - check
spheres3d(mirrored_curves_raw[c(1:9, 27:58),,16], radius = 2, col = 'yellow') #midline 
spheres3d(mirrored_curves_raw[c(10:26),,16], radius = 2, col = 'blue') #LHS landmarks 
spheres3d(mirrored_curves_raw[c(62:101),,16], radius = 2, col = 'pink') #LHS curves 
spheres3d(mirrored_curves_raw[c(102:141),,16], radius = 2, col = 'black') #RHS curves 
spheres3d(mirrored_curves_raw[c(142:158),,16], radius = 2, col = 'gold') #RHS LMs

#Procrustes the whole data set (EC curves)
#Procrustes the whole skull
Y.gpa=gpagen(mirrored_curves_raw) #Procrustes
whole_EC_proc=Y.gpa$coords #Subset out the coords 
size=Y.gpa$Csize #Centroid size

plot(whole_EC_proc) #raw data


#Remove the extra landmarks from the raw data - leave the midline 
#LMs only 
EC_lhs_midline_proc <- whole_EC_proc[-c(102:158),,]
spheres3d(EC_lhs_midline_proc[,,11], radius = 0.009)


#RUN PCA ON PROCRUSTED DATA WITH MIRRORED LMS REMOVED
#This runs a PC and pulls out PC scores (you can also use this for the mesh)
PCA.ALL <- gm.prcomp(EC_lhs_midline_proc)
plot(PCA.ALL) #raw morphospace

#to look at specimen position (numbers)
text(PCA.ALL$x[,1],PCA.ALL$x[,2], pos = 4) #these numbers are the same as dimnames(arraylm)[3]


#Pull out the max and min shapes 
spheres3d(PCA.ALL$shapes$shapes.comp1$min, radius=.003,color = "grey")
spheres3d(PCA.ALL$shapes$shapes.comp1$max, radius=.003,color = "grey")
spheres3d(PCA.ALL$shapes$shapes.comp2$min, radius=.003,color = "grey")
spheres3d(PCA.ALL$shapes$shapes.comp2$max, radius=.003,color = "grey")


#Building a nicer morphospace with ecological data 

species_data <- read.csv('species_data_diadema.csv', fileEncoding="UTF-8-BOM")

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


#######################
#
#  Missy's curves
#
#



#Pull out Missy's curves 
#LMs only 
#Missy's curves 
#Pull out Missy's original curves 
whole_MH_raw <- mirrored_curves_raw[c(1:26, 32:41, 62:91, 102:131, 142:158),,]
spheres3d(whole_MH_raw[c(1:26, 97:113),,11], radius = 1, col = 'yellow')# Check LMs
spheres3d(whole_MH_raw[c(27:96),,11], radius = 1, col = 'blue') # Check curves
spheres3d(whole_MH_raw[,,11])

#Procrustes the whole data set (EC curves)
#Procrustes the whole skull
Y.gpa=gpagen(whole_MH_raw) #Procrustes
whole_MH_proc=Y.gpa$coords #Subset out the coords 
size=Y.gpa$Csize #Centroid size

plot(whole_MH_proc) #raw data

#REMOVE RHS OF MH DATA BEFORE PCA
MH_lhs_midline_proc <- whole_MH_proc[-c(67:113),,] #remove RHS
spheres3d(MH_lhs_midline_proc[,,11], radius = 0.004)

#RUN PCA ON PROCRUSTED DATA WITH MIRRORED LMS REMOVED
#This runs a PC and pulls out PC scores (you can also use this for the mesh)
PCA.ALL <- gm.prcomp(MH_lhs_midline_proc)
plot(PCA.ALL) #raw morphospace

#to look at specimen position (numbers)
text(PCA.ALL$x[,1],PCA.ALL$x[,2], pos = 4) #these numbers are the same as dimnames(arraylm)[3]


#Pull out the max and min shapes 
spheres3d(PCA.ALL$shapes$shapes.comp1$min, radius=.003,color = "grey")
spheres3d(PCA.ALL$shapes$shapes.comp1$max, radius=.003,color = "grey")
spheres3d(PCA.ALL$shapes$shapes.comp2$min, radius=.003,color = "grey")
spheres3d(PCA.ALL$shapes$shapes.comp2$max, radius=.003,color = "grey")


#Building a nicer morphospace with ecological data 

species_data <- read.csv('species_data_diadema.csv', fileEncoding="UTF-8-BOM")

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




