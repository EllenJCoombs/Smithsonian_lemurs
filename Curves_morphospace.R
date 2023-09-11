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
midline<-as.integer(c(1:8))
left.curves<-c(1:7)
left.lm <- c(10:26)
right.lm <- c(27:43)
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


#Double check the extra landmarks
open3d();
spheres3d(newarray[,,13],radius=1, col = 'green')
spheres3d(newarray[,,14],radius=1, col = 'brown')
spheres3d(newarray[bilats[,1],,11],col='red',radius=1)
spheres3d(newarray[bilats[,2],,11],col='blue',radius=1)


#Remove the extra landamrks after Procrutes 

#Procrustes the data (not phylo correction, see above)
Y.gpa=gpagen(newarray) #Procrustes
shape_data=Y.gpa$coords #Subset out the coords 
size=Y.gpa$Csize #Centroid size

plot(shape_data) #raw data

#Remove additional lms 
shape_data_lhs <- shape_data[c(1:101),,]

#This runs a PC and pulls out PC scores (you can also use this for the mesh)
PCA.ALL <- gm.prcomp(shape_data_lhs)
plot(PCA.ALL) #raw morphospace

#to look at specimen position (numbers)
text(PCA.ALL$x[,1],PCA.ALL$x[,2], pos = 4) #these numbers are the same as dimnames(arraylm)[3]


spheres3d(shape_data[,,13], radius = 0.002, col = 'blue')
spheres3d(shape_data[,,3], radius = 0.002, col = 'pink')


spheres3d(shape_data[c(1:101),,13], radius = 0.002, col = 'blue') #LHS
spheres3d(shape_data[c(102:210),,13], radius = 0.002, col = 'pink') #RHS

spheres3d(shape_data_lhs[c(1:26),,13], radius = 0.002, col = 'blue') #LMs
spheres3d(shape_data_lhs[c(27:101),,13], radius = 0.002, col = 'pink') #curves

spheres3d(shape_data_lhs[,,13], radius = 0.002, col = 'blue')
spheres3d(shape_data_lhs[,,22], radius = 0.002, col = 'pink')


#Pull out the max and min shapes 
spheres3d(PCA.ALL$shapes$shapes.comp1$min, radius=.002,color = "grey")
spheres3d(PCA.ALL$shapes$shapes.comp1$max, radius=.002,color = "brown")
spheres3d(PCA.ALL$shapes$shapes.comp2$min, radius=.002,color = "grey")
spheres3d(PCA.ALL$shapes$shapes.comp2$max, radius=.002,color = "brown")

