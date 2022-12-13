
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
arraylm<-ptsarray
arraylm


##### CHANGING MISSING LANDMARKS BEFORE MIRRORING #########
#Otherwise you get weird numbers that aren't 9999 mirroring
#estimate.missing from Geomorph is used for estimating missing landmarks 

arraylm[which(arraylm==9999)] <- NA
arraylm <- estimate.missing(arraylm,method="TPS")


#MIRROR THESE LANDMARKS over the central line plane of the skull 

########## SYMMETRISATION TO IMPROVE THE SHAPE ANALYSES #########################
#Ellen's code 

midline<-as.integer(c(1:9)) # LM that are on the midline + parasphenoid curve points + NO patch point
#got length(midline)= 9 points on the midline

left.lm <- c(10:26)
#exclude midline points. Last number = last number of newpts 

lengmatrice=dim(arraylm)[1]*2-length(midline)#-length(nasalfrontal) #should be the length with the both sides, 1 is the column and 2 
#just means that we are duplicating the data to be on both sides of the skull 

Matrice=array(NA,dim = c(lengmatrice,3,36)) #3 is the dimensions (x, y, z), 2 is specimen number 
Matrice[1:dim(arraylm)[1],,]=arraylm

#left.lm <- c(1:37,39,41:47,50,52,53,57:60,62:66)
#left.lm <- c(2,3,5:18,21:37,39,41:47,50,52,53,57:60,62:66)
#exclude midline points. Last number = last number of newpts 

#Check left.lm and midline [left.lm,,x] = species number
spheres3d(arraylm[left.lm,,10],radius=0.5) #left LMs
spheres3d(arraylm[midline,,10],radius=0.5,col='red') #midline

right.lm <- c(27:43) #left.lm +1:lenmatrice

bilat.landmarks <- cbind(left.lm, right.lm) #one column is the rows of the right side LM and the other column the rows of the left side

MirroredAC=mirrorfill(A=Matrice,  l1=midline, l2=bilat.landmarks) # the NA rows are now filled so that the numbers are the same on both
#sides of the skull 
MirroredAC
#deformGrid3d(MirroredAC[67:123,,2], Matrice[,,2], ngrid=0) #This shows you the new mirroed landmarks 

#These visualisations are done before Procrusted data 
#This shows the original landmarks
spheres3d(MirroredAC[c(27:43),,1],col='green',radius=0.75)
spheres3d(MirroredAC[c(10:26),,1],col='red',radius=0.75)
spheres3d(MirroredAC[c(1:9),,1],col='black',radius= 0.75)
#check dimensions


#Procrustes the data 
arraylm_proc <- MirroredAC
Y.gpa=gpagen(arraylm_proc) #Procrustes
data=Y.gpa$coords #Subset out the coords 
size=Y.gpa$Csize #centroid size

plot(data)

#PCA
PCA.ALL <- gm.prcomp(data)
plot(PCA.ALL)

#to look at specimen position 
text(PCA.ALL$x[,1],PCA.ALL$x[,2], pos = 4)

