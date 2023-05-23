
#load your libraries
library(Rvcg)
library(Morpho)
library(rgl)
library(geomorph)
library(paleomorph)
library(geomorph)


#Read in AP's landmarks (landmark data set number 1)
landmarks <- geomorph::read.morphologika('AP_lemurs_LHS.txt') #this is how to read in morphologika data, see below for .pts from checkpoint

arraylm <- landmarks #first data set

#Check LM numbers match your desciptions and are placed where you would expect
spheres3d(arraylm[,,2]) #number denotes specimen number
text3d(arraylm[c(1:26),,3], text = 1:52) #This numbers the landmarks in the order they were placed

#There are NAs - sort with as.missing 
landmarks[which(landmarks==9999)] <- NA
landmarks <- estimate.missing(landmarks,method="TPS")


#RHS only - or use numbers (below)

RHS_only <- landmarks[,,c("Propithecus diadema USNM_63350", "Propithecus diadema MCZ_BOM-5016", "Propithecus verreauxi MCZ_16719",
                          "Propithecus coronatus MCZ_8040", "Propithecus coquereli MCZ_8041", "Propithecus edwardsi MCZ_8042", "Propithecus verreauxi MCZ_44873",
                          "Propithecus verreauxi MCZ_44875", "Propithecus verreauxi MCZ_44876", "Propithecus diadema MCZ_18741",
                          "Propithecus verreauxi MCZ_16452", "Propithecus verreauxi MCZ_16375", "Propithecus verreauxi MCZ_16377", "Propithecus verreauxi MCZ_16379", 
                          "Propithecus verreauxi MCZ_16380", "Propithecus verreauxi MCZ_16381", "Propithecus verreauxi MCZ_16383", "Propithecus diadema MCZ_44853", "Propithecus diadema MCZ_44854", 
                          "Propithecus diadema MCZ_44855", "Propithecus verreauxi MCZ_44860", )]


RHS_only <- landmarks[,,c(35, 39, 40, 41, 42, 43, 44, 45, 46, 47, 
                          48, 49, 51, 52, 53, 54, 55, 56, 57, 58, 60, 
                          61, 62, 64, 65, 66, 67, 68, 69, 70, 71, 72, 105)]

#for mirroring
arraylm <- RHS_only

#Mirror the RHS 

########## SYMMETRISATION TO IMPROVE THE SHAPE ANALYSES #########################
#Ellen's code 

midline<-as.integer(c(1:9)) # LM that are on the midline + parasphenoid curve points + NO patch point
#got length(midline)= 9 points on the midline

right.lm <- c(10:26)
#exclude midline points. Last number = last number of newpts 

lengmatrice=dim(arraylm)[1]*2-length(midline)#-length(nasalfrontal) #should be the length with the both sides, 1 is the column and 2 
#just means that we are duplicating the data to be on both sides of the skull 

Matrice=array(NA,dim = c(lengmatrice,3,33)) #3 is the dimensions (x, y, z), 2 is specimen number 
Matrice[1:dim(arraylm)[1],,]=arraylm

#left.lm <- c(1:37,39,41:47,50,52,53,57:60,62:66)
#left.lm <- c(2,3,5:18,21:37,39,41:47,50,52,53,57:60,62:66)
#exclude midline points. Last number = last number of newpts 

#Check left.lm and midline [left.lm,,x] = species number
spheres3d(arraylm[right.lm,,1],radius=0.5) #left LMs
spheres3d(arraylm[midline,,1],radius=0.5,col='red') #midline

left.lm <- c(27:43) #left.lm +1:lenmatrice

bilat.landmarks <- cbind(right.lm, left.lm) #one column is the rows of the right side LM and the other column the rows of the left side

MirroredAC=mirrorfill(A=Matrice,  l1=midline, l2=bilat.landmarks) # the NA rows are now filled so that the numbers are the same on both
#sides of the skull 
MirroredAC

spheres3d(MirroredAC[right.lm,,10],radius=1) #left LMs
spheres3d(MirroredAC[midline,,10],radius=1,col='red') #midli
spheres3d(MirroredAC[left.lm,,10],radius=1,col='green') 


#add names
dimnames(MirroredAC)[3]<-list(
  substr(dir("./pts",pattern=".pts"),1,(nchar(dir("./pts",pattern=".pts"))-4)))


#Procrustes and then remove the RHS
mirrored_PROC=geomorph::gpagen(MirroredAC) #Remove non-shape aspects 
COORDS=mirrored_PROC$coords #Subset

#Pull out the LHS only 
EC_lemurs_LHS <- COORDS[c(1:9, 27:43),,]

#PCA - Principle components analysis
pca_res <- geomorph::gm.prcomp(EC_lemurs_LHS)


species_data <- read.csv('EC_cleaned_data.csv')


#Read in Ellen's LMs (landmark data set number 2)

###=========== LOADING DATA SET 1:LE LANDMARKED SKULL 
#Read in landmarks manually placed on the whole skull 

ntaxa<-92 ## number of specimens in the test data set - NB you also put this in the code below (x,y,z) 

#data set .pts from Checkpoint - read in 
ptslist<-dir(pattern='.pts',recursive=T)
ptsarray<-array(dim=c(26,3,92)) #dim=c(number of landmarks and semilandmarks, number of dimensions, number of specimens)
for(i in 1:length(ptslist))
{
  ptsarray[,,i]<-as.matrix(read.table(file=ptslist[i],skip=2,header=F,sep="",row.names=1))
}

#Pull the names from the .pts filesso they match the coords
#[3] stays the same 
dimnames(ptsarray)[3]<-list(
  substr(dir("./pts",pattern=".pts"),1,(nchar(dir("./pts",pattern=".pts"))-4)))

#Fix missing landmarks 
ptsarray[which(ptsarray==9999)] <- NA
ptsarray <- estimate.missing(ptsarray,method="TPS")


arraylm <- ptsarray #this is your second data set

#Have a look at the two data sets 
spheres3d(landmarks[,,1], radius = 12)
spheres3d(array2[,,2], col = 'red')

########## SYMMETRISATION TO IMPROVE THE SHAPE ANALYSES #########################
#Ellen's code 

midline<-as.integer(c(1:9)) # LM that are on the midline + parasphenoid curve points + NO patch point
#got length(midline)= 9 points on the midline

left.lm <- c(10:26)
#exclude midline points. Last number = last number of newpts 

lengmatrice=dim(arraylm)[1]*2-length(midline)#-length(nasalfrontal) #should be the length with the both sides, 1 is the column and 2 
#just means that we are duplicating the data to be on both sides of the skull 

Matrice=array(NA,dim = c(lengmatrice,3,92)) #3 is the dimensions (x, y, z), 2 is specimen number 
Matrice[1:dim(arraylm)[1],,]=arraylm

#left.lm <- c(1:37,39,41:47,50,52,53,57:60,62:66)
#left.lm <- c(2,3,5:18,21:37,39,41:47,50,52,53,57:60,62:66)
#exclude midline points. Last number = last number of newpts 

#Check left.lm and midline [left.lm,,x] = species number
spheres3d(arraylm[right.lm,,1],radius=0.5) #left LMs
spheres3d(arraylm[midline,,1],radius=0.5,col='red') #midline

right.lm <- c(27:43) #left.lm +1:lenmatrice

bilat.landmarks <- cbind(right.lm, left.lm) #one column is the rows of the right side LM and the other column the rows of the left side

MirroredAC=mirrorfill(A=Matrice,  l1=midline, l2=bilat.landmarks) # the NA rows are now filled so that the numbers are the same on both
#sides of the skull 
MirroredAC

spheres3d(MirroredAC[right.lm,,10],radius=1) #left LMs
spheres3d(MirroredAC[midline,,10],radius=1,col='red') #midli
spheres3d(MirroredAC[left.lm,,10],radius=1,col='green') 


#add names
dimnames(MirroredAC)[3]<-list(
  substr(dir("./pts",pattern=".pts"),1,(nchar(dir("./pts",pattern=".pts"))-4)))


#Procrustes and then remove the RHS
mirrored_PROC=geomorph::gpagen(MirroredAC) #Remove non-shape aspects 
COORDS=mirrored_PROC$coords #Subset

#Pull out the LHS only 
#EC_lemurs_RHS <- MirroredAC[c(1:9, 27:43),,]
EC_lemurs_LHS <- COORDS[c(1:26),,]



species_data <- read.csv('EC_cleaned_data.csv')

View(cbind(dimnames(EC_lemurs_LHS)[[3]], species_data$full))

#PCA - Principle components analysis
pca_res <- geomorph::gm.prcomp(EC_lemurs_LHS)

plot(pca_res, axis1 = 1, axis2 = 2) #to plot 
#to look at specimen position 
text(pca_res$x[,1],pca_res$x[,2], pos = 4)


#Building a nicer morphospace with ecological data

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

species_data <- read.csv('EC_cleaned_data.csv', fileEncoding="UTF-8-BOM") #this is your species excel csv with the new coloum you added 

pca_scores <- as_tibble(pca_res$x)   #this is your PCA data as before 
pca_scores <- pca_scores %>% mutate(., full=rownames(pca_res$x))

pca_scores2 <- left_join(species_data, pca_scores, by = c("full")) 
pca_scores2 <- pca_scores2 %>% convert(fct(age)) #pulling out landmarking_side as the factor (could also pull out diet, sex, age, etc)

#pca_res <- PCA.ALL #the Procrusted, aligned data

#Set the palette
cbPalette <- c("#481567FF", "#20A387FF", "green", "pink", "black", "blue", "yellow", "brown", "orange",
               'red', 'grey') #You only need two colours, for 'left' and 'right' 

g <- ggplot(pca_scores2, aes(x=Comp1, y=Comp2, color=age))+
  geom_point()
g <-ggplot(pca_scores2, aes(x=Comp1, y=Comp2, fill=genetic_id))+
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




species_mesh <- vcgImport( paste( "D:/Missy lemurs/Diadema/plys/ply ASCII/32_HawkinsMCZ44853_Pdiadema.ply", sep = "" ) )


species_mesh #that you already defined to vizualise the landmarks
Species_landmark <-EC_lemurs_LHS[,,3]#same
#PC min and max
PC1max<-pca_res$shapes$shapes.comp1$max
PC1min<-pca_res$shapes$shapes.comp1$min
PC2max<-pca_res$shapes$shapes.comp2$max
PC2min<-pca_res$shapes$shapes.comp2$min

colors()
# max PC1

warpskull <- tps3d(species_mesh,Species_landmark, PC1max)
shade3d(warpskull, col='white', specular=1)

warpskull <- tps3d(species_mesh,Species_landmark, PC1min)
shade3d(warpskull, col='white', specular=1)

warpskull <- tps3d(species_mesh,Species_landmark, PC2max)
shade3d(warpskull, col='white', specular=1)

warpskull <- tps3d(species_mesh,Species_landmark, PC2min)
shade3d(warpskull, col='white', specular=1)



spheres3d(pca_res$shapes$shapes.comp1$min, radius=.008,color = "grey") #PC1 min (left side of x axis)
spheres3d(pca_res$shapes$shapes.comp1$max, radius=.01,color = "grey") #PC2 max (right side of x axis)
spheres3d(pca_res$shapes$shapes.comp2$min, radius=.008,color = "grey") #PC2 min (bottom of y axis)
spheres3d(pca_res$shapes$shapes.comp2$max, radius=.01,color = "grey") #PC2 max (top of y axis)


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



g <- ggplot(pca_scores2, aes(x=Comp1, y=Comp2, color=species)) +
  geom_point(size=4) +
  stat_chull(aes(fill=species), alpha=9) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  scale_alpha(range=c(0.4, 1), guide='none') +
  xlab(paste0("PC Axis 1 (", signif((pca_res$pc.summary$importance[2,1]*100),3), "% of Total Variance)")) +
  ylab(paste0("PC Axis 2 (", signif((pca_res$pc.summary$importance[2,2]*100),3), "% of Total Variance)")) +
  labs(fill = "Species") +
  scale_color_manual(values= cbPalette)

g



specimen_names <- dimnames(RHS_only)[[3]]


# Add the specimen names to the third dimension of the new array
dimnames(AP_lemurs_RHS)[[3]] <- specimen_names




