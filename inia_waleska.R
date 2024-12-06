
#morphometric code for Waleska
install.packages('geomorph') #example of how to install a package

#load the packages 
library(geomorph)
library(paleomorph)
library(Morpho)
library(rgl)

setwd("D:/Inia skulls/pts") #set to your pts folder
#NB the plys have been transformed to ASCII (from binary) so that landmarks can be visualised on the mesh 
#See code 'Binary_ASCII_ply.R" if needed for this step


ntaxa<-40 #number of specimens - NB can also put this in the code below (x,y,z) 
#data set is pulled from the .pts from Checkpoint

#load the pts data 
ptslist<-dir(pattern='.pts',recursive=T)
ptsarray<-array(dim=c(121,3,40)) #dim=c(number of landmarks and semilandmarks, number of dimensions, number of specimens)
for(i in 1:length(ptslist))
{
  ptsarray[,,i]<-as.matrix(read.table(file=ptslist[i],skip=2,header=F,sep="",row.names=1))
}

setwd("D:/Inia skulls") 
#Set to the whole project directory for this - the code wants to look for 'ply' folder 
dimnames(ptsarray)[3]<-list(
  substr(dir("./pts",pattern=".pts"),1,(nchar(dir("./pts",pattern=".pts"))-4))) #pulls the specien name from the file
arraylm <- ptsarray #make a duplicate of the data and call it arraylm
arraylm #run to check how this looks

dimnames(arraylm)[3] #check the final list order by pulling out the dimension names from dimension [3] where the species names/number sits

##### CHANGING MISSING LANDMARKS #########
#estimate.missing from Geomorph is used for estimating missing landmarks

arraylm[which(arraylm==9999)] <- NA
arraylm[which(arraylm==Inf)] <- NA
arraylm[which(arraylm==-Inf)] <- NA
arraylm <- estimate.missing(arraylm,method="TPS")

#check manually without the mesh first 
spheres3d(arraylm[,,14],radius=4, col = 'red') #check lms
spheres3d(arraylm[,,14],radius=5, col = 'blue') 

text3d(arraylm[,,14], text = 1:121) #check the number of each landmark and landmark ordering with text3d

#can also look at specimens on the mesh 
atarfa=ply2mesh(file="D:/Inia skulls/test/ply/AMNH209106.ply")
shade3d(atarfa,col='white')
spheres3d(arraylm[,,6], col = 'blue', radius = 4) #add specimen number 

#this code runs through each of the specimens in turn (you have to hit return in the console to see the next specimen)
checkLM(arraylm,path="./ply/", pt.size = 2,suffix=".ply",render="s", begin = 6)

#If you're just reading in the curves data set, then do that here (they've been resampled and slid already)
#load('slidedlms_final.R')

#arraylm <- slidedlms #FOR CURVES ONLY 

#Procrustes the data - remove translation, rotation, and size info
arraylm_proc <- arraylm # create another copy of the data and call it something new 
Y.gpa=gpagen(arraylm_proc) #Procrustes
data=Y.gpa$coords #Subset out the coords 
size=Y.gpa$Csize #centroid size

plot(data) #check the data are distributed as you would expect 

#save the centroid size data (make sure you have set the working directory to the folder you want this to be saved in
write.csv(size, file = 'centroid_size.csv') 

#PCA - run the pca
pca.all <- gm.prcomp(data)
plot(pca.all) #simple morphospace

#to look at specimen position 
text(pca.all$x[,1],pca.all$x[,2], pos = 4)

#check what the min and max PCs look like
spheres3d(pca.all$shapes$shapes.comp1$min, radius = .002,color = "grey") #PC1 min (left side of x axis)
spheres3d(pca.all$shapes$shapes.comp1$max, radius=.002,color = "red") #PC2 max (right side of x axis)
spheres3d(pca.all$shapes$shapes.comp2$min, radius=.002,color = "grey") #PC2 min (bottom of y axis)
spheres3d(pca.all$shapes$shapes.comp2$max, radius=.002,color = "red") #PC2 max (top of y axis)

#build a nicer morphospace w/ ecological data if you have

library(ggplot2)
library(hablar)
library(tidyr)
library(dplyr)
library(ggrepel)

#Read in your ecological data spreadsheet (this should be an excel file saved as a .csv file).  
#Change the name of the species column (with the species names in) to 'full', this is how we join the shape data with the species data. 
#Make sure the name matches your .pts names 

pca_res <- pca.all #the Procrusted, aligned data

species_data <- read.csv('Morpho_Inia_data_reordered.csv', fileEncoding="UTF-8-BOM") #this is your species excel csv with the new coloum you added 

pca_scores <- as_tibble(pca_res$x)   #this is your PCA data as before 
pca_scores <- pca_scores %>% mutate(., number=rownames(pca_res$x))

pca_scores2 <- left_join(species_data, pca_scores, by = c("number")) 
#pca_scores2 <- pca_scores2 %>% convert(fct(age)) #pulling out landmarking_side as the factor (could also pull out diet, sex, age, etc)

#Set the palette
cbPalette <- c("#9566BF", "#F4E916", "#20A387FF", "#375689", "#BEBEBE") #Change depending on how many colours you need for your categories 

#plot with ggplot 
g <- ggplot(pca_scores2, aes(x=Comp1, y=Comp2, color=species))+
  geom_point()
g <-ggplot(pca_scores2, aes(x=Comp1, y=Comp2, fill=species))+
  geom_point(size=6,aes(colour=species, fill= species))+
  theme_bw()+
  theme(aspect.ratio = 1) +
  scale_alpha(range=c(0.4, 1),guide='none')+
  #ggtitle("Principal Components Analysis") +
  xlab(paste0("PC Axis 1 (", signif((pca_res$pc.summary$importance[2,1]*100),3), "% of Total Variance)")) +
  ylab(paste0("PC Axis 2 (",signif((pca_res$pc.summary$importance[2,2]*100),3),"% of Total Variance)")) +
  labs(fill = "species") +
  scale_color_manual(values= cbPalette)

g #plot
#To add specimen names if you like! Always helpeful to see what's where (this pulls the names from the 'full' column 
g2<- g + geom_text_repel(inherit.aes = FALSE, aes(fontface=1,x=Comp1, y=Comp2,label=sub("_", " ", full), size = 0.1))
g2


#Repel species names 
g2 <- g + 
  geom_text_repel(
    inherit.aes = FALSE, 
    aes(x = Comp1, y = Comp2, label = sub("_", " ", full)), 
    fontface = 1, 
    size = 2  # Adjust this value to make the text smaller; for very small, try size = 1 or 0.5
  )

g2

#add convex hull to plot 

# Install ggforce if not already installed
# install.packages("ggforce")

library(ggplot2)
library(ggforce)
library(dplyr)

# Remove rows with NA or Inf values in Comp1 or Comp2
pca_scores2_clean <- pca_scores2 %>%
  filter(is.finite(Comp1) & is.finite(Comp2))

# Create convex hull data for each species
hulls <- pca_scores2_clean %>%
  group_by(species) %>%
  slice(chull(Comp1, Comp2))

# Plot with ggplot
g <- ggplot(pca_scores2_clean, aes(x=Comp1, y=Comp2, color=species)) +
  geom_point(size=8, aes(fill=species)) +
  geom_polygon(data=hulls, aes(x=Comp1, y=Comp2, fill=species, group=species), 
               alpha=0.3, color=NA) + # Adjust alpha for transparency of hulls
  theme_bw() +
  theme(aspect.ratio = 1) +
  scale_alpha(range=c(0.4, 1), guide='none') +
  xlab(paste0("PC Axis 1 (", signif((pca_res$pc.summary$importance[2,1]*100),3), "% of Total Variance)")) +
  ylab(paste0("PC Axis 2 (", signif((pca_res$pc.summary$importance[2,2]*100),3), "% of Total Variance)")) +
  labs(fill = "species") +
  scale_color_manual(values=cbPalette) +
  scale_fill_manual(values=cbPalette) # Ensure hull colors match point colors

g # Plot


#Stats - Function performs Procrustes ANOVA with permutation procedures to assess statistical hypotheses describing patterns of shape variation and covariation for a set of Procrustes shape variables
gdf <- geomorph.data.frame(Y.gpa)
#This is looking at habitat type and coords 
gdf <- geomorph.data.frame(gdf, location = species_data$location, species = species_data$species, 
                           sex = species_data$sex, number = species_data$number)

model1 <- procD.lm(coords ~ location, iter=999, data = gdf, RRPP = TRUE)
summary( model1)
model2 <- procD.lm(coords ~ species, iter=999, data = gdf, RRPP = TRUE)
summary(model2)


#check for collinrarity between variables. If not sig, you 
colinearity_test <- chisq.test(table(gdf$species, gdf$location))


#ADDING MORPHED MESHES TO YOUR PLOT
library(Morpho)
########################getting the minimum and maximum shape variation at the extreme of each axes
####################Vizualisation of landmarks onto 


species_mesh <- vcgImport( paste( "D:/Inia skulls/Renamed ply files/ply binary/CAS16631.ply", sep = "" ) )


species_mesh #that you already defined to vizualise the landmarks
Species_landmark <-arraylm[,,13]#same
#PC min and max
PC1max<-pca.all$shapes$shapes.comp1$max
PC1min<-pca.all$shapes$shapes.comp1$min
PC2max<-pca.all$shapes$shapes.comp2$max
PC2min<-pca.all$shapes$shapes.comp2$min

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
