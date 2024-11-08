
#morphometric code for Waleska
install.packages('geomorph') #example of how to install a package

#load the packages 
library(geomorph)
library(paleomorph)
library(Morpho)
library(rgl)

setwd("D:/Inia skulls/test/pts") #set to your pts folder
#NB the plys have been transformed to ASCII (from binary) so that landmarks can be visualised on the mesh 
#See code 'Binary_ASCII_ply.R" if needed for this step

ntaxa<-20 #number of specimens - NB can also put this in the code below (x,y,z) 
#data set is pulled from the .pts from Checkpoint

#load the pts data 
ptslist<-dir(pattern='.pts',recursive=T)
ptsarray<-array(dim=c(121,3,20)) #dim=c(number of landmarks and semilandmarks, number of dimensions, number of specimens)
for(i in 1:length(ptslist))
{
  ptsarray[,,i]<-as.matrix(read.table(file=ptslist[i],skip=2,header=F,sep="",row.names=1))
}

#Set to the whole project directory for this - the code wants to look for 'ply' folder 
dimnames(ptsarray)[3]<-list(
  substr(dir("./pts",pattern=".pts"),1,(nchar(dir("./pts",pattern=".pts"))-4))) #pulls the specien name from the file
arraylm <- ptsarray #make a duplicate of the data and call it arraylm
arraylm #run to check how this looks

dimnames(arraylm)[3] #check the final list order by pulling out the dimension names from dimension [3] where the species names/number sits

##### CHANGING MISSING LANDMARKS #########
#estimate.missing from Geomorph is used for estimating missing landmarks 

arraylm[which(arraylm==9999)] <- NA
arraylm <- estimate.missing(arraylm,method="TPS")

#check manually without the mesh first 
spheres3d(arraylm[,,6],radius=4, col = 'red') #check lms
spheres3d(arraylm[,,3],radius=5, col = 'blue') 

text3d(arraylm[,,12], text = 1:121) #check the number of each landmark and landmark ordering with text3d

#can also look at specimens on the mesh 
atarfa=ply2mesh(file="D:/Inia skulls/test/ply/AMNH209106.ply")
shade3d(atarfa,col='white')
spheres3d(arraylm[,,6], col = 'blue', radius = 4) #add specimen number 

#this code runs through each of the speimens in turn (you have to hit return in the console to see the next specimen)
checkLM(arraylm,path="./ply/", pt.size = 2,suffix=".ply",render="s", begin = 6)

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

species_data <- read.csv('EC_cleaned_data.csv', fileEncoding="UTF-8-BOM") #this is your species excel csv with the new coloum you added 

pca_scores <- as_tibble(pca_res$x)   #this is your PCA data as before 
pca_scores <- pca_scores %>% mutate(., full=rownames(pca_res$x))

pca_scores2 <- left_join(species_data, pca_scores, by = c("full")) 
pca_scores2 <- pca_scores2 %>% convert(fct(age)) #pulling out landmarking_side as the factor (could also pull out diet, sex, age, etc)

#Set the palette
cbPalette <- c("#481567FF", "#20A387FF", "green", "pink", "black", "blue", "yellow", "brown", "orange",
               'red', 'grey') #Change depending on how many colours you need for your categories 

#plot with ggplot 
g <- ggplot(pca_scores2, aes(x=Comp1, y=Comp2, color=age))+
  geom_point()
g <-ggplot(pca_scores2, aes(x=Comp1, y=Comp2, fill=age))+
  geom_point(size=6,aes(colour=diet, fill= diet))+
  theme_bw()+
  theme(aspect.ratio = 1) +
  scale_alpha(range=c(0.4, 1),guide='none')+
  #ggtitle("Principal Components Analysis") +
  xlab(paste0("PC Axis 1 (", signif((pca_res$pc.summary$importance[2,1]*100),3), "% of Total Variance)")) +
  ylab(paste0("PC Axis 2 (",signif((pca_res$pc.summary$importance[2,2]*100),3),"% of Total Variance)")) +
  labs(fill = "genetic_id") +
  scale_color_manual(values= cbPalette)

g #plot

#To add specimen names if you like! Always helpeful to see what's where (this pulls the names from the 'full' column 
g2<- g + geom_text_repel(inherit.aes = FALSE, aes(fontface=1,x=Comp1, y=Comp2,label=sub("_", " ", label), size = 1 ))
g2

