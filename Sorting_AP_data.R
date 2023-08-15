#Sorting Anna's data
#Making erroneous data 'missing' and dealing with it accordingly 
#Pulling out AP and EC data separately to compare 
#Changing names in AP's data to correct them 
#Plotting AP and EC data 
#T.test comparing the distances in the data sets 

#EXTRAS - changing names within a shape data set 
#       - pulling specimens out of the .csv and the shape data 

#Adding missing lms

AP_EC_coords[26,,145] <- 9999

AP_EC_coords[which(AP_EC_coords==9999)] <- NA
fixed_missing <- estimate.missing(AP_EC_coords,method="TPS")


RHS_wrong <- fixed_missing[,,c(37,59)] #pull out the wrong RHS
fixed_missing <- fixed_missing[,,-c(37,59)] #delete them from the original data set 


#reflect the wrong specimens 
#reflect the LHS to the RHS
arr.refl <-(-1)*(RHS_wrong[1:length(1:nrow(RHS_wrong)),,])


#Bind the 2 datasets 
binded_data=abind::abind(fixed_missing, arr.refl, along = 3)


#Procrustes and then remove the RHS
binded_proc=geomorph::gpagen(binded_data) #Remove non-shape aspects 
binded_coords=binded_proc$coords #Subset out coords 
size=binded_proc$Csize #Centroid size

plot(binded_coords) #raw data

#This runs a PC and pulls out PC scores (you can also use this for the mesh)
PCA.ALL <- gm.prcomp(sorted_data_matrix)
plot(PCA.ALL) #raw morphospace

text(PCA.ALL$x[,1],PCA.ALL$x[,2], pos = 4)


species_data <- read.csv('AP_MH_ALL.csv')

#pull the phylo names from the species data 
Full_names_phylo=species_data$full


sorted_indices <- order(Full_names_phylo)
# Step 3: Reorder the names and corresponding data using the sorted indices
sorted_names <- Full_names_phylo[sorted_indices]
sorted_data_matrix <- binded_coords[, , sorted_indices]


#This runs a PC and pulls out PC scores (you can also use this for the mesh)
PCA.ALL <- gm.prcomp(sorted_data_matrix)
plot(PCA.ALL) #raw morphospace

text(PCA.ALL$x[,1],PCA.ALL$x[,2], pos = 4)


#Plot to see where mine and APs data lie 

library(geomorph)
library(rgl)
library(ggplot2)
library(hablar)
library(dplyr)
library(ggrepel)
library(ggConvexHull)


#Building a nicer morphospace with ecological data 
pca_scores <- as_tibble(PCA.ALL$x)
pca_scores <- pca_scores %>% mutate(., full=rownames(PCA.ALL$x))

pca_scores2 <- left_join(species_data, pca_scores, by = c("full"))
#pca_scores2 <- pca_scores2 %>% convert(fct(sex)) #or whichever factor you are keen to look at 


pca_res <- PCA.ALL #the Procrusted, aligned data 

#Set the palette
cbPalette <- c("#756BB1", "#E69F00", "#CC79A7", "#0072B2", '#F0E442', '#009E73', '#D55E00', '#D1D1CF')

g <- ggplot(pca_scores2, aes(x=Comp1, y=Comp2, color=LMer))+
  geom_point()
g <-ggplot(pca_scores2, aes(x=Comp1, y=Comp2))+ 
  
  geom_point(size=8,aes(colour=LMer, fill = LMer))+
  theme_bw()+
  theme(aspect.ratio = 1) +
  scale_alpha(range=c(0.4, 1),guide="none")+
  #ggtitle("Principal Components Analysis") +
  
  xlab(paste0("PC Axis 1 (", signif((pca_res$pc.summary$importance[2,1]*100),3), "% of Total Variance)")) + 
  ylab(paste0("PC Axis 2 (",signif((pca_res$pc.summary$importance[2,2]*100),3),"% of Total Variance)")) +
  
  labs(fill = "LMer") + 
  geom_convexhull(alpha = 0.3,aes(fill = LMer))+
  
  
  
  #geom_text_repel(aes(label = species_data$full), 
  #box.padding = unit(0.45, "lines"),  max.overlaps = 30) + #uncomment this if you want names on the plot
  scale_color_manual(values= cbPalette)

g #plot


#To add specimen names to a seperate plot
#With max.overlaps (so names don't overlap)
ggplot(data = pca_scores2, aes(x = Comp1, y = Comp2)) + theme_bw() + 
  geom_text_repel(aes(label = species_data$full), 
                  box.padding = unit(0.45, "lines"),  max.overlaps = 200) +
  geom_point(colour = 'grey', size = 3)



#Pulling out all AP and EC data to look at the morphospace

# Using grep to filter rows with "_AP" or "_EC" at the end of the specimen names
selected_rows <- grepl("_AP$|_EC$", pca_scores2$full)

# Creating a new data set with the selected rows
new_dataset <- pca_scores2[selected_rows, ]


g <- ggplot(new_dataset, aes(x=Comp1, y=Comp2, color=LMer))+
  geom_point()
g <-ggplot(new_dataset, aes(x=Comp1, y=Comp2))+ 
  
  geom_point(size=8,aes(colour=LMer, fill = LMer))+
  theme_bw()+
  theme(aspect.ratio = 1) +
  scale_alpha(range=c(0.4, 1),guide="none")+
  #ggtitle("Principal Components Analysis") +
  
  xlab(paste0("PC Axis 1 (", signif((pca_res$pc.summary$importance[2,1]*100),3), "% of Total Variance)")) + 
  ylab(paste0("PC Axis 2 (",signif((pca_res$pc.summary$importance[2,2]*100),3),"% of Total Variance)")) +
  
  labs(fill = "LMer") + 
  geom_convexhull(alpha = 0.3,aes(fill = LMer))+
  
  
  
  #geom_text_repel(aes(label = species_data$full), 
  #box.padding = unit(0.45, "lines"),  max.overlaps = 30) + #uncomment this if you want names on the plot
  scale_color_manual(values= cbPalette)

g #plot



#To add specimen names to a seperate plot
#With max.overlaps (so names don't overlap)
ggplot(data = new_dataset, aes(x = Comp1, y = Comp2)) + theme_bw() + 
  geom_text_repel(aes(label = new_dataset$full), 
                  box.padding = unit(0.45, "lines"),  max.overlaps = 200) +
  geom_point(colour = 'grey', size = 3)



#Pulling out all AP and EC data to run a t-test 

# Using grep to filter rows with "_AP" or "_EC" at the end of the specimen names
# Assuming your data is a 3-dimensional array with dimensions 26 x 3 x 303
# and you want to pull names from dimnames(data)[[3]]

# Step 1: Get the names from the third dimension of the array
names_list <- dimnames(sorted_data_matrix)[[3]]

# Step 2: Use grep to find the indices of names ending in "_AP" or "_EC"
selected_indices <- grep("_AP$", names_list)

# Step 3: Subset the array to keep only the selected names and all dimensions
new_data_AP <- sorted_data_matrix[, , selected_indices]

dimnames(new_data_AP)


# Step 1: Get the names from the third dimension of the array
names_list <- dimnames(sorted_data_matrix)[[3]]

# Step 2: Use grep to find the indices of names ending in "_AP" or "_EC"
selected_indices <- grep("_EC$", names_list)

# Step 3: Subset the array to keep only the selected names and all dimensions
new_data_EC <- sorted_data_matrix[, , selected_indices]

dimnames(new_data_EC)

#To REMOVE specimens with _EC on the end 
#Use grep to find the indices of names not ending in "_EC"
selected_indices <- which(!grepl("_EC$", names_list))

#Name changing some of AP's specimens - old name first then new name 

species_dim <- dimnames(new_data_AP)[[3]]

dimnames(new_data_AP)[[3]][dimnames(new_data_AP)[[3]] == "Propithecus deckenii AMNH_100554_AP"] <- "Propithecus coronatus AMNH_100554_AP"
dimnames(new_data_AP)[[3]][dimnames(new_data_AP)[[3]] == "Propithecus diadema AMNH_100556_AP"] <- "Propithecus candidus AMNH 100556_AP"
dimnames(new_data_AP)[[3]][dimnames(new_data_AP)[[3]] == "Propithecus diadema AMNH_100557_AP"] <- "Propithecus candidus AMNH_100557_AP"
dimnames(new_data_AP)[[3]][dimnames(new_data_AP)[[3]] == "Propithecus diadema MCZ_44855_AP"] <- "Propithecus coronatus MCZ_44855_AP"
dimnames(new_data_AP)[[3]][dimnames(new_data_AP)[[3]] == "Propithecus dickeni USNM_83964_AP"] <- "Propithecus coquereli USNM_83964_AP"
dimnames(new_data_AP)[[3]][dimnames(new_data_AP)[[3]] == "Propithecus verreauxi MNHN_1961-56_AP"] <- "Propithecus coquereli MNHN_1961-56_AP"


library(tidyr)
library(dplyr)

#T.test 

distances_AP <- dist(two.d.array(new_data_AP)) %>%
  broom::tidy(.) %>%
  mutate(., tog.or.sep = "AP")


distances_EC <- dist(two.d.array(new_data_EC)) %>%
  broom::tidy(.) %>%
  mutate(., tog.or.sep = "EC")


tibble2 <- tibble(Distances.AP = distances_AP$distance, Distances.EC = distances_EC$distance)


#t.test
t.test(x = tibble2$Distances.AP, y = tibble2$Distances.EC)


###############################################################
#                                                             #
#   EXTRAS! Change names in data set + match with shape data  #
#                                                             #
###############################################################

#Add AP suffix to AP's data 

# Row numbers (indices) of specimens to update
specimen_indices_to_update <- c(166, 168, 170)  # Replace with the actual indices

# Append "_AP" to the specimen names in the identified rows
dimnames(sorted_data_matrix)[[3]][specimen_indices_to_update] <- paste0(
  dimnames(sorted_data_matrix)[[3]][specimen_indices_to_update], "_AP"
)

# Print the updated specimen names for verification
print(dimnames(sorted_data_matrix)[[3]])



# Load the shape data matrix (replace this with your actual data loading code)
# For the purpose of this example, let's assume the matrix is named 'shape_data_matrix'
shape_data_matrix <- sorted_data_matrix

# Load the CSV file with specimen names
csv_data <- read.csv('Ap_data.csv')

# Extract the specimen names from the CSV
csv_specimen_names <- csv_data$full

# Find the intersection of specimen names
matching_indices <- which(dimnames(shape_data_matrix)[[3]] %in% csv_specimen_names)

# Subset the shape data matrix to include only the matched specimens
matched_shape_data <- shape_data_matrix[, , matching_indices]

# Perform additional processing on the matched data if needed
# ...

# Print the dimensions of the matched data
dim(matched_shape_data)

# Renamed to AP_final

#Procrustes the two data sets together (important!)

AP_PROC=geomorph::gpagen(AP_final) #Remove non-shape aspects 
COORDS=AP_PROC$coords #Subset out the coords 

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



species_data <- AP_data #this is your species excel csv with the new coloum you added 

pca_scores <- as_tibble(pca_res$x)   #this is your PCA data as before 
pca_scores <- pca_scores %>% mutate(., full=rownames(pca_res$x))

pca_scores2 <- left_join(species_data, pca_scores, by = c("full")) 
pca_scores2 <- pca_scores2 %>% convert(fct('genus')) #pulling out landmarking_side as the factor (could also pull out diet, sex, age, etc)

#pca_res <- PCA.ALL #the Procrusted, aligned data

#Set the palette
cbPalette <- c("#481567FF", "#20A387FF", "darkgreen", "pink", "black", "blue", "grey", "brown", "orange",
               'red', 'grey') #You only need two colours, for 'left' and 'right' 

g <- ggplot(pca_scores2, aes(x=Comp1, y=Comp2, color=genus))+
  geom_point()
g <-ggplot(pca_scores2, aes(x=Comp1, y=Comp2, fill=genus))+
  geom_point(size=6,aes(colour=genus, fill= genus))+
  theme_bw()+
  theme(aspect.ratio = 1) +
  scale_alpha(range=c(0.4, 1),guide='none')+
  #ggtitle("Principal Components Analysis looking at landmarking side") +
  
  xlab(paste0("PC Axis 1 (", signif((pca_res$pc.summary$importance[2,1]*100),3), "% of Total Variance)")) +
  ylab(paste0("PC Axis 2 (",signif((pca_res$pc.summary$importance[2,2]*100),3),"% of Total Variance)")) +
  
  labs(fill = "genus") +
  scale_color_manual(values= cbPalette)

g #plot

#To add specimen names if you like! Always helpful to see what's where (this pulls the names from the 'full' column 
g2<- g + geom_text_repel(inherit.aes = FALSE, aes(fontface=1,x=Comp1, y=Comp2,label=sub("_", " ", label), size = 1 ))
g2






