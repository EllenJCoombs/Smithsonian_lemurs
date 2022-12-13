


library(phytools)
library(ggplot2)
library(musculusColors)
library(gridExtra)
library(ggplot2)
library(ggrepel)


#classifer.whales is the species data
#Add species data 

classifier.lemurs <- read.csv('species_data_lemurs.csv')

# Firstly, perform Procrustes on shape data
Y.gpa <-Y.gpa$coords
PCA.3D.W <- geomorph::gm.prcomp( Y.gpa )

PCW <-data.frame(cbind(PCA.3D.W$x[,c(1:22)],classifier.lemurs,dimnames(Y.gpa)[[3]]))

###############################################################################
#### Get tree data and convert to geom_segment to create phylomorphospace #####
###############################################################################

## create matrix for use in phytools::fastAnc()
mat1 <- cbind(eval(substitute(PCW$Comp1), PCW),eval(substitute(PCW$Comp2), PCW))
mat2 <- cbind(eval(substitute(PCW$Comp1), PCW),eval(substitute(PCW$Comp3), PCW))
rownames(mat1) <- eval(substitute(PCW$dimnames.Y.gpa...3..), PCW)
rownames(mat2) <- eval(substitute(PCW$dimnames.Y.gpa...3..), PCW)
stopifnot(length(setdiff(treeMAND$tip.label, rownames(mat1))) == 0)
stopifnot(length(setdiff(treeMAND$tip.label, rownames(mat2))) == 0)

xAnc <- fastAnc(treeMAND, mat1[,1])
yAnc <- fastAnc(treeMAND, mat1[,2])
zAnc <- fastAnc(treeMAND, mat2[,2])

all_node_coords.1 <-
  data.frame(
    #put PC values for all nodes and tips in a dataframe
    #tips go first in order of tip labels, then numerical order for nodes
    x=c(mat1[treeMAND$tip.label,1], xAnc),
    y=c(mat1[treeMAND$tip.label,2], yAnc),
    nodeid=1:(treeMAND$Nnode + length(treeMAND$tip.label))
  )
all_node_coords.2 <-
  data.frame(
    #put PC values for all nodes and tips in a dataframe
    #tips go first in order of tip labels, then numerical order for nodes
    x=c(mat2[treeMAND$tip.label,1], xAnc),
    y=c(mat2[treeMAND$tip.label,2], zAnc),
    nodeid=1:(treeMAND$Nnode + length(treeMAND$tip.label))
  )

#get edge list from tree object
edges.1 <- data.frame(treeMAND$edge)
names(edges.1) <- c("node1", "node2")
#translate tip/node numbers into PC coordinates in all_node_coords dataframe
edgecoords.1 <- merge(
  merge(edges.1, all_node_coords.1, by.x="node1", by.y="nodeid"),
  all_node_coords.1, by.x="node2", by.y="nodeid")

#get edge list from tree object
edges.2 <- data.frame(treeMAND$edge)
names(edges.2) <- c("node1", "node2")
#translate tip/node numbers into PC coordinates in all_node_coords dataframe
edgecoords.2 <- merge(
  merge(edges.2, all_node_coords.2, by.x="node1", by.y="nodeid"),
  all_node_coords.2, by.x="node2", by.y="nodeid")


PhM1 <-ggplot()+ #PCs 1 and 2
  #geom_segment(data=edgecoords.1,aes(x=x.x,xend=x.y, y=y.x, yend=y.y), size=0.7, alpha=0.5) +
  geom_point(data = PCW, aes(x= Comp1, y = Comp2), colour = "black", shape = 21,  size=3.5)+
  theme_bw() +
  geom_text_repel(label= 'specimen_id', size = 1 )
  
plot(PhM1)


#Add specimen names 


g2<- ggplot(data = PCW, aes(x = Comp1, y = Comp2)) + theme_bw() + 
  geom_text_repel(aes(label = specimen_id), 
                  box.padding = unit(0.45, "lines"),
                  max.overlaps = Inf) +
  geom_point(size = 3.5, color = "#7795A2") 

g2


#Coloured by eco variable 


PhM1 <-ggplot()+ #PCs 1 and 2
  #geom_segment(data=edgecoords.1,aes(x=x.x,xend=x.y, y=y.x, yend=y.y), size=0.7, alpha=0.5) +
  geom_point(data = PCW, aes(x=Comp1, y = Comp2, colour = age, fill = age), colour = "black", shape = 21,  size=3.5)+
  scale_fill_manual(values=musculus_palette("Bmpoop", 5))+
  theme_bw()
plot(PhM1) 



PhM2 <-ggplot()+ #PCs 1 and 2
  #geom_segment(data=edgecoords.1,aes(x=x.x,xend=x.y, y=y.x, yend=y.y), size=0.7, alpha=0.5) +
  geom_point(data = PCW, aes(x=Comp1, y = Comp2, colour = sex, fill = sex), colour = "black",shape = 21,  size=3.5)+
  scale_fill_manual(values=musculus_palette("Bmpoop", 4))+
  theme_bw()
plot(PhM2) 


PhM3 <-ggplot()+ #PCs 1 and 2
  #geom_segment(data=edgecoords.1,aes(x=x.x,xend=x.y, y=y.x, yend=y.y), size=0.7, alpha=0.5) +
  geom_point(data = PCW, aes(x=Comp1, y = Comp2, colour = institution, fill = institution), colour = "black",shape = 21,  size=3.5)+
  scale_fill_manual(values=musculus_palette("Bmsurface", 6))+
  theme_bw()
plot(PhM3) 


grid.arrange(PhM1, PhM2, PhM3, ncol=2)
