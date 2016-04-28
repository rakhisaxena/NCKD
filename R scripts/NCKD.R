library('igraph')
setwd("/home/stratuser/graph-algos/NetworkComparison/data")

labelnames = c("B1","W2","F1","A-","E1","M-","WC","CA","FW")
labelnames = c("CA","B1","W2","F1","A-","E1","M-")
labelnames = c("B1","W2","F1","M-","CA")

labelnames = c("W3-", "W4-","W5-","W6-","W7-")
labelnames=c("A-", "B-", "E-")
labelnames=c("W1","W2","W3","W4","W5")
labelnames=c("DC","DY","CA")
labelnames = c("W10", "W20","W30","W40","W50")
labelnames = c("F10", "F20","F30","F40","F50")
labelnames = c("foo","cel","san","kar","adj","dol","pol","les")


colors <- rainbow(length(labelnames))
numclasses = length(labelnames)



#NCKD Dendro
nckddata <- as.matrix(read.csv(file="./MetaboliteGraphs/BarabasiGraphs/AEData/newNCKD/AEGraphs_NewNodeDistnProbabilities.txt", sep=" ",header=T))
nckddata <- as.matrix(read.csv(file="./NetSimileGraphs/NS-Data-LBD/NS-Data_NodeEdgeDistanceMatrix.txt", sep=" ",header=T))
nckddata <- as.matrix(read.csv(file="./NetSimileGraphs/NS-Data_NodeEdgeArithmeticDistanceMatrix.txt", sep=",",header=T))
nckddata <- as.matrix(read.csv(file="./NetSimileGraphs/NS-Data_NSdists.txt", sep=",",header=T))
nckddata <- as.matrix(read.csv(file="./Sensitivity/FFGraphs-changeSize/FF-data_NodeEdgeDistanceMatrix.txt", sep=" ",header=T))

nckddata <- as.matrix(read.csv(file="./NetSimileGraphs/NS-Data-Less6/data_NodeEdgeDistanceMatrix.txt", sep=",",header=T))
nckddata <- as.matrix(read.csv(file="./NetSimileGraphs/NS-Data-Less6/data_NSdists.txt", sep=",",header=T))

nckddata <- as.matrix(read.csv(file="./NetSimileGraphs/NS-Data-LBD/NS-Data-less/Less-FWWC/NS-Data_NSdists.txt", sep=",",header=T))
nckddata <- as.matrix(read.csv(file="./NetSimileGraphs/NS-Data-LBD/NS-Data-less/Less-FWWC/NS-Data_LBDdists.txt", sep=",",header=T))
nckddata <- as.matrix(read.csv(file="./NetSimileGraphs/NS-Data-LBD/NS-Data-less/Less-FWWC/NS-Data_NodeEdgeDistanceMatrix.txt", sep=",",header=T))



meandistanc <- as.dist(nckddata)
hc <- hclust(meandistanc) 
dendro <- as.dendrogram(hc)
plot(dendro)
d <- dendrapply(dendro,  labelCol)
plot(d)
classificationnckd <-cutree(hc,k=numclasses)

y <- c(9,3,24,26,32,7,28,30,20,19,17,18,4,5,1,2,14,16,13,15,10,11,8,12,21,6,29,23,25,22,27,31)
rx <- reorder(d, y, agglo.FUN=mean)
order.dendrogram(rx) 
plot(rx,ylab="Distance",dLeaf=0.05)

par(cex=0.8,lty=1,lwd=2,font=1,font.axis=2)
plot(dendro)
d <- dendrapply(dendro,  labelCol)
plot(d,ylab="Distance",dLeaf=0.05)
par(lty=3,lwd=2)
rect.hclust(hc, k=numclasses, border="blue")



classificationNode <- cutree(hc1,k=numclasses)
classificationEdge <- cutree(hc2,k=numclasses)

names(classificationActual) <- names(classificationnckd)
compare(classificationActual, classificationnckd, method=c("nmi"))

classificationActual <- c(4,4,3,4,3,1,1,1,2,4,5,4,3,3,2,5,3,5,5,1,1,2,5,2,2)
names(classificationActual) <- names(classificationNode)
compare(classificationActual, classificationNode, method=c("nmi"))

classificationActual <- c(3,3,3,3,3,4,4,4,4,4,5,5,5,5,5,2,2,2,2,2,1,1,1,1,1)
names(classificationActual) <- names(classificationEdge)
compare(classificationActual, classificationEdge, method=c("nmi"))
#NetSimile Dendro
nsdata <- as.matrix(read.csv(file="./NetSimileGraphs/NS-Data_NSdists.txt", sep=",",header=T))
nsdata <- as.matrix(read.csv(file="./MetaboliteGraphs/BarabasiGraphs/AEData/AEGraphs_NSdists.txt", sep=",",header=T))
nsdata <- as.matrix(read.csv(file="./NetSimileGraphs/NS-Data-LBD/NS-Data_NSdists.txt", sep=",",header=T))

distanc0 <- as.dist(nsdata)
hc0 <- hclust(distanc0) 
dendro0 <- as.dendrogram(hc0)
plot(dendro0)
d0 <- dendrapply(dendro0,  labelCol)
plot(d0)
classificationnetsimile <-cutree(hc0,k=numclasses)
compare(classificationActual, classificationnetsimile, method=c("nmi"))

#LBD Dendro
lbddata <- as.matrix(read.csv(file="./NetSimileGraphs/NS-Data-LBD/NS-Data_LBDdists.txt", sep=",",header=T))
lbddata <- as.matrix(read.csv(file="./NetSimileGraphs/NS-Data-LBD/NS-Data_NSdists.txt", sep=",",header=T))
lbddata <- as.matrix(read.csv(file="./NetSimileGraphs/NS-Data-LBD/NS-Data_NodeDistanceMatrix.txt", sep=",",header=T))

distanc0 <- as.dist(lbddata)
hc0 <- hclust(distanc0) 
dendro0 <- as.dendrogram(hc0)
plot(dendro0)
d0 <- dendrapply(dendro0,  labelCol)
plot(d0)
rect.hclust(hc, k=numclasses, border="blue")
classificationnlbd <-cutree(hc0,k=numclasses)
compare(classificationActual, classificationnlbd, method=c("nmi"))


#NodeDistribution Dendro
nodedatafile <- "./NetSimileGraphs/NS-Data_NodeDistn.txt"
nodedatafile <- "./MetaboliteGraphs/BarabasiGraphs/AEData/newNCKD/AEGraphs_NewNodeDistnProbabilities.txt"

maxcols <- max(count.fields(nodedatafile, sep = ','))
data <- read.table(nodedatafile, header = FALSE, sep = ",", col.names = paste0("K",seq_len(maxcols)), fill = TRUE)
dmat <- as.matrix(data[2:ncol(data)])
rownames(dmat) <- data[,1]
pseudocount <- 0.0000001
dmat[is.na(dmat)] <- pseudocount
data.dist=dist.JSD(dmat)
write.table(as.matrix(data.dist),"./MetaboliteGraphs/BarabasiGraphs/AEData/newNCKD/AEGraphs_New_NodeDistanceMatrix.txt")

nodedist <- read.table(file="./NetSimileGraphs/NS-Data-LBD/NS-Data_NodeDistanceMatrix.txt",sep=",",header=F)
nodedist <- read.table(file="./NetSimileGraphs/NS-Data_NodeDistanceMatrix.txt",sep=",",header=T)

nodedist <- read.table(file="./Sensitivity/FFGraphs-changeSize/FF-data_NodeDistanceMatrix.txt",sep=" ",header=T)
nodedist <- as.matrix(read.csv(file="./NetSimileGraphs/NS-Data-Less6/data_NodeDistnProbabilities_JSDistanceMatrix.txt", sep=",",header=T))


data.dist <- as.dist(nodedist)
hc1 <- hclust(data.dist)   #node distn
classificationnode <-cutree(hc1,k=numclasses)
dendro1 <- as.dendrogram(hc1)
plot(dendro1)
d1 <- dendrapply(dendro1,  labelCol)
plot(d1) #, main = "Node Distn ",xlab="")

#EdgeDistribution
data <- read.table(file="./NetSimileGraphs/NS-Data_EdgeProbabilities.txt",sep=",",header=T)
data <- read.table(file="./MetaboliteGraphs/BarabasiGraphs/AEData/newNCKD/AEGraphs_NewIntraEdgeDistnProbabilities.txt",sep=",",header=F)

datamatrix <- as.matrix(data[,2:ncol(data)])
rownames(datamatrix) <- data[,1]
distanc <- dist.JSD(datamatrix)
write.table(as.matrix(distanc),"./MetaboliteGraphs/BarabasiGraphs/AEData/newNCKD/AEGraphs_New_IntraEdgeDistanceMatrix.txtEdgeDistanceMatrix.txt")

edgedist <- read.table(file="./NetSimileGraphs/NS-Data_EdgeDistanceMatrix.txt",sep=",",header=T)
edgedist <- read.table(file="./Sensitivity/FFGraphs-changeSize/FF-data_EdgeDistanceMatrix.txt",sep=" ",header=T)
edgedist <- as.matrix(read.csv(file="./NetSimileGraphs/NS-Data-Less6/data_EdgeProbabilities_JSDistanceMatrix.txt", sep=",",header=T))

distanc <- as.dist(edgedist)
hc2 <- hclust(distanc)   #edge distn
classificationedge <-cutree(hc2,k=numclasses)
dendro2 <- as.dendrogram(hc2)
#plot(dendro1)
d2 <- dendrapply(dendro2,  labelCol)
plot(d2) #, main = "Edge Distn ",xlab="")
            

data.dist <- as.matrix(read.csv(file="./NetSimileGraphs/NS-Data-LBD/NS-Data-less/Less-FWWC/data_NodeDistnProbabilities_JSDistanceMatrix.txt", sep=",",header=T))
distanc <- as.matrix(read.csv(file="./NetSimileGraphs/NS-Data-LBD/NS-Data-less/Less-FWWC/data_EdgeProbabilities_JSDistanceMatrix.txt", sep=",",header=T))


meandistanc <- (data.dist + distanc)/2
write.table(as.matrix(meandistanc),"./NetSimileGraphs/NS-Data-LBD/NS-Data-less/Less-FWWC/NS-Data_NodeEdgeDistanceMatrix.txt")

nodedist <- read.table(file="./TestClassification/TwoClassData/TC2-Data_NodeDistnProbabilities_JSDistanceMatrix.txt",sep=",",header=T)
edgedist <- read.table(file="./TestClassification/TwoClassData/TC2-Data_EdgeProbabilities_JSDistanceMatrix.txt",sep=",",header=T)
meandistanc <- (nodedist + edgedist)/2
write.table(as.matrix(meandistanc),"./TestClassification/TwoClassData/TC2-Data_NodeEdgeProbabilitie_JSDistanceMatrix.txt")


labelCol <- function(x) {
  if (is.leaf(x)) {
    ## fetch label
    labelstring <- attr(x, "label") 
    label <- substr(labelstring,1,7)
    ## set label color to red for A and B, to blue otherwise
    #attr(x, "nodePar") <- list(lab.col=ifelse(label %in% c("A", "B"), "red", "blue"))
    if (grepl(labelnames[1],label,)){
      attr(x, "nodePar") <- list(lab.col=colors[1])
    } else if (grepl(labelnames[2],label)) {
      attr(x, "nodePar") <- list(lab.col=colors[2])
    }else if (grepl(labelnames[3],label)){
      attr(x, "nodePar") <- list(lab.col=colors[3])
    } else if (grepl(labelnames[4],label)) {
      attr(x, "nodePar") <- list(lab.col=colors[4])
    } else if (grepl(labelnames[5],label)) {
      attr(x, "nodePar") <- list(lab.col=colors[5])
    } else if (grepl(labelnames[6],label)) {
      attr(x, "nodePar") <- list(lab.col=colors[6])
    } else if (grepl(labelnames[7],label)) {
      attr(x, "nodePar") <- list(lab.col=colors[7])
    } else if (grepl(labelnames[8],label)) {
      attr(x, "nodePar") <- list(lab.col=colors[8])
    } else if (grepl(labelnames[9],label)) {
      attr(x, "nodePar") <- list(lab.col=colors[9])
    } else {
      attr(x, "nodePar") <- list(lab.col="black")
    }
    
    #attr(x, "nodePar") <- labelcolor
  }
  return(x)
}
