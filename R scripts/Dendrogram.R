library('igraph')
library('clue')
setwd("~/graph-algos/NetworkComparison/data")
#NS
classificationActual <- c(1,1,1,1,1,2,2,2,2,3,3,3,3,3,6,6,6,6,4,4,4,4,5,5,5,5,5,5,7,7,7,9,9,9,9,8,8,8)
#nNCKD
classificationActual <- c(9,9,9,9,9,2,2,2,2,3,3,3,3,3,1,1,1,1,4,4,4,4,8,8,8,8,8,8,5,5,5,6,6,6,6,7,7,7)
names(classificationActual) <- names(classificationNode)
numclasses= 2

#data <- as.matrix(read.csv(file=file.choose(), sep=",",header=T))
#data <- as.matrix(read.table(file=file.choose(), sep=",",header=T))

#0.8564313 0.8896605  0.8787436
data <- as.matrix(read.csv(file="./NetSimileGraphs/NS-Data-LBD/NS-Data-less/NS-Data_LBDdists.txt", sep=",",header=T))
data <- as.matrix(read.csv(file="./NetSimileGraphs/NS-Data-LBD/NS-Data-less/NS-Data_NSdists.txt", sep=",",header=T))
data <- as.matrix(read.csv(file="./NetSimileGraphs/NS-Data-LBD/NS-Data-less/NS-Data_JSdists.txt", sep=",",header=T))
data <- as.matrix(read.csv(file="./NetSimileGraphs/NS-Data-LBD/PreviousOutputs/NS-Data_NodeDistanceMatrix.txt", sep=",",header=T))

data <- as.matrix(read.csv(file="./MetaboliteGraphs/BarabasiGraphs/AEData/NCKD/AEGraphs_NodeEdgeDistanceMatrix.txt", sep=",",header=T))
data <- as.matrix(read.csv(file="./TestClassification/TCData_NodeEdgeMaxDistanceMatrix.txt", sep=",",header=T))

numclasses=2

distanc <- as.dist(data)
hc <- hclust(distanc) 
dendro <- as.dendrogram(hc)
par(cex=0.8,lty=1,lwd=2,font=1,font.axis=2)
plot(dendro)
d <- dendrapply(dendro,  labelCol)
plot(d,ylab="Distance",dLeaf=0.05)
par(lty=3,lwd=2)
rect.hclust(hc, k=numclasses, border="blue")

classificationNode <-cutree(hc,k=numclasses)
classificationNode<-classificationNode[names(classificationActual)]
nmi = compare(classificationActual, classificationNode, method=c("nmi"))
nmi
legend("topright", "NMI=0.2746372", bty="n") 
#581 X 317 pixels

data <- read.table(file="./MetaboliteGraphs/BarabasiGraphs/AEData/AEGraphs_NSdists.txt",sep=",",header=T)
data <- read.table(file="./MetaboliteGraphs/BarabasiGraphs/AEData/AEGraphs_NodeDistanceMatrix.txt",sep=",",header=T)

data <- as.matrix(read.csv(file="./KShellTest-KCore-DistanceMatrix.csv", sep="",header=T))
data <- as.matrix(read.csv(file="./KShellTest-Eigen-DistanceMatrix.csv", sep="",header=T))
data <- as.matrix(read.csv(file="./KShellTest-EIG-DistanceMatrix.csv", sep="",header=T)


class(d)
order.dendrogram(d)
y <- c(14,15, 13,9,12,10,7,6,8,11,2,3,5,1,4)
y <- c(1,4 ,2,3,5, 10,7,6  ,8  ,9 ,13, 11, 12, 14, 15  )
y <- c(10,7,9,8,11,4,5,3,6,1,2)

y <- c(9,3,24,26,32,7,28,30,20,19,17,18,4,5,1,2,14,16,13,15,10,11,8,12,21,6,29,23,25,22,27,31)

distanc <- as.dist(data)
hc <- hclust(distanc) 

dendro <- as.dendrogram(hc)
par(cex=0.8,lty=1,lwd=2,font=2,font.axis=2)
plot(dendro)
d <- dendrapply(dendro,  labelCol)
plot(d,ylab="Distance")
par(lty=2,lwd=2)
rect.hclust(hc, k=numclasses, border="blue")
rx <- reorder(d, y, agglo.FUN=mean)
order.dendrogram(rx) 
plot(rx,,ylab="Distance")





#NetSimile Dendrogram
data <- as.matrix(read.csv(file="./MetaboliteGraphs/BarabasiGraphs/AEData/AEGraphs_NSdists.txt", sep=",",header=T))

#Node Distribution Dendrogram
datafile <- "./NetSimileGraphs/NS-Data_NodeDistn.txt"
maxcols <- max(count.fields(datafile, sep = ','))
data <- read.table(datafile, header = FALSE, sep = ",", col.names = paste0("K",seq_len(maxcols)), fill = TRUE)
dmat <- as.matrix(data[2:ncol(data)])
rownames(dmat) <- data[,1]
pseudocount <- 0.0000001
dmat[is.na(dmat)] <- pseudocount
data.dist=dist.JSD(dmat)
write.table(as.matrix(data.dist),"./DeletedEdgesGraphs/SecondExperiment/Run1/AS_NodeDistanceMatrix.txt")

datamatrix <- as.matrix(data[,2:ncol(data)])
rownames(datamatrix) <- data[,1]
distanc <- dist.JSD(datamatrix)
write.table(as.matrix(data.dist),"./DeletedEdgesGraphs/SecondExperiment/Run1/AS_EdgeDistanceMatrix.txt")

write.table(as.matrix(distanc),"./TestClassification/TCData_EdgeDistanceMatrix.txt")
data <- as.matrix(read.csv(file="./TestClassification/TCData_EdgeDistanceMatrix.txt", sep=",",header=T))
distanc <- as.dist(data)

hc1 <- hclust(distanc)   #edge distn
dendro1 <- as.dendrogram(hc1)
plot(dendro1)
d1 <- dendrapply(dendro1,  labelCol)
plot(d1) #, main = "Edge Distn ",xlab="")


data1 <- as.matrix(read.csv(file="./Sensitivity/WSGraphs-changeSize/WS-data_NodeDistanceMatrix.txt", sep=" ",header=T))
dist1 <- as.dist(data1)


mean <- (distanc + data.dist)/2
write.table(as.matrix(data.dist),"./Sensitivity/WSGraphs-changeSize/WS-data_NodeEdgeDistanceMatrix.txt")
write.table(as.matrix(mean),"./DeletedEdgesGraphs/SecondExperiment/Run1/AS_NodeEdgeDistanceMatrix.txt")

data <- as.matrix(read.csv(file="./ScalabilityGraphs/SGData_NodeEdgeDistanceMatrix.txt", sep=",",header=T))

data <- as.matrix(read.csv(file="./Sensitivity/WS_changeD_NodeEdgeDistanceMatrix.txt", sep="",header=T))

mean <- as.dist(data)
 
hc3 <- hclust(mean) # node+edge
classificationnodeedge <- cutree(hc3,k=2)
dendro3<- as.dendrogram(hc3)
plot(dendro3)
d3 <- dendrapply(dendro3,  labelCol)
plot(d3)
classificationActual <- c(rep(1,5), rep(2,6))
names(classificationActual) <- names(classificationnetsimile)
compare(classificationActual, classificationnodeedge, method=c("nmi"))

class(d3)
order.dendrogram(d3)
plot(d3)
y <- c(13,18,20,15,19,24,16,21,23,25,17,22,7,8,9,11,14,12,6,10,4,5,2,1,3)

y <- c(13,18,20,15,19,24,16,21,23,25,17,22,7,8,9,6,10,11,14,12,4,5,2,1,3)
rx <- reorder(d3, y, agglo.FUN=mean)
order.dendrogram(rx) 
plot(rx)







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



classification3<-cutree(hc,k=3)

classificationGiant <- classification1
classificationEdge <- classification2
classificationKCore <- classification3
classificationActual <- c(rep(1,6),rep(2,32),rep(3,5))
names(classificationActual) <- names(classificationKCore)
rev(classificationActual)
compare(classificationcombined, classificationActual,method=c("nmi"))

m1 <- as.numeric(datamatrix[1,])
m1 <- matrix(m1, nrow=6, ncol=6)
m2 <- as.numeric(datamatrix[2,])
m2 <- matrix(m2, nrow=6, ncol=6)
m3 <- as.numeric(datamatrix[7,])
m3 <- matrix(m3, nrow=6, ncol=6)
m4 <- as.numeric(datamatrix[8,])
m4 <- matrix(m4, nrow=6, ncol=6)
m5 <- as.numeric(datamatrix[41,])
m5 <- matrix(m5, nrow=6, ncol=6)
m6 <- as.numeric(datamatrix[42,])
m6 <- matrix(m6, nrow=6, ncol=6)
m7 <- as.numeric(datamatrix[43,])
m7 <- matrix(m7, nrow=6, ncol=6)

sum(dist(m1,m4, pairwise=T))

sum(dist(m1,m2,pairwise=T))
sum(dist(m1,m4,pairwise=T))
sum(dist(m3,m4,pairwise=T))
sum(dist(m5,m6,pairwise=T))
sum(dist(m6,m7,pairwise=T))
