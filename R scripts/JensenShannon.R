paste0 <- function(..., collapse = NULL) {
  paste(..., sep = "", collapse = collapse)
}


dist.JSD <- function(inMatrix, pseudocount=0.000001, ...) {
  KLD <- function(x,y) sum(x *log(x/y))
  JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
  matrixColSize <- length(colnames(inMatrix))
  matrixRowSize <- nrow(inMatrix)
  resultsMatrix <- matrix(0, matrixRowSize, matrixRowSize)
  
  inMatrix[is.na(inMatrix)] <- pseudocount
  inMatrix[inMatrix==0] <- pseudocount
  
  for(i in 1:matrixRowSize) {
    for(j in i:matrixRowSize) { 
      resultsMatrix[i,j]=JSD(as.vector(inMatrix[i,]), as.vector(inMatrix[j,]))
      resultsMatrix[j,i]=JSD(as.vector(inMatrix[i,]), as.vector(inMatrix[j,]))
    }
  }
  rownames <- rownames(inMatrix)
  rownames -> rownames(resultsMatrix) -> colnames(resultsMatrix)
  as.dist(resultsMatrix)->resultsMatrix
  attr(resultsMatrix, "method") <- "dist"
  return(resultsMatrix) 
}


setwd("~/graph-algos/NetworkComparison/data/DeletedEdgesGraphs")
datafile <- "./MetabolicNetworks-KCoreSignature"
datafile <- "./TestClassification/data_KCoreSignatures.txt"
datafile <- "./EdgeMatrixTest"
datafile <- "./KShellProbabilityDistributionTest.txt"
datafile <- "./NetSimileVsKCore-Outputs/EdgesDeleted_KCoreFeatures.txt"
datafile <- "./NetSimileVsKCore-Outputs/EdgesDeleted-1_KCoreFeatures.txt"
datafile <- "./DBLP-C-Rewired-KCore-ProbDistn.txt"
datafile <- "./ER_KCoreSignatures.txt"
datafile <- "./ER_KCoreSignatures.txt"
datafile <- "./MetaboliteGraphs/BarabasiGraphs/BarabasiGraphsData_NodeDistn.txt"
datafile <- "./MetaboliteGraphs/Barabasi-NodeDistnClassificationData"
datafile <- "./MetaboliteGraphs/Barabasi-NodeDistnClassificationData-Giant"
datafile <- "./MetaboliteGraphs/MakeGraphs/MakeGraphsData_NodeDistn.txt"
datafile <- "./MetaboliteGraphs/BarabasiGraphs/AEData/AEGraphs_NodeDistn.txt"
datafile <- "./MetaboliteGraphs/MakeGraphsSimple/MGSimpleData_giant_mCorePD.txt"
datafile <- "./MetaboliteGraphs/MakeGraphsSimple/MGSimpleData_dCorePD.txt"
datafile <- "./TestClassification/TCData_NodeDistn.txt"
datafile <- "./NetSimileGraphs/NS-Data_mCorePD.txt"
datafile <-"./NetSimileGraphs/NS-Data/data_Weighted_NodeDistn.txt"


#datafile <- "./TestGraphEigenValues.txt"
maxcols <- max(count.fields(datafile, sep = ','))
data <- read.table(datafile, header = FALSE, sep = ",", col.names = paste0("K",seq_len(maxcols)), fill = TRUE)
dmat <- as.matrix(data[2:ncol(data)])
rownames(dmat) <- data[,1]
pseudocount <- 0.0000001
dmat[is.na(dmat)] <- pseudocount
#dmat[dmat==0] <- pseudocount
data.dist=dist.JSD(dmat)
#data.dist
#write.table(as.matrix(data.dist),"KShellTest-Eigen-DistanceMatrix.csv")
#write.table(as.matrix(data.dist),"./ER-KCore-DistanceMatrix.txt")
write.table(as.matrix(data.dist),"./NetSimileGraphs/NS-Data/data_Weighted_DistanceMatrix.txt")
library(ape)
hc = hclust(data.dist)
dendro <- as.dendrogram(hc)
#plot(dendro)
d <- dendrapply(dendro,  labelCol)
plot(d)

windows()
par(cex=0.7)
plot(hcd)
plot(as.phylo(hc), cex = 0.9, label.offset = 1)
png("ERGraphs-KCore-Dendrogram.png",width=1600,height=800)
#par(cex=0.7,font=3)
plot(hcd, main= "KCoreAlgorithm: Comparison of ER Networks", xlab="Datasets", axes=T)
dev.off()


labelCol <- function(x) {
  if (is.leaf(x)) {
    ## fetch label
    label <- attr(x, "label") 
    ## set label color to red for A and B, to blue otherwise
    #attr(x, "nodePar") <- list(lab.col=ifelse(label %in% c("A", "B"), "red", "blue"))
    if (grepl("A-",label)){
      attr(x, "nodePar") <- list(lab.col="red")
    } else if (grepl("E-",label)){
      attr(x, "nodePar") <- list(lab.col="green")
    } else {
      attr(x, "nodePar") <- list(lab.col="blue")
    }
    #attr(x, "nodePar") <- labelcolor
  }
  return(x)
}

d <- dendrapply(hcd,  labelCol)
png("MetaboliteDendrogram.png",width=1200)
plot(d, main= "Comparison of Metabolic Networks", xlab="Organisms",)
dev.off()


kc <- kmeans(mat, 3)
kc



dist.Cosine <- function(inMatrix, pseudocount=0.000001, ...) {
  mycosine <- function(x,y){ sum(x*y) / (sqrt(sum(x*x)) * sqrt(sum(y*y))) }
  matrixColSize <- length(colnames(inMatrix))
  matrixRowSize <- nrow(inMatrix)
  resultsMatrix <- matrix(0, matrixRowSize, matrixRowSize)
  
  inMatrix[is.na(inMatrix)] <- pseudocount
  inMatrix[inMatrix==0] <- pseudocount
  
  for(i in 1:matrixRowSize) {
    for(j in i:matrixRowSize) { 
      resultsMatrix[i,j]=mycosine(as.vector(inMatrix[i,]), as.vector(inMatrix[j,]))
      resultsMatrix[j,i]=mycosine(as.vector(inMatrix[i,]), as.vector(inMatrix[j,]))
    }
  }
  rownames <- rownames(inMatrix)
  rownames -> rownames(resultsMatrix) -> colnames(resultsMatrix)
  as.dist(resultsMatrix)->resultsMatrix
  attr(resultsMatrix, "method") <- "dist"
  return(resultsMatrix) 
}

