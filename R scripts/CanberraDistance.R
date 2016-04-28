setwd("~/graph-algos/NetworkComparison/data/RewiredGraphs/")
datafile <- "./TestGraphEigenValues.txt"
datafile <- "./NetSimileVsKCore-Outputs/EdgesDeleted_NetSimileFeatures.txt"
datafile <- "./NetSimileVsKCore-Outputs/EdgesDeleted-1_NetSimileFeatures.txt"
datafile <- "./NetSimileVsKCore-Outputs/NSvsKC-Sorted-NetSimileFeatures.txt"
datafile <- "./DBLP-C_NetSimileFeatures.txt"
maxcols <- max(count.fields(datafile, sep = ','))
data <- read.table(datafile, header = FALSE, sep = ",", col.names = paste0("K",seq_len(maxcols)), fill = TRUE)
dmat <- as.matrix(data[2:maxcols])
rownames(dmat) <- data[,1]
simmat <- 1- dmat
pseudocount <- 0.0000001
dmat[is.na(dmat)] <- pseudocount
#dmat[dmat==0] <- pseudocount
ncol(dmat)
data.dist= dist(dmat, method = "canberra", diag=T)
data.dist = data.dist/ncol(dmat)
#data.dist
#write.table(as.matrix(data.dist),"KShellTest-Eigen-DistanceMatrix.csv")
write.table(as.matrix(data.dist),"./DBLP-C-Rewired-NetSimile-DistanceMatrix.txt")
library(ape)
hc = hclust(data.dist)
hcd = as.dendrogram(hc)
plot(hcd)

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
