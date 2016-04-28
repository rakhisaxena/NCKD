library('emdist')

dist.EMD <- function(inMatrix, pseudocount=0.000001, ...) {
  matrixColSize <- length(colnames(inMatrix))
  matrixRowSize <- nrow(inMatrix)
  resultsMatrix <- matrix(0, matrixRowSize, matrixRowSize)
  
  for(i in 1:matrixRowSize) {
    for(j in i:matrixRowSize) { 
      resultsMatrix[i,j]=emd(as.matrix(inMatrix[i,]), as.matrix(inMatrix[j,]))
      resultsMatrix[j,i]=emd(as.matrix(inMatrix[i,]), as.matrix(inMatrix[j,]))
    }
  }
  rownames <- rownames(inMatrix)
  rownames -> rownames(resultsMatrix) -> colnames(resultsMatrix)
  as.dist(resultsMatrix)->resultsMatrix
  attr(resultsMatrix, "method") <- "dist"
  return(resultsMatrix) 
}

data <- read.table(file="../Python-Programs/TestData_EdgeProbabilities.txt",sep=",",header=F)

datamatrix <- as.matrix(data[,2:ncol(data)])
rownames(datamatrix) <- data[,1]
distanc <- dist.EMD(datamatrix)
hc1 <- hclust(distanc)   # LBD Distribution
dendro1 <- as.dendrogram(hc1)
d1 <- dendrapply(dendro1,  labelCol)
plot(d1) #, main = "Edge Distn ",xlab="")
