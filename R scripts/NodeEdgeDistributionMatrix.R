setwd("/home/stratuser/graph-algos/NetworkComparison/data")
#Node Distribution Matrix
dataname <- "./RewiredGraphs/SecondExperiment/SameGraphs/Run3/DY1"
dataname <- "./RewiredGraphs/SecondExperiment/SameGraphs/Run3/DY3"
dataname <- "./RewiredGraphs/SecondExperiment/SameGraphs/Run3/CA4"

dataname <- "./RewiredGraphs/SecondExperiment/Run3/ER"

nodedatafile <- paste0(dataname,"_NodeDistnProbabilities.txt")
maxcols <- max(count.fields(nodedatafile, sep = ','))
data <- read.table(nodedatafile, header = FALSE, sep = ",", col.names = paste0("K",seq_len(maxcols)), fill = TRUE)
dmat <- as.matrix(data[2:ncol(data)])
rownames(dmat) <- data[,1]
pseudocount <- 0.0000001
dmat[is.na(dmat)] <- pseudocount
data.dist=dist.JSD(dmat)
outputfile <- paste0(dataname,"_NodeDistanceMatrix.txt")
write.table(as.matrix(data.dist),outputfile)

edgedatafile <- paste0(dataname,"_EdgeProbabilities.txt")
data <- read.table(file=edgedatafile,sep=",",header=F)
datamatrix <- as.matrix(data[,2:ncol(data)])
rownames(datamatrix) <- data[,1]
distanc <- dist.JSD(datamatrix)
outputfile <- paste0(dataname,"_EdgeDistanceMatrix.txt")
write.table(as.matrix(distanc),outputfile)

mean <- (distanc + data.dist)/2
mean <- as.matrix(mean)
mean <- as.dist(mean[order(rownames(mean)),order(colnames(mean))])
outputfile <- paste0(dataname,"_NodeEdgeDistanceMatrix.txt")
write.table(as.matrix(mean),outputfile)

