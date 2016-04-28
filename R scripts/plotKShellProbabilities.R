library('igraph')
setwd("~/graph-algos/NetworkComparison/data")

paste0 <- function(..., collapse = NULL) {
  paste(..., sep = "", collapse = collapse)
}
datafile <- '../Python-Programs/TestData/SCKCom/SC.nodedistribution'

datafile <- "../Python-Programs/TestData/data/weighted_NodeDistn.txt"
datafile <- "./MetaboliteGraphs/MakeGraphsSimple/MGSimpleData_mCorePD.txt"
datafile <- "./DeletedEdgesGraphs/SecondExperiment/Run1/CA_NodeDistnProbabilities.txt"
datafile <- "./LineGraphsData/BarabasiGraphsData_LineGraph_NodeDistn.txt"
datafile <- "./LineGraphsData/MGSimpleData_LineGraph_NodeDistn.txt"
datafile <- "./WeightedKCore/CA_Weighted_NodeDistn.txt"
datafile <- "./NetSimileGraphs/NS-Data-M1-M6-NodeDistn"
datafile <- "./SocialCapitalGraphs/SCDistn_CA.txt"
datafile <- "./NetSimileGraphs/NodeDistributionPlots/CA-NodeDistn"

count.fields(datafile, sep = ",")
maxcols <- max(count.fields(datafile, sep = ','))   #maxcols <- 239
#data <- read.csv(datafile,header=F)
data <- read.table(datafile, header = FALSE, sep = ",", col.names = paste0("V",seq_len(maxcols)), fill = TRUE)

#data <- read.table(datafile, header = FALSE, sep = ",", row.names=1,col.names = paste0("V",seq_len(maxcols)), fill = TRUE)
outputfile <- '../Python-Programs/TestData/SCKCom/tobesent/coauthor-nodedistn.png'
png(filename=outputfile,width=360,height=350,,res=72)
mat <- as.matrix(data[])
tmat <- t(mat)
par(mar=c(5, 6, 4, 2) + 0.1)
matplot(1:maxcols, tmat, type="l", xlab="Coreness (K)", ylab='Normalized Frequency of \n nodes with coreness K',axes=F, cols=c('red','blue','black'),lwd=2)
axis(side=1, at=seq(1,maxcols, by=10),font=2)
axis(side=2, at=seq(0, 0.3, by=0.01),font=2)
legend("topright",legend=c('DBLP','USPTO','AMiner'), lty=1,lwd=2,col=c("red","blue","darkgreen"), ncol=3,bty="n",cex=0.8,text.col=c("red","blue","darkgreen")) #, inset=0.01)
box()
dev.off()
matplot(1:maxcols, tmat, type="l",xlab="Shell Number", ylab="P", axes=FALSE)
axis(side=1, at=seq(1,10, by=1))
axis(side=2, at=seq(0, 1.2, by=0.1))
box()
dev.off()

matplot(1:maxcols, tmat, type="l", xlab="Shell Number", ylab=expression('P(S'[K]*'(G))'))
matplot(1:maxcols, tmat[, 2:3], type="l", xlab="K-Shell Number", ylab=expression('P(S'[K]*'(G))'))

matplot(1:maxcols, tmat, type="l", xlab="K-Shell Number", ylab="P(num nodes in Kth Shell)", main="Autonomous Systems K-Shell Probability Distribution")
matplot(1:maxcols, tmat, type="l", xlab="K-Shell Number", ylab="P(num nodes in Kth Shell)", main="Forest Fire Networks K-Shell Probability Distribution")
matplot(1:maxcols, tmat, type="l",  xlim= c(1,5),xlab="K-Shell Number", ylab="P(num nodes in Kth Shell)", main="Erdos Reyni Networks K-Shell Probability Distribution")
matplot(1:maxcols, tmat, type="l",  xlab="K-Shell Number", ylab="P(num nodes in Kth Shell)", main="Watts Strogatz Networks K-Shell Probability Distribution")
matplot(1:maxcols, tmat, type="l",  xlab="K-Shell Number", ylab="P(num nodes in Kth Shell)", main="Metabolic Networks K-Shell Probability Distribution")
matplot(1:maxcols, tmat, type="l",  xlab="K-Shell Number", ylab="P(num nodes in Kth Shell)", main="CoAuthorship Networks K-Shell Probability Distribution")
matplot(1:maxcols, tmat, type="l",  xlab="K-Shell Number", ylab="P(num nodes in Kth Shell)", main="DBLP-Yearwise K-Shell Probability Distribution")
matplot(1:maxcols, tmat, type="l",  xlim = c(1,5),xlab="K-Shell Number", ylab="P(num nodes in Kth Shell)", main="Barabassi Albert K-Shell Probability Distribution")

matplot(1:maxcols, tmat, type="l",xlab="Degree", ylab="P(degree)", main="Erdos Reyni Networks Degree Distribution")
matplot(1:maxcols, tmat, type="l",xlab="Degree", ylab="P(degree)", main="Erdos Reyni Networks Degree Distribution", axes=FALSE)
axis(side=1, at=seq(1,10, by=1))

lab <- c("G(n,2n)","G(n,3n)","G(n,4n)", "G(n,5n)","G(n,6n)") 
text(locator(5), lab, adj=0) 
#legend("bottomright", inset=.05, legend=c("b", "c"), pch=1, col=c(2,4), horiz=TRUE)

matplot(1:maxcols, tmat, type="l",  xlab="Degree", ylab="P(degree)", main="Erdos Reyni Networks Degree Distribution",axes=FALSE)
axis(side=1, at=seq(0, 30, by=2))
axis(side=2, at=seq(0, 0.3, by=0.1))
lab <- c("G(n,2n)","G(n,3n)","G(n,4n)", "G(n,5n)","G(n,6n)") 
text(locator(5), lab, adj=0) 
box()

matplot(1:maxcols, tmat[, 7:21], type="l",xlab="KShellNumber", ylab="Probability numnodes", main="ERModel-KShells Prob Distn")
matplot(1:maxcols, tmat[, 22:31], type="l", xlab="KShellNumber", ylab="Probability numnodes", main="BAModel-KShells Prob Distn")
matplot(1:maxcols, tmat[, 32:41], type="l", xlab="KShellNumber", ylab="Probability numnodes", main="FFModel-KShells Prob Distn")
matplot(1:maxcols, tmat[, 42:51], type="l", xlab="KShellNumber", ylab="Probability numnodes", main="WSModel-KShells Prob Distn")
matplot(1:maxcols, tmat[, 52:56], ylim=c(0,0.3),type="l", xlab="KShellNumber", ylab="Probability numnodes", main="CA-KShells Prob Distn")
matplot(1:maxcols, tmat[, 57:65], type="l", xlab="KShellNumber", ylab="Probability numnodes", main="AS-KShells Prob Distn")
matplot(1:maxcols, tmat[, 65:108], type="l", xlab="KShellNumber", ylab="Probability numnodes", main="MetaboliteModel-KShells Prob Distn")



categories <- c("Archae", "Bacteria", "Eukaryote")
colors <- c("green", "blue", "red")
markers <- 1:3
# Plot
matplot(1:maxcols, tmat[, 1:6], type="l", col=c("green"), lty=1, pch=markers, bty="n", las=1, xlab="KShellNumber", ylab="Probability numnodes", main="MetaboliteNetworks-KShells Prob Distn")
par(new=TRUE)
matplot(1:maxcols, tmat[, 7:30], type="l", col=c("blue"), lty=1, pch=markers, bty="n", las=1, xlab="KShellNumber", ylab="Probability numnodes", main="MetaboliteNetworks-KShells Prob Distn")
par(new=TRUE)
matplot(1:maxcols, tmat[, 31:35], type="l", col=c("red"), lty=1, pch=markers, bty="n", las=1, xlab="KShellNumber", ylab="Probability numnodes", main="MetaboliteNetworks-KShells Prob Distn")
legend("topright", col=colors, categories, bg="white", lwd=1, pch=markers)
