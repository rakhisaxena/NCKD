library(igraph)
setwd("/home/stratuser/graph-algos/NetworkComparison/data")

#rewire edges
datafile =read.csv("./ER10K-00", sep = "")
g=graph.data.frame(datafile, directed = FALSE)
g <- erdos.renyi.game(n=1000, p.or.m=5000, type="gnm")
write.graph(g, "./DBLP-2014-0.2" ,format=c("edgelist"))
g1 <- rewire.edges(g, p=0.2)
write.graph(g1, "./ER10K-0.2" ,format=c("edgelist"))

#delete 5% edges
datafile =read.csv("./AllGraphs/ERGraphs/ER3-1K-1.txt", sep = "")
g=graph.data.frame(datafile, directed = FALSE)
g1 <- erdos.renyi.game(n=1000, p.or.m=5000, type="gnm")
write.graph(g1, "./AllGraphs/ERGraphs/EdgesDeleted-1/ER5-1K-d00.txt" ,format=c("edgelist"))
g <-g1
numedgestodelete <- 0.50* ecount(g)
for (i in 1:numedgestodelete) {
  e <- sample(E(g),1)
  g <- delete.edges(g,E(g)[e])
}
write.graph(g, "./AllGraphs/ERGraphs/EdgesDeleted-1/ER5-1K-d50.txt" ,format=c("edgelist"))
ecount(g)
plot(degree.distribution(g), xlab="node degree", cex=0.1)
lines(degree.distribution(g))

#Rewire a directory
inputgraphDirectory <- "~/graph-algos/NetworkComparison/data/RewiredGraphs/SecondExperiment/SameGraphs/Run3/DY1"
setwd(inputgraphDirectory)
inputfile <- "./DY-1"
datafile =read.csv(inputfile, sep = "")
g=graph.data.frame(datafile, directed = FALSE)
fname <- gsub("./|.txt|.dat|.csv" , "" , inputfile)
for (prob in seq(from=0.02, to=0.2, by=0.02)) {
  g1 <- rewire.edges(g, p=prob)
  outputfile <- paste(fname,"-",prob*10,sep="")
  write.graph(g1, outputfile ,format=c("edgelist"))
}


#Edge Delete a directory
inputgraphDirectory <- "~/graph-algos/NetworkComparison/data/DeletedEdgesGraphs/SecondExperiment/SameGraphs/Run3/DY3"
setwd(inputgraphDirectory)
inputfile <- "./DY-3"
datafile =read.csv(inputfile, sep = "")
g=graph.data.frame(datafile, directed = FALSE)
fname <- gsub("./|.txt|.dat|.csv" , "" , inputfile)
for (num  in seq(from=0.02, to=0.2, by=0.02)) {
  numedgestodelete <- num * ecount(g)
  e <- sample(E(g),numedgestodelete)
  g1 <- delete.edges(g,E(g)[e])
  outputfile <- paste(fname,"-",num*10,sep="")
  write.graph(g1, outputfile ,format=c("edgelist"))
}



ns <- c(0,4.91441590659022,4.40622741316813,4.21928126617076,4.32660708938112,3.17977136157065,8.57972942495426,7.52352773366857,8.08942630241684,2.18131956161667,12.6525331411439)
plot(ns)
lines(ns)
ns <- c(0, 0.154913725090265,0.247297744179949,0.349158718652402,0.425055293603848,0.458963104020677)
sim <- 1 - ns
plot(sim)
lines(sim)

sim


data <- read.csv("./RewiredGraphs/PlotData-KCore", sep="",header=TRUE)
dmat <- as.matrix(data)
smat <- 1 - dmat
matplot(1:ncol(smat), t(smat), type="l", xlab="Rewiring parameter",ylab="NetSimile Similarity Score", axes=F)
axis(1,at=1:6,labels=seq(0,1.0,0.2))
axis(2,at=1:6,labels=seq(0,1.0,0.2))
smat[1,]

data <- read.csv("./RewiredGraphs/SecondExperiment/Plotdata-NCKD-Averaged", sep="",header=TRUE)
data <- read.csv("./RewiredGraphs/SecondExperiment/Plotdata-NetSimile-Averaged", sep="",header=TRUE)
dmat <- as.matrix(data)
smat <- 1 - dmat
par(mar=c(4.0, 4.0, 1.5, 1.5))
#NCKD
plot(smat[1,],type="b",lty=1,lwd=2, xaxt="n",ylim=c(0.2,1.0),col="black", xlab="Percentage of edges rewired",ylab="NCKD SS")
#Netsimile
plot(smat[1,],type="b",lty=1,lwd=2, xaxt="n",ylim=c(0.2,1.0),col="black", xlab="Percentage of edges rewired",ylab="NetSimile SS")

axis(1,at=1:11,labels=seq(0,20,2))
lines(smat[2,],col="red",type="b",lty=2,lwd=2)
lines(smat[3,],col="mediumvioletred",type="b",lty=3,lwd=2)
lines(smat[4,],col="blue",type="b",lty=4,lwd=2)
legend("bottomleft",legend=c("A-1","CA-1","DC-1","E10K-1"), lty=c(1,2,3,4),lwd=2,pch=21,col=c("black","red","mediumvioletred", "blue"), ncol=4,bty="n",cex=1.0,text.col=c("black","red","mediumvioletred","blue"), inset=0.01)

axis(1,at=1:11,labels=seq(0,0.2,0.02))
lines(smat[2,],type="b",lty=2,lwd=2)
lines(smat[3,],type="b",lty=3,lwd=2)
lines(smat[4,],type="b",lty=4,lwd=2)
legend("bottomleft",legend=c("A-1","CA-1","DC-1","E10K-1"), lty=c(1,2,3,4),lwd=2,pch=21,col=c("black","red","green", "blue"), ncol=4,bty="n",cex=0.8,text.col=c("black","red","green","blue"), inset=0.01)

data <- read.csv("./DeletedEdgesGraphs/SecondExperiment/Plotdata-KCore-Averaged", sep="",header=TRUE)
data <- read.csv("./DeletedEdgesGraphs/SecondExperiment/Plotdata-NetSimile-Averaged", sep="",header=TRUE)



data <- read.csv("./DeletedEdgesGraphs/SecondExperiment/Plotdata-SameGraphs-KCore-Averaged", sep="",header=TRUE)
data <- read.csv("./DeletedEdgesGraphs/SecondExperiment/Plotdata-SameGraphs-NetSimile-Averaged", sep="",header=TRUE)

dmat <- as.matrix(data)
smat <- 1 - dmat
par(mar=c(4.0, 4.0, 1.5, 1.5))
#NetSimile
plot(smat[1,],type="b",lty=1,lwd=2, xaxt="n",ylim=c(0.6,1.0),col="black", xlab="Percentage of edges deleted",ylab="NetSimile SS")
#NCKD
plot(smat[1,],type="b",lty=1,lwd=2, xaxt="n",ylim=c(0.2,1.0),col="black", xlab="Percentage of edges deleted",ylab="NCKD SS")

axis(1,at=1:11,labels=seq(0,20,2))
#axis(2,at=1:6,labels=seq(0,1.0,0.2))
lines(smat[2,],col="red",type="b",lty=2,lwd=2)
lines(smat[3,],col="mediumvioletred",type="b",lty=3,lwd=2)
lines(smat[4,],col="blue",type="b",lty=4,lwd=2)
legend("bottomleft",legend=c("CA-1","CA-4","DY-1","DY-3"), lty=c(1,2,3,4),lwd=2,pch=21,col=c("black","red","mediumvioletred", "blue"), ncol=4,bty="n",text.col=c("black","red","mediumvioletred","blue"), inset=0.01)

legend("bottomleft",legend=c("A-1","CA-1","DC-1","E10K-1"), lty=c(1,2,3,4),lwd=2,pch=21,col=c("black","red","mediumvioletred", "blue"), ncol=4,bty="n",cex=1.0,text.col=c("black","red","mediumvioletred","blue"), inset=0.01)


inputfile <- "./kshelltest - 0.4 .txt"
datafile =read.csv(inputfile, sep = "")
g=graph.data.frame(datafile, directed = FALSE)
g
plot(g)


