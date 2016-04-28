library(igraph)
setwd("/home/stratuser/graph-algos/NetworkComparison/data/")
datafile =read.csv("./AllGraphs/BA10K-1.txt", sep = "")
datafile =read.csv("./AllGraphs/BAGraphs/BA4-1K-1.txt", sep = "")
datafile = read.csv("./AllGraphs/ERGraphs/Rewired/ER3-1K-r05.txt")
datafile = read.csv("../Python-Programs/TestData/kedge/kedgetest", sep="",head=F)

g=graph.data.frame(datafile, directed = FALSE)
tkplot(g, vertex.size=3,vertex.label=NA, layout=layout.kamada.kawai)

tkplot(g, layout=layout.fruchterman.reingold)
plot(degree.distribution(g), xlab="node degree", cex=0.1)
lines(degree.distribution(g))
degree(g)
tkplot(g, vertex.size=3,vertex.label=NA, layout=layout.fruchterman.reingold)
tkplot(g, vertex.size=3, layout=layout.fruchterman.reingold)

tkplot(barabasi.game(1000,1, directed=F), vertex.size=3,vertex.label=NA, layout=layout.kamada.kawai)

datafile =read.csv("./AllGraphs/ERGraphs/ER3-1K-1.txt", sep = "")
g=graph.data.frame(datafile, directed = FALSE)
g
plot(degree.distribution(g), xlab="node degree", cex=0.1)
lines(degree.distribution(g))

g1 <- rewire.edges(g, p=0.5)
write.graph(g, "./AllGraphs/ERGraphs/Rewired/ER3-1K-r50.txt" ,format=c("edgelist"))
plot(degree.distribution(g1), xlab="node degree", cex=0.1)
lines(degree.distribution(g1))

betweenness(g, v=V(g))



