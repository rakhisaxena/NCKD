library(nnet)
library('igraph')
setwd("~/graph-algos/NetworkComparison/data")
setwd("~/graph-algos/NetworkComparison/Python-Programs/TestData/data")
datafile <- "./DBLP/icdm-union-ids/2013-ids.txt"
datafile <- "./TestEdgeMatrix/kshelltest.txt"
datafile <- "./TestEdgeMatrix/kedgetest.txt"
datafile <- "../Python-Programs/TestData/kedge/kedgetest"
datafile <- "../Python-Programs/TestData/fpgraph/FingerprintGraphV2"
datafile <- "./MetaboliteGraphs/MakeGraphsSimple/MGSimpleData/A-AG.smet"
datafile <- "./MetaboliteGraphs/Ecolimetabolite.csv"
datafile <- "./MetaboliteGraphs/Bacillusmetabolite.csv"
datafile <- "../Python-Programs/TestData/data/karate.csv"
datafile <- "./TrialGraphs/linetestgraph"
datafile <- "./NetSimileGraphs/EdgeDistributionPlots/CA/CA-data/CA-HepTh.txt"

datafile <- "../data/kedgetest.txt"
datafile <- "./karate.csv"

el =read.csv(datafile, sep = "", head=F)
g=graph.data.frame(el, directed = FALSE)
vcount(g)
ecount(g)
coreness <- graph.coreness(g)
maxCoreness <- max(coreness)
colors <- rainbow(max(coreness+1))
tkplot(g,vertex.color=colors[coreness+1], vertex.size=10, layout=layout.kamada.kawai)



plot(g)
l1 <- line.graph(g)
plot(l1)
l2 <- line.graph(l1)
plot(l2)
l3 <- line.graph(l2)
plot(l3)
lp <- label.propagation.community(g)
tkplot(g)

degree <- degree(g)

coreness <- graph.coreness(g)
maxCoreness <- max(coreness)
# if you just need to know the vertices and not to build the subgraph 
# you can use this variable
verticesHavingMaxCoreness <- which(coreness == maxCoreness) 
kcore <- induced.subgraph(graph=g,vids=verticesHavingMaxCoreness)

plot(kcore, 
     vertex.label=get.vertex.attribute(kcore,name='vert.names',index=V(kcore)))

g2 <- graph.complementer(g)
coreness2 <- graph.coreness(g2)

mycoreness <- c(rep(4,5),rep(3,6),1,2,1,2,rep(1,2),rep(2,2),rep(1,4),2,rep(1,8))
mycoreness <- c(5,5,5,5,5,2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
names(mycoreness) <- names(coreness)

sg <- simplify(g)
coreness <- graph.coreness(g)
maxCoreness <- max(coreness)
colors <- rainbow(max(coreness+1))
plot.igraph(g,layout=layout.kamada.kawai, vertex.color=colors[coreness])

tkplot(g,vertex.color=colors[coreness+1], vertex.size=10, layout=layout.kamada.kawai)
plot.igraph(g,vertex.size=5,vertex.label=NA,layout=layout.kamada.kawai, vertex.color=colors[coreness])
