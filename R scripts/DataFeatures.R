library('igraph')
setwd("~/graph-algos/NetworkComparison/data")
file <- "../Python-Programs/TestData/data/myAMiner-CoauthorWeighted.txt"
plfit <- function(g) {
  d <- degree(g)
  fit <- power.law.fit(d+1, 10)
  return(fit$alpha)
}

metrics <- data.frame(
  deg=degree(graph),
  bet=betweenness(graph),
  clo=closeness(graph),
  eig=evcent(graph)$vector,
  tra=transitivity(graph,type=c("local"),
  alpha = plfit(graph)
    )
)

filenames <- list.files( "./MetaboliteGraphs/BarabasiGraphs/AEData/NCKD/AEGraphs", full.names=T, recursive=FALSE)
filenames <- list.files("../data/FingerPrintGraphs/LBD-Data", full.names=T, recursive=FALSE)
filenames <- list.files("../data/FingerPrintGraphs/NS-Data/", full.names=T, recursive=FALSE)

length(filenames)
file <- 
properties <- list()
for (file in filenames) {
  aae = read.csv(file, sep = "",header=F)
  gr = graph.data.frame(aae, directed = FALSE)
  g = simplify(gr)
  properties <- rbind(properties,c(Name=file,
                                   N=length(V(g)), 
                                   E=length(E(g)), 
                                   dia=diameter(g),
                                   C=clusters(g)$no,
                                   GCC=transitivity(g, type="global"),
                                   alpha=plfit(g))

                                   )
}
# d=graph.density(g),
#m=max(graph.coreness(g))
properties
print(properties, quote=F,row.names=FALSE)
sink("./MetaboliteGraphs/BarabasiGraphs/AEData/NCKD/AEGraphs-Features")
print(properties, quote=FALSE,row.names=FALSE)
sink()


file <- "./ScalabilityTest/BA/BA-300K"
aae = read.csv(file, sep = "",header=F)
g = graph.data.frame(aae, directed = FALSE)
length(V(g))
length(E(g)) 



file <- "./AllGraphs/CA-GrQc.txt"
aae = read.csv(file, sep = "")
g = graph.data.frame(aae, directed = FALSE)
clusters(g)$no
