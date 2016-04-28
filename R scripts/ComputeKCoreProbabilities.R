library('igraph')
library(data.table)
setwd("~/graph-algos/NetworkComparison/data")
options(scipen=999)   #to supress scientific notation for printing 
#options(scipen=0)

inputgraphDirectory <- "~/graph-algos/NetworkComparison/data/RewiredGraphs/DBLP-C"
outputFile <- "../DBLP-C-Rewired-NetSimile-ProbDistn.txt"
setwd(inputgraphDirectory)
unlink(outputFile)   # delete previous output file
filenames <- list.files( ".", full.names=T, recursive=FALSE)
length(filenames)
alldata <- rbind(lapply(filenames, function(x) {
  #x <- "./NetSimileGraphs/EdgeDistributionPlots/CA/CA-data/CA-4"
  aae = read.csv(x, sep = "")
  g<- graph.data.frame(aae, directed = FALSE)
  coreness <- graph.coreness(g)
  #max(coreness)
  kshellprobability <- hist(coreness, -1:max(coreness), plot = FALSE)$intensities
  fn <- gsub("./|.txt|.dat|.csv" , "" , x)
  result <- c(fn, kshellprobability)
}))
lapply(alldata, write,outputFile , append=TRUE, sep=",", ncol=1000)


x <- "./CitationNetworks/Cit-Patents.txt
x <- "./TestGraphs/kshelltest.txt""
aae = read.csv(x, sep = "")
g<- graph.data.frame(aae, directed = FALSE)
coreness <- graph.coreness(g)
kshellprobability <- hist(coreness, -1:max(coreness), plot = FALSE)$intensities
fn <- gsub("./|.txt|.dat|.csv" , "" , x)
result <- c(fn, kshellprobability)
tr <- t(result)
outputFile <- "./NetSimileVsKCore-Outputs/AllGraphs_KCoreFeatures.txt"
lapply(tr, write,outputFile , append=TRUE, sep=",", ncol=2000)


#write.table(t(alldata),col.names = F, row.names = filenames, quote = F, "./KShellProbabilityDistributionTest.txt")
#sink("./KShellProbabilityDistributionTest.txt")
#alldata
#sink()