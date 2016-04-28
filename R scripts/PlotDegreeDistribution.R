library('igraph')
setwd("~/graph-algos/NetworkComparison/data")
options(scipen=999)   #to supress scientific notation for printing 
#options(scipen=0)
datafile =read.csv("./Oregon/oregon1_010331.txt", sep = "")
datafile =read.csv("./arXiv/CA-AstroPh.txt", sep = "")
datafile = read.csv("./RealWorld/DBLP-CIKM.txt")
datafile = read.csv("./RewiredGraphs/ER/ER10K-0.0.txt", sep="")

g=graph.data.frame(datafile, directed = FALSE)


g
d <- degree.distribution(g, cumulative=T)
plot(d, xlab="node degree", cex=0.1, log="xy")
lines(degree.distribution(g))

# plot and fit the power law distribution
fit_power_law = function(graph) {
  # calculate degree
  d = degree(graph, mode = "all")
  dd = degree.distribution(graph, mode = "all", cumulative = FALSE)
  degree = 1:max(d)
  probability = dd[-1]
  # delete blank values
  nonzero.position = which(probability != 0)
  probability = probability[nonzero.position]
  degree = degree[nonzero.position]
  reg = lm(log(probability) ~ log(degree))
  cozf = coef(reg)
  power.law.fit = function(x) exp(cozf[[1]] + cozf[[2]] * log(x))
  alpha = -cozf[[2]]
  R.square = summary(reg)$r.squared
  print(paste("Alpha =", round(alpha, 3)))
  print(paste("R square =", round(R.square, 3)))
  # plot
#  plot(probability ~ degree, log = "xy", xlab = "Degree (log)", ylab = "Probability (log)", 
#       col = 1, main = "Degree Distribution")
  plot(probability ~ degree, cex=0.3,log = "xy", xlab = "Degree (log)", ylab = "Probability (log)", col = 1)

#curve(power.law.fit, col = "red", add = T, n = length(d))
}

fit_power_law(g)


cumhist <- function(x, plot=FALSE, ...){
  h <- hist(x, plot=FALSE, ...)
  h$counts <- cumsum(h$counts)
  h$density <- cumsum(h$density)
  h$itensities <- cumsum(h$itensities)
  
  if(plot)
    plot(h)
  h
}



plot_distribution = function (cs, cumulative = FALSE, ...) 
{
  #hi <- hist(cs, -1:max(cs), plot = FALSE)$density
  hi <- hist(cs, -1:max(cs), plot = FALSE)$intensities
  if (!cumulative) {
    res <- hi
  }
  else {
    res <- rev(cumsum(rev(hi)))
  }
  res
}

plot_distribution = function(d) {
  #dd = hist(d,-1:max(d), plot = FALSE)$intensities
  dd = hist(d,breaks=seq(0,1.0,0.05),plot = FALSE)$intensities
  degree = 1:max(d)
  probability = dd[-1]
  # delete blank values
  nonzero.position = which(probability != 0)
  probability = probability[nonzero.position]
  degree = degree[nonzero.position]
  plot(probability ~ degree, cex=0.3,log = "xy", xlab = "Ci (log)", ylab = "Probability (log)", col = 1)
  #curve(power.law.fit, col = "red", add = T, n = length(d))
}

data = read.csv("./RewiredGraphs/AS-ci", header=F)
d = as.numeric(data[1,])
cd <- cumhist(d)$intensities
cd
plot(cd, log="xy")
plot(cd, xlab="Local Clustering Coefficient (c)", ylab="P(x <= c)", axes=F)
axis(1,at=0:20,labels=seq(0,1.0,0.05))
axis(2,at=1:6,labels=seq(0,1.0,0.2))
plot_distribution(d1)
dd = hist(d,breaks=seq(0,1.0,0.05),plot = FALSE)$intensities
dd
