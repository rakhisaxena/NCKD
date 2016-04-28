library('igraph')
setwd("~/graph-algos/NetworkComparison/data")


labelCol <- function(x) {
  labelnames = c("CA-", "B1","A-", "FW","M-","WC","F1", "E1", "B1", "W2")
  colors <- rainbow(length(labelnames))
  if (is.leaf(x)) {
    ## fetch label
    labelstring <- attr(x, "label") 
    label <- substr(labelstring,1,7)
    ## set label color to red for A and B, to blue otherwise
    #attr(x, "nodePar") <- list(lab.col=ifelse(label %in% c("A", "B"), "red", "blue"))
    if (grepl(labelnames[1],label,)){
      attr(x, "nodePar") <- list(lab.col=colors[1])
    } else if (grepl(labelnames[2],label)) {
      attr(x, "nodePar") <- list(lab.col=colors[2])
    }else if (grepl(labelnames[3],label)){
      attr(x, "nodePar") <- list(lab.col=colors[3])
    } else if (grepl(labelnames[4],label)) {
      attr(x, "nodePar") <- list(lab.col=colors[4])
    } else if (grepl(labelnames[5],label)) {
      attr(x, "nodePar") <- list(lab.col=colors[5])
    } else if (grepl(labelnames[6],label)) {
      attr(x, "nodePar") <- list(lab.col=colors[6])
    } else if (grepl(labelnames[7],label)) {
      attr(x, "nodePar") <- list(lab.col=colors[7])
    } else if (grepl(labelnames[8],label)) {
      attr(x, "nodePar") <- list(lab.col=colors[8])
    } else if (grepl(labelnames[9],label)) {
      attr(x, "nodePar") <- list(lab.col=colors[9])
    } else {
      attr(x, "nodePar") <- list(lab.col="black")
    }
    
    #attr(x, "nodePar") <- labelcolor
  }
  return(x)
}


data <- as.matrix(read.csv(file="./NetSimileGraphs/NS-Data-New_NodeDistanceMatrix.txt", sep=" ",header=T))

#distanc <- 1 - as.dist(data)   # if you read in a similarity matrix
distanc <- as.dist(data) 
hc <- hclust(distanc) 
dendro <- as.dendrogram(hc)
#plot(dendro)
d <- dendrapply(dendro,  labelCol)
plot(d)

classificationNode <- cutree(hc,k=9)
classificationActual <- c(rep(1,5),rep(2,4),rep(3,5),rep(4,4),rep(5,4),rep(6,6),rep(7,3),rep(8,4),rep(9,3))
names(classificationActual) <- names(classificationNode)
compare(classificationActual, classificationNode, method=c("nmi"))

