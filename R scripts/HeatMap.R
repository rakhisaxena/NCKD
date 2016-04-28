library(gplots)
library(RColorBrewer)
library(lattice) 
library(raster)

data <- read.csv("./heatmapTestdata.csv",  sep=",")
data <- read.csv("./ER_NSdists.txt", sep="")
data <- read.csv("./NSvsKC-Sorted-KCore-DistanceMatrix.csv",  sep="", header=T)
data <- read.csv("./NetSimileVsKCore-Outputs/Test-1-NetSimile-DistanceMatrix.txt",  sep="", header=T)


mat_data <- as.matrix(data)
png("./PlotsForPaper/NSvsKC-KCore-HeatMap2.png",    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size
#heatmap(mat_data, symm=TRUE, Rowv=NA, Colv=NA, 
#        col = heat.colors(256), scale="column", margins=c(5,10))
levelplot(mat_data[1:ncol(mat_data),ncol(mat_data):1], 
          xlab=NULL, ylab=NULL, col.regions = heat.colors(256), 
          scale=list(x=list(rot=90, cex=.5),y=list(cex=.5)),margins=c(5,10))

levelplot(mat_data[1:ncol(mat_data),ncol(mat_data):1], 
          xlab=NULL, ylab=NULL, col.regions = heat.colors(256), 
          scale=list(x=list(rot=90, cex=1),y=list(cex=1)),margins=c(5,10))
dev.off()

#ncol(data)
#nrow(data)
#rnames <- data[1,]                            # assign labels in column 1 to "rnames"
#mat_data <- data.matrix(data[,2:ncol(data)])  # transform column 2-5 into a matrix
#mat_data <- as.matrix(data[,2:ncol(data)])

#rownames(mat_data) <- rnames                  # assign row names 
$rownames(mat_data)


levelplot(mat_data[1:ncol(mat_data),ncol(mat_data):1], 
          xlab=NULL, ylab=NULL, scale=list(x=list(rot=90)))

# creates a own color palette from red to green
my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)

# (optional) defines the color breaks manually for a "skewed" color transition
col_breaks = c(seq(-1,0,length=100),  # for red
               seq(0,0.8,length=100),              # for yellow
               seq(0.8,1,length=100))              # for green

# creates a 5 x 5 inch image
png("./PlotsForPaper/NSvsKC-KCore-HeatMap3.png",    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

heatmap.2(mat_data,
          #cellnote = mat_data,  # same data set for cell labels
          #main = "Correlation", # heat map title
          symm=TRUE,             #tells whether it is symmetrical or not
          #symkey=TRUE,symbreaks=TRUE, # Tells it to keep things symmetric around 0
          scale="none",         # Tells it not to scale the data
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(15,15),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="row",     # only draw a row dendrogram
          #distfun =function(c) {as.dist(c)}, hclustfun=hclust, 
          Colv="NA")            # turn off column clustering

dev.off()               # close the PNG device
# don't name an object dist --- clashes with dist function
rm(list=ls()) #will remove ALL objects 

mat_data <- matrix( c(0,1,1,1,1,0,1,1,1,1,0,1,1,1,1,0), nrow=4, ncol=4, byrow=FALSE)
mat_data <- matrix( c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1), nrow=4, ncol=4, byrow=FALSE)
mat_data
heatmap(mat_data, symm=TRUE, Rowv=NA, Colv=NA, 
        col = heat.colors(256), scale="column", margins=c(5,10))
heatmap(mat_data, symm=TRUE, Rowv=NA, Colv=NA, 
        col=my_palette, scale="column", margins=c(5,10))

rd <- dist(data)
rc<-hclust(rd)
cd<-dist(t(data))
cc<-hclust(cd)
heatmap(data, Rowv=as.dendrogram(rc), Colv=as.dendrogram(cc))