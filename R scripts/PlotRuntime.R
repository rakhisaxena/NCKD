datafile <- "./NetSimileGraphs/runtimeNSvsNCKD"
data <- read.csv(datafile,header=T)
plot(data$NetSimile,type="b",xpd=T,lwd=2, xaxt="n",ylim=c(0,max(data$NetSimile)),col="blue", xlab="Dataset Size",ylab="Runtime(seconds)")
axis(1,at=1:length(data$Dataset),labels=data$Dataset,las=1, cex.lab=)
lines(data$NCKD,col="red",type="b",lwd=2)
legend("topleft",legend=c("NetSimile","NCKD"), lty=1,lwd=2,pch=21,col=c("blue","red"), ncol=2,bty="n",cex=0.8,text.col=c("blue","red"), inset=0.01)


data <- structure(list(data$NCKD, data$NetSimile)) # , .Names = data$Dataset)
colours <- c("blue","green")
barplot(as.matrix(data),axis.lty=1,cex.lab = 1.5, beside=TRUE, col=colours)
legend("topleft", c("KarateClub Data","Football Data"), cex=0.8, bty="n", fill=colours)



x <- c(data$NetSimile)
y <- c(data$NCKD)
height <- rbind(x, y)     # create a two row matrix with x and y
colours <- c("blue","green")
# Use height and set 'beside = TRUE' to get pairs , save the bar midpoints in 'mp'
mp <- barplot(height, ylim = c(seq(0,max(data$NetSimile), by=50)), axis.lty=1, beside = TRUE, names.arg = data$Dataset, col=colours,las=2)
# Draw the bar values above the bars
text(mp, height, labels = format(height, 4),pos = 3, cex = .75)
