# new chloro pca

library(devtools)
library(QstFstComp)
library(psych) # for pairs.panels
library(car) # for dataEllipse
library(MASS) # for boxcox transformation
library(adegenet)
library(hierfstat)
library(pegas)

setwd("~/Dropbox/QSTFST/chloro/")

mydata <- read.csv("chloro_df.csv")
head(mydata)
data <- mydata[,2:6]
is.data.frame(data)
chloro <- df2genind(data, missing=NA, pop=mydata[,7], ploidy=1, ind.names=(mydata[,1]))

chloro.na <- na.replace(chloro, met=0)
summary(mydata.na)
chloro.na

pca1 <- dudi.pca(chloro.na, scannf=FALSE,scale=FALSE)
temp<- as.integer(pop(chloro.na))
plot(pca1$li, col=temp, pch=19)

chl <- cbind(temp, pca1$li$Axis1, pca1$li$Axis2)
chlx <- tapply(chl[,2], chl[,1], mean)
chly <- tapply(chl[,3], chl[,1], mean)

pdf("chloro.pdf", height=5, width=5)
par(mar=c(5,5,4,2), oma=c(0,0,0,0))
plot(chlx, chly, xlim=c(-4.5,3.5), ylim=c(-3.5, 4.5), pch=19,xlab="Principle Component 1", 
     ylab="Principle Component 2", main="Chloroplast Microsatellite PCA")
dataEllipse(pca1$li$Axis1, pca1$li$Axis2, groups=as.factor(temp), levels=.95, 
            plot.points=FALSE, col=rep("#d3d3d3",length(mydata[,1])), 
            group.label=rep("",length(mydata[,1])),center.pch=FALSE)
points(chlx, chly,pch=19)
dev.off()