## map for merged populations
library(maps)
setwd("~/Dropbox/QSTFST/match/")
map.data <- read.csv("matchedpopsmap.csv")

head(map.data)

fst <- map.data[which(map.data$type=="fst"),]
qst <- map.data[which(map.data$type=="qst"),]

north <- map.data[which(map.data$lat>39),]
fst.n <- north[which(north$type=="fst"),]
qst.n <- north[which(north$type=="qst"),]

north.zoom <- map.data[which(map.data$lat>39.4),]
fst.nz <- north.zoom[which(north.zoom$type=="fst"),]
qst.nz <- north.zoom[which(north.zoom$type=="qst"),]




##### Color by collection matched site #####

## Full
pdf("matchmap.pdf")
map('county', "Colorado",xlim=extendrange(map.data$long, f=2), ylim=extendrange(map.data$lat, f=.1))
points(qst$long, qst$lat, pch=15, cex=1, col=rainbow(qst$number)) 
points(fst$long, fst$lat, pch=17, cex=1, col=rainbow(fst$number))
text(fst$long, fst$lat, labels=fst$site,pos=4, cex=.5) 
text(qst$long, qst$lat, labels=qst$site,pos=2, cex=.5) 
#dev.off()

## North
#pdf("matchmap_north.pdf")
map('county', "Colorado",xlim=extendrange(north$long, f=.3), ylim=extendrange(north$lat, f=.05))
points(qst.n$long, qst.n$lat, pch=15, cex=1, col=rainbow(qst.n$number))
points(fst.n$long, fst.n$lat, pch=17, cex=1, col=rainbow(qst.n$number)) 
text(fst.n$long, fst.n$lat, labels=fst.n$site,pos=4, cex=.8) 
text(qst.n$long, qst.n$lat, labels=qst.n$site,pos=2, cex=.8) 
#dev.off()

## North Zoom
#pdf("matchmap_northzoom.pdf")
map('county', "Colorado",xlim=extendrange(north.zoom$long, f=.3), ylim=extendrange(north.zoom$lat, f=.05))
points(qst.nz$long, qst.nz$lat, pch=15, cex=1, col=rainbow(qst.nz$number))
points(fst.nz$long, fst.nz$lat, pch=17, cex=1, col=rainbow(fst.nz$number)) 
text(fst.nz$long, fst.nz$lat, labels=fst.nz$site,pos=4, cex=.8) 
text(qst.nz$long, qst.nz$lat, labels=qst.nz$site,pos=2, cex=.8) 
dev.off()

