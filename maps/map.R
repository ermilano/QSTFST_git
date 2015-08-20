## Maping maps
library(maps)
library(mapdata)

setwd("~/Dropbox/QSTFST/maps/")
mydata <- read.csv("qstpopmerge_fst_map.csv")
data(stateMapEnv)


fst <- mydata[which(mydata$type=="Fst"),]
qst <- mydata[which(mydata$type=="Qst"),]

north <- mydata[which(mydata$lat>39),]
fst.n <- north[which(north$type=="Fst"),]
qst.n <- north[which(north$type=="Qst"),]

north.zoom <- mydata[which(mydata$lat>39.4),]
fst.nz <- north.zoom[which(north.zoom$type=="Fst"),]
qst.nz <- north.zoom[which(north.zoom$type=="Qst"),]

## Full
pdf("qstfstmap_full.pdf")
map('county', "Colorado",xlim=extendrange(mydata$long, f=1), ylim=extendrange(mydata$lat, f=.1))
points(qst$long, qst$lat, pch=15, cex=1, col=as.character(qst$color)) 
points(fst$long, fst$lat, pch=19, cex=1, col=as.character(fst$color))
text(fst$long, fst$lat, labels=fst$site,pos=2, cex=.5) 
text(qst$long, qst$lat, labels=qst$site,pos=4, cex=.5) 
dev.off()

## North
pdf("qstfstmap_north.pdf")
map('county', "Colorado",xlim=extendrange(north$long, f=.3), ylim=extendrange(north$lat, f=.05))
points(qst.n$long, qst.n$lat, pch=15, cex=1, col=as.character(qst.n$color))
points(fst.n$long, fst.n$lat, pch=19, cex=1, col=as.character(fst.n$color)) 
text(fst.n$long, fst.n$lat, labels=fst.n$site,pos=2, cex=.8) 
text(qst.n$long, qst.n$lat, labels=qst.n$site,pos=4, cex=.8) 
dev.off()

## North Zoom
pdf("qstfstmap_northzoom.pdf")
map('county', "Colorado",xlim=extendrange(north.zoom$long, f=.3), ylim=extendrange(north.zoom$lat, f=.05))
points(qst.nz$long, qst.nz$lat, pch=15, cex=1, col=as.character(qst.nz$color))
points(fst.nz$long, fst.nz$lat, pch=19, cex=1, col=as.character(fst.nz$color)) 
text(fst.nz$long, fst.nz$lat, labels=fst.nz$site,pos=2, cex=.8) 
text(qst.nz$long, qst.nz$lat, labels=qst.nz$site,pos=4, cex=.8)
dev.off()


#############################
## Terrain maps
library(RgoogleMaps)
long<- range(mydata$long)
lat<-range(mydata$lat)


head(mydata)
markerdat <- cbind(mydata[,3:4])
names(markerdat) <- c("lat", "lon")

bb <- qbbox(lat = markerdat[,"lat"], lon = markerdat[,"lon"])

terrain_close <- GetMap.bbox(bb$lonR, bb$latR, destfile= "terrtest.png", 
                             markers= markerdat, maptype="terrain")

lat <- c(37, 40.5)
long <- c(-102, -109)
center = c(mean(lat), mean(long))
zoom <- 7
test <- GetMap.bbox(long, lat, destfile= "co.png", maptype="hybrid", GRAYSCALE = TRUE)



library(RgoogleMaps)
lat <- c(48,64) #define our map's ylim
lon <- c(-140,-110) #define our map's xlim
center = c(mean(lat), mean(lon))  #tell what point to center on
zoom <- 5  #zoom: 1 = furthest out (entire globe), larger numbers = closer in
terrmap <- GetMap(center=center, zoom=zoom, maptype= "terrain", destfile = "terrain.png") 
# lots of visual options, just like google maps: maptype = c("roadmap",
# "mobile", "satellite", "terrain", "hybrid", "mapmaker-roadmap",
# "mapmaker-hybrid")


#############################
# admixture map 
#Initialize packages
# from http://www.molecularecologist.com/2014/11/admixture-maps-in-r-for-dummies/

setwd("~/Dropbox/QSTFST/maps/")
library(maps)
library(plotrix)

gps<-read.csv("fst.csv") #Read input files
gps.fix<-read.csv("piemapfix.csv")

# === start tv 
gps.lines <- gps
gps.lines[1,2] <- gps.lines[1,2] + 1
gps.lines[2,2] <- gps.lines[2,2] - 1
gps.lines[5,3] <- gps.lines[5,3] - 0.5

# pdf("kmeansmap.pdf")
map('county', "Colorado",xlim=extendrange(gps$Lon, f=1), ylim=extendrange(gps$Lat, f=.1))
map.axes() #Add axes
#To add admixture plots – here I used K = 2.
for(it in 1:nrow(gps.fix)) {
 #it <- 1
  lines(c(gps.lines[it,3], gps.lines[it,3]), c(gps[it,2], gps[it,2]), col = "red" )
}

for (x in 1:nrow(gps)) { 
  floating.pie(gps$Lon[x],gps$Lat[x], #You can modify your loop to reflect this
               c(admix$K1[x],admix$K2[x]),radius=.08,         
               col=c("red","blue") )
}

# ======== end tv

admix<-read.csv("kmeansassign.csv")

# pdf("kmeansmap.pdf")
map('county', "Colorado",xlim=extendrange(gps$Lon, f=1), ylim=extendrange(gps$Lat, f=.1))
map.axes() #Add axes
#To add admixture plots – here I used K = 2.
for (x in 1:nrow(gps)) { 
  floating.pie(gps$Lon[x],gps$Lat[x], #You can modify your loop to reflect this
               c(admix$K1[x],admix$K2[x]),radius=.08,         
               col=c("red","blue") )
}
title(main="DAPC K-means population assignments", xlab="Longitude", ylab="Latitude",  cex.lab=0.75)


# plot just pop number to figure out pie placement
ylim=c(36.9, extendrange(gps$Lat, f=.1)[2])
map('county', "Colorado",xlim=extendrange(gps$Lon, f=1), ylim=ylim)
map.axes() #Add axes
text(gps$Lon,gps$Lat, gps[,1])
text(gps.fix$pie.lon, gps.fix$pie.lat, gps.fix[,1])

# testing
pdf("fixpie.pdf")
map('county', "Colorado",xlim=extendrange(gps.fix$Lon, f=1), ylim=extendrange(gps$Lat, f=.1))
map.axes() #Add axes
#To add admixture plots – here I used K = 2.
for(i in 1:nrow(gps.fix)) {
  #it <- 1
  lines(c(gps.fix[i,3], gps.fix[i,5]), c(gps.fix[i,2], gps.fix[i,4]), col = "black", lwd=1.5, pch=20 )
}
points(gps.fix$Lon, gps.fix$Lat, pch=20)
for (x in 1:nrow(gps.fix)) { 
  floating.pie(gps.fix$pie.lon[x],gps.fix$pie.lat[x], #You can modify your loop to reflect this
               c(admix$K1[x],admix$K2[x]),radius=.08,         
               col=c("red","blue") )
}
title(main="DAPC K-means population assignments", xlab="Longitude", ylab="Latitude",  cex.lab=0.75)
dev.off()


par(fig=c(0, 1, 0, 1), new=TRUE)
map('state', 'Colorado', lwd=3)
par(fig=c(0, 1, 0, 1), new=TRUE)
pdf("colorado.pdf")
map('state', "Colorado", lwd=2)
rect(extendrange(gps.fix$Lon, f=1)[1], 36.95, extendrange(gps.fix$Lon, f=1)[2], 
     extendrange(gps.fix$Lat, f=.1)[2], border="red", lwd=5)
dev.off()

title(main="DAPC K-means population assignments", xlab="Longitude", ylab="Latitude",  cex.lab=0.75)



#################
## testing inset
library(lattice)
library(maptools)
library(ggplot2)
# aa <- rep(1,5); ab <- c(2,3,1,1,1)
# mat <- rbind(aa,aa,aa,ab)
# layout(mat)
# layout(mat)
# map('county', "Colorado",xlim=extendrange(gps$Lon, f=1), ylim=extendrange(gps$Lat, f=.1))
# map.axes()
# par(mar=c(0,0,0,0))
# map('worldHires', region='USA:Colorado', xlim=c(-104.5, -103.5), col="red")

mat <- matrix(c(1,2),2, byrow=TRUE)
layout(mat, c(10, 10), c(2,2))
par(mar=c(0.05,0.05,0.05,0.05))
map('county', "Colorado",xlim=extendrange(gps$Lon, f=1), ylim=extendrange(gps$Lat, f=.1))
par(mar=c(0.05,0.05,0.05,0.05))
map('state', region='Colorado', col="red")



#################
# trying to get all sites onto map
library(maps)

setwd("~/Dropbox/QSTFST/maps/")
gps.all <- read.csv("newmap.csv")
fst <- gps.all[which(gps.all$type=="Fst"),]
qst <- gps.all[which(gps.all$type=="Qst"),]

## Full
# pdf("qstfstmap_allpops.pdf")
map('county', "Colorado",xlim=extendrange(gps$Lon, f=1), ylim=extendrange(gps$Lat, f=.1))
points(qst$long, qst$lat, pch=15, cex=1, col=as.character(qst$color)) 
points(fst$long, fst$lat, pch=19, cex=1, col=as.character(fst$color))
text(fst$long, fst$lat, labels=fst$num,pos=2, cex=.5) 
text(qst$long, qst$lat, labels=qst$num,pos=2, cex=.5) 
map.axes()
title(main="All collection sites", xlab="Longitude", ylab="Latitude",  cex.lab=0.75)
# dev.off()

co.map <- get_map("colorado springs", #first argument is basically querrying google maps 
                  maptype="terrain", #you can specify roadmaps or terrain
                  zoom=7) #play with zoom to get the extent you want 1-21
co.map.r <- get_map("colorado springs", #first argument is basically querrying google maps 
                  maptype="road", #you can specify roadmaps or terrain
                  zoom=7)
co.map.h <- get_map("colorado springs", #first argument is basically querrying google maps 
                    maptype="terrain-background", #you can specify roadmaps or terrain
                    zoom=7)
pdf("allpointsco_nooverlap.pdf")
gps.all <- gps.all[c(1:22,35:39),]
ggmap(co.map.r)+
  geom_point(data=gps.all, aes(x=long, y=lat), fill=gps.all$color,
             alpha=.8,size=5, shape=21) +
  
  ggtitle("All collection sites") +
 # geom_text(data=gps.all, aes(x=long, y=lat, label=num), size=3, color="white") + #adding IDs
  ylim(36.9, 40.25) +
  xlim(-107,-103.5) +
  xlab("Longitude")+
  ylab("Latitude")
dev.off()

pdf("allpointsco_new.pdf")
ggmap(co.map.h)+
  geom_point(data=gps.all, aes(x=long, y=lat), fill=gps.all$color,
             alpha=.8,size=7.5, shape=22) +
  ggtitle("All collection sites") +
  #geom_text(data=gps.all, aes(x=long, y=lat, label=num), size=3, color="white") + #adding IDs
  ylim(36.9, 40.25) +
  xlim(-107,-103.5) +
  xlab("Longitude")+
  ylab("Latitude")
dev.off()





