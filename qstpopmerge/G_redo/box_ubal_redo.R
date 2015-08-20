##### Checking all analyses for discrepenses 
# August 4, 2015
#####

#########################################################
#### Load libraries and data


## New Qst-Fst comparison tool from Whitlock
## https://github.com/kjgilbert/QstFstComp
## install.packages("devtools")
## install_github("kjgilbert/QstFstComp")

library(devtools)
library(QstFstComp)
library(psych) # for pairs.panels
library(car) # for dataEllipse
library(MASS) # for boxcox transformation
library(adegenet)
library(hierfstat)
library(pegas)

setwd("~/Dropbox/QSTFST/qstpopmerge/G_redo")

MarkerInfo <- read.csv("geno.csv")
TraitInfo <- read.csv("qstmerge_unbal.csv")


#########################################################
#### BoxCox transformations

box <- TraitInfo[1:10] # raw phenotype data
head(box)

# Boxcox steps:
tmp <- boxcox(box$Corollawidth ~ box$dam) # calculate lambda
lam <- tmp$x[max(tmp$y)==tmp$y] # optimal lambda 
new <- box$Corollawidth^lam # apply transformation

# boxcox function
funbox <- function(trait, dam){
  tmp <- boxcox(trait ~ dam)
  lam <- tmp$x[max(tmp$y)==tmp$y]
  return(lam)
}

par(mfrow=c(2,4))
lamda <- as.vector(apply(box[3:10], 2, funbox, dam=box$dam))
lamda[3]
head(box)

b1 <- box$Corollawidth^lamda[1]
b2 <- box$Corollalength^lamda[2]
b3 <- box$Petallength^lamda[3]
b4 <- box$color492^lamda[4]
b5 <- box$nectaramt^lamda[5]
b6 <- box$nectarconc^lamda[6]
b7 <- box$Root^lamda[7]
b8 <- box$fdate^lamda[8]

all.box <- cbind(box[,1:2], b1, b2, b3, b4, b5, b6, b7, b8)

par(mfrow=c(2,4))
apply(all.box[3:10], 2, hist)

all.box1 <- scale(all.box[3:10])
box.scale <- cbind(all.box[1:2], all.box1)
names(box.scale) <- names(box)

par(mfrow=c(2,4))
apply(box.scale[3:10], 2, hist)
apply(box.scale[3:10], 2, qqnorm)

box.scale.na <- na.exclude(box.scale)

## save files
write.csv(box.scale, "box.scale.csv")
write.csv(box.scale.na, "box.scale.na.csv")

box.scale.na <- read.csv("box.scale.na.csv")
box.scale <- read.csv("box.scale.csv")
names(box.scale.na) <- c("id","pop","dam","Corollawidth","Corollalength","Petallength",
                        "color492","nectaramt","nectarconc","Root","fdate")
names(box.scale) <- c("id","pop","dam","Corollawidth","Corollalength","Petallength",
                      "color492","nectaramt","nectarconc","Root","fdate")


#########################################################
#### Scatterplot matrix, correlations and PCA for ***floral traits only***

### Scatterplot matrix
pairs.panels(box.scale.na[4:9], pch=20, cex = 1)

### Basic PCA
# means and standard deviations
sapply(box.scale.na[4:9],mean) 
sapply(box.scale.na[4:9],sd) 

# Perforn the PCA on scaled data
x <- box.scale.na[4:9]
pca.x <- princomp(x,scale=TRUE, center=TRUE)
summary(pca.x)
names(pca.x)
screeplot(pca.x, type="barplot")

# Trait loadings on each PC
pca.x$loadings

### Barplots for trait loadings on the first 4 PCs
xname <- c("CW", "CL", "PL", "COL", "NA", "NC")
color <- c("darkred", "darkred", "darkred", "darkred", "yellow", "yellow")
summary(pca.x)
# pdf("pheno_pcs-flower.pdf")
par(mfrow=c(2,2), mar=c(3,3,3,1), oma=c(0,0,3,0))
barplot(pca.x$loadings[,1], names.arg=xname, col=color, main="PC1, 45.1%")
box(which="figure")
barplot(pca.x$loadings[,2], names.arg=xname, col=color, main="PC2, 17.3%")
box(which="figure")
barplot(pca.x$loadings[,3], names.arg=xname, col=color, main="PC3, 15.7%")
box(which="figure")
barplot(pca.x$loadings[,4], names.arg=xname, col=color, main="PC4, 9.7%")
box(which="figure")
title(main="Trait Loadings for Flower Phenotype PCA", cex.main=2, outer=TRUE)
# dev.off()
par(mfrow=c(1,1))

# Plot the PCA
# pdf("pheno_flower_pca_box.pdf")
plot(pca.x$scores[,1], pca.x$scores[,2], pch="", cex=2, xlim=c(-3, 6), ylim=c(-5.5, 5.5),
     xlab="Principle Component 1", ylab="Principle Component 2", 
     main="Phenotype PCA for Floral Traits") 
pops <- c("westchicagocreek","guanellapass","Golden","Smiths","Cuchara","noe",
          "hwy86","72med","72high","caribou_rd","springcreek","Lefthand","MRS",
          "wilkersonpass","TopCaribou","hwy105","Georgetown")
pops.code <- c("WCC","BLY","BLH","COI","CUC","NOE",
          "H86","GDR","R72","RVR","SPC","LFC","MRS",
          "WKP","RVR","H105","GTN")
newcol <- c("#C7D34F","#CA53D5","#D65033","#89C9CA","#56335B","#587C3D","#CB98C0",
            "#C6AB8E","#7B4627","#796ECD","#83D598","#CF9A3F","#C84792","#3F4239",
            "#6BD249","#B84B5C","#6286AE")
dataEllipse(pca.x$scores[,1], pca.x$scores[,2], groups=as.factor(box.scale.na$pop), levels=.95, 
            plot.points=FALSE, col=newcol, group.labels=pops.code)
# dev.off()

#biplot of first two principal components
biplot(pca.x)


#########################################################
#### Qst-Fst

# This script only runs one phenotype at a time.
# Need to create objects that pair each phenotype with "pop" and "dam" columns
trait.col <- list()
for (i in c(4:length(box.scale))){
  trait.col[[i]] <- box.scale[,c(2,3,i)]
}

# Name list items
names(trait.col) <- c("ID", "Pop", "Dam", "CorollaWidth",  "CorollaLength", "PetalLength", "Color", 
                      "NectarAmt",  "NectarConc", "Root", "Fdate")

# Loop through QstFstComp and save to list
qst.out <- list()
for (i in names(trait.col[4:length(trait.col)])){
  qst.out[[i]] <- QstFstComp(fst.dat = MarkerInfo, qst.dat = trait.col[[i]], numpops = 22, 
              breeding.design = "half.sib.dam", nsim = 100000, , output="concise")
  }


### Plotting the Qst-Fst distribution from the output file
# test <- read.table("/Users/elizabethmilano/Dropbox/QSTFST/QminusFvalues_2015-03-16_12-07-18.txt")
# hist(test[,1])


### Saving values
qst.out
Q.F<-as.data.frame(sapply(qst.out,function(x) as.numeric(x[[1]][1])))  # Calculated Qst-Fst 
pval <-as.data.frame(sapply(qst.out,function(x) as.numeric(x[[3]][3]))) #Two-tailed p value
fst<- as.data.frame(sapply(qst.out,function(x) as.numeric(x[[4]][1]))) #Fst
qst<- as.data.frame(sapply(qst.out,function(x) as.numeric(x[[5]][1]))) #Qst
qst.low <- as.data.frame(sapply(qst.out,function(x) as.numeric(x[[5]][2]))) #Qst lower bound
qst.up <- as.data.frame(sapply(qst.out,function(x) as.numeric(x[[5]][3]))) #Qst upper bound

results<-cbind(Q.F,pval,fst, qst, qst.low, qst.up)
colnames(results) <- c("Qst_Fst", "pval", "Fst", "Qst", "Qst CI lower", "Qst CI upper")
results
write.csv(results, "results.csv")



############################
#### Plotting Qst

results.all <- read.csv("results.csv")
results <- results.all[1:8,2:7]
rownames(results) <- results.all$X[1:8]

res.name <- c("Corolla Width",  "Corolla Length", "Petal Length", "Color", "Nectar Volume",  
              "Nectar Volume", "Root Diameter", "Flowering Date")
col <- c("darkred", "darkred", "darkred", "darkred", "yellow", "yellow", "darkolivegreen3", 
         "darkolivegreen3")
pdf("qst_bar.pdf")
par(mfrow=c(1,1), mar=c(8,5,4,1))
# barplot(results[,4], names.arg=res.name,las=2, cex.names=1, col=col, main="test")
# error bars
# could not get the confidence intervals out of list format so I created a work around
barx <- barplot(results[,4], names.arg=res.name,las=2, cex.names=1, col=col, 
                main="Qst with 95% Confidence Interval", ylab="Qst", ylim=c(0,1.3))
low <- results$Qst.CI.lower
up<- results$Qst.CI.upper
arrows(barx,up, barx, low, angle=90, code=3, length=.1)
abline(h=0.032, lwd=1.5)
abline(h=0.023, lty=2)
abline(h=0.042, lty=2)
dev.off()




#########################################################
### Neutral markers

geo <- read.csv("xy.csv")
mydata <- df2genind(MarkerInfo[,2:12], pop=MarkerInfo$pop)
summary(mydata)
mydata@other$xy <- geo

mydata.na <- na.replace(mydata, met=0)
summary(mydata.na)
mydata.na

##################Data Summary#############
head(mydata)
toto <- summary(mydata)

par(mfrow=c(2,2))
plot(toto$pop.eff, toto$pop.nall, xlab="Pops sample size", ylab="Number of alleles", main="Alleles numbers and Sample Size")
text(toto$pop.eff,toto$pop.nall,lab=names(toto$pop.eff))
barplot(toto$loc.nall,ylab="Number of alleles", main="number alleles per locus")
barplot(toto$Hexp-toto$Hobs, main= "hetero, exp-obsv", ylab="Hexp-Hobs")
barplot(toto$pop.eff,main="sample sizes per pop", ylab="number of genotypes",las=3)

#### allele stats
# function for standard error
std <- function(x) sd(x)/sqrt(length(x))
# alleles per locus
mean(toto$loc.nall)
std(toto$loc.nall)
# observed heterozygosity
mean(toto$Hobs)
std(toto$Hobs)


#### hwe
toto <- HWE.test.genind(mydata, res="matrix")
dim(toto)
colnames(toto)
idx <- which(toto<.0001, TRUE)
idx
toto <- HWE.test.genind(mydata, res="full")
mapply(function(i,j) toto[[i]][[j]], idx[,2], idx[,1], SIMPLIFY=FALSE)
toto
################

pca1 <- dudi.pca(mydata.na, scannf=FALSE,scale=FALSE)
temp<- as.integer(pop(mydata.na))
plot(pca1$li, col=temp, pch=19)
# how much to the PCs explain
eig <- data.frame(pca1$eig)
eig$percentage1 <- (eig[, 1]/sum(eig$pca1.eig))
sum(eig$percentage1)
sum(eig$percentage1[1:2])
eig$percentage2 <- (eig[, 2]/sum(eig$pca1.eig))
sum(eig$percentage2)
sum(eig$percentage2[1:2])

fst.col <- c("#404520",  "#BF5BD5",	"#68DA48",	"#DA485A",	"#77A8D0",	"#D1A540",	
             "#623C72",	"#CCC79E",	"#649036",	"#C94B92",	"#C07C7E",	"#6C75CC",	
             "#7ED992",	"#C9DA4E",	"#7AD2CA",	"#506372",	"#D45728",	"#7C2B2C",	
             "#412833",	"#5F8565",	"#D0A9CB",	"#A36E3D")
pops <- c("17.5E.Bailey",  "Bailey",	"Blackhawk",	"Copperdale Inn",	"Cuchara",
          "Elbert",	"Foxton",	"Gross Dam Rd",	"La Veta",	"Mosca",	"North Fork",
          "Pine Valley rd",	"Rt. 72",	"Rollinsville Rd", "Spring Creek",
          "San Francisco",	"Sugar Loaf",	"Sawmill Rd",	"Scaffer'sCrossing",
          "N.LaVetaPass",	"Nathrop",	"Wilkerson Pass")
pops.code <- c("EBY",  "BLY",  "BLH",	"COI",	"CUC", "ELB",	"FOX",	"GDR",	"LAV",	"MOS",	"NFK",
          "PVR",	"R72",	"RVR", "SPC","SFR",	"LFC",	"MRS",	"SCR","NLV",	"WNR",	"WKP")
pdf("fst_pca_full.pdf")
par(mar=c(5,5,4,2), oma=c(0,0,0,0))
test1 <- pca1$li
plot(pca1$li$Axis1, pca1$li$Axis2,col=fst.col, pch="", cex=2,xlim=c(-2,2), ylim=c(-2, 2),
     xlab="Principle Component 1", ylab="Principle Component 2", 
     main="Microsatellite PCA")
dataEllipse(pca1$li$Axis1, pca1$li$Axis2, groups=as.factor(temp), levels=.95, 
            plot.points=FALSE, col=fst.col, group.labels=pops.code)
dev.off()







##############################
# Matching PCA

# Quantitative PCA
x <- box.scale.na[4:9]
pca.x <- princomp(x,scale=TRUE, center=TRUE)
quant <- cbind(box.scale.na$pop,pca.x$scores[,1], pca.x$scores[,2])
quant.pop <- levels(as.factor(box.scale.na$pop))
quantx <- tapply(quant[,2], quant[,1], mean)
quanty <- tapply(quant[,3], quant[,1], mean)

# Microsat PCA
pca1 <- dudi.pca(mydata.na, scannf=FALSE,scale=FALSE)
temp<- as.integer(pop(mydata.na))
micro <- cbind(temp, pca1$li$Axis1, pca1$li$Axis2)
micro.pop <- levels(as.factor(pop(mydata.na)))
microx <- tapply(micro[,2], micro[,1], mean)
microy <- tapply(micro[,3], micro[,1], mean)

pdf("bothpca_bw.pdf", height=5, width=10)
par(mfrow = c(1,2)) 
par(mar=c(5,5,4,2), oma=c(0,0,0,0))

plot(microx, microy, xlim=c(-2,2), ylim=c(-2, 2), pch=19,xlab="Principle Component 1", 
     ylab="Principle Component 2", main="Microsatellite PCA")
dataEllipse(pca1$li$Axis1, pca1$li$Axis2, groups=as.factor(temp), levels=.95, 
            plot.points=FALSE, col=rep("grey68",length(micro.pop)), 
            group.label=rep("",length(micro.pop)),center.pch=FALSE)
points(microx, microy,pch=19)
plot(quantx, quanty, xlim=c(-3, 6), ylim=c(-5.5, 5.5), pch=19, xlab="Principle Component 1",
     ylab="Principle Component 2", main="Phenotype PCA for Floral Traits")
dataEllipse(pca.x$scores[,1], pca.x$scores[,2], groups=as.factor(box.scale.na$pop), 
            levels=.95, plot.points=FALSE, center.pch=FALSE, 
            col=rep("grey68",length(quant.pop)), group.label=rep("",length(quant.pop)))
points(quantx, quanty,pch=19)
dev.off()



##############################
#### isolation by distance
par(mfrow=c(1,1))
mydata2 <- genind2genpop(mydata)
# Dgen <- dist.genpop(mydata2, method=3)

# f-statistics

p.fst <- pairwise.fst(mydata)
p.fst
fstat(mydata)
fsttab <- Fst(as.loci(mydata))
fsttab
apply(fsttab, 2, mean)
apply(fsttab, 2, sd)


Dgen <- quasieuclid(pairwise.fst(mydata, res.type="dist"))
Dgeo <-dist(mydata2@other$xy)
ibd <- mantel.randtest(Dgen, Dgeo, nrepet=9999)
ibd
plot(ibd)
plot(Dgeo, Dgen)
abline(lm(Dgen~Dgeo))
dens <- kde2d(Dgeo, Dgen, n=100, lims=c(-.1, 3.1, -.01, .1))
library(grDevices)
graypal <- colorRampPalette(gray.colors(10, start=.1, end=.9))
mypal <- colorRampPalette(c("white", "#E0BD7F", "#C2E195", "#A5D9CA", "#DCBFCB", "#D5CE71"))
pdf("ibyd_bw.pdf")
plot(Dgeo, Dgen, pch=20, cex=.5, xlab="Geographic Distance (Decimal Degrees)", ylab="Genetic Distance (Pairwise Fst)")
#image(dens, col=transp(mypal(100), .7), add=TRUE)
image(dens, col=transp(greypal(100), .8), add=TRUE)
points(Dgeo, Dgen, pch=20, cex=.5)
abline(lm(Dgen~Dgeo))
title(main="Isolation by Distance")
dev.off()







###################
# save.image("inprogress_april7.RData")

# Figureing out corolla length and width problems
cl.cw.test <- cbind(TraitInfo[,c(3:4,11:12)],box.scale[,4:5])
pairs.panels(cl.cw.test, cex=1)
pairs.panels(box.scale.na[,4:5], cex=1)

test1 <- cbind(TraitInfo[,c(3:6)],box.scale[,4:7])
pairs.panels(test1, cex=1)

