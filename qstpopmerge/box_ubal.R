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

setwd("~/Dropbox/QSTFST/qstpopmerge/")

MarkerInfo <- read.csv("geno.csv")
TraitInfo <- read.csv("qstmerge_unbal.csv")

# 
# ## normalize differently than natural log
# ## quantile normalization function from John's github
# quantnorm<-function(x) {
#   n=sum(!is.na(x),na.rm=T)
#   x=rank(x)/(n+1)
#   x=qnorm(x)
#   x[is.infinite(x)]=NA
#   x
# }
# 
# # Testing difference in scaling
# test<-quantnorm(TraitInfo[,3])
# test2<-scale(TraitInfo[,3])
# hist(TraitInfo[,3])
# hist(test)
# hist(test2)
# 
# ## Apply quantnorm function
# quant.data <- sapply(TraitInfo[,3:10],quantnorm)
# head(quant.data)
# quant.data<-cbind(TraitInfo[,1:2],quant.data)
# head(quant.data)
# quant.data<-quant.data[,1:10]
# names(quant.data)<- c("pop", "dam","qn_Corollawidth","qn_Corollalength","qn_Petallength",  
#                       "qn_color492","qn_nectaramt",	"qn_nectarconc","qn_Root", "qn_fdate")
# 
# ## Rename original trait info, if necessary
# # head(TraitInfo)
# # names(TraitInfo) <- c("pop", "dam", "fdate",  "Corollawidth",  "Corollalength",	"Petallength",	"nectaramt",	
# #                       "nectarconc",	"color492",	"Root",	"ln_Corollawidth",	"ln_Corollalength",	"ln_Petallength",	
# #                       "ln_nectaramt",	"ln_nectarconc",	"ln_color_adj",	"ln_Root",	"ln_fdate")
# 
# ## Merge quant norm traits with original file
# all.data<-cbind(TraitInfo[,1:18], quant.data[,3:10])
# head(all.data)
# 
# ## Look at trait distributions
# pdf("trait_dist.pdf", height=5, width=10)
# for (i in c(3:length(all.data))){
#   hist(all.data[,i], main=paste(names(all.data)[[i]]))
# }
# dev.off()
# # 
# ## QQ visualization plots
# par(mfcol=c(2,4))
# 
# # raw
# qqnorm(all.data[,3]);qqline(all.data[,3])
# qqnorm(all.data[,4]);qqline(all.data[,4])
# qqnorm(all.data[,5]);qqline(all.data[,5])
# qqnorm(all.data[,6]);qqline(all.data[,6])
# qqnorm(all.data[,7]);qqline(all.data[,7])
# qqnorm(all.data[,8]);qqline(all.data[,8])
# qqnorm(all.data[,9]);qqline(all.data[,9])
# qqnorm(all.data[,10]);qqline(all.data[,10])
# # ln transformed
# qqnorm(all.data[,11]);qqline(all.data[,11])
# qqnorm(all.data[,12]);qqline(all.data[,12])
# qqnorm(all.data[,13]);qqline(all.data[,13])
# qqnorm(all.data[,14]);qqline(all.data[,14])
# qqnorm(all.data[,15]);qqline(all.data[,15])
# qqnorm(all.data[,16]);qqline(all.data[,16])
# qqnorm(all.data[,17]);qqline(all.data[,17])
# qqnorm(all.data[,18]);qqline(all.data[,18])
# # quantile norm.
# qqnorm(all.data[,19]);qqline(all.data[,19])
# qqnorm(all.data[,20]);qqline(all.data[,20])
# qqnorm(all.data[,21]);qqline(all.data[,21])
# qqnorm(all.data[,22]);qqline(all.data[,22])
# qqnorm(all.data[,23]);qqline(all.data[,23])
# qqnorm(all.data[,24]);qqline(all.data[,24])
# qqnorm(all.data[,25]);qqline(all.data[,25])
# qqnorm(all.data[,26]);qqline(all.data[,26])
# 
# ## scale and center around the mean
# 
# test <-scale(all.data[,3:26])
# par(mfcol=c(2,4))
# # raw
# qqnorm(test[,3]);qqline(test[,3])
# qqnorm(test[,4]);qqline(test[,4])
# qqnorm(test[,5]);qqline(test[,5])
# qqnorm(test[,6]);qqline(test[,6])
# qqnorm(test[,7]);qqline(test[,7])
# qqnorm(test[,8]);qqline(test[,8])
# qqnorm(test[,9]);qqline(test[,9])
# qqnorm(test[,10]);qqline(test[,10])
# # ln transformed
# qqnorm(test[,11]);qqline(test[,11])
# qqnorm(test[,12]);qqline(test[,12])
# qqnorm(test[,13]);qqline(test[,13])
# qqnorm(test[,14]);qqline(test[,14])
# qqnorm(test[,15]);qqline(test[,15])
# qqnorm(test[,16]);qqline(test[,16])
# qqnorm(test[,17]);qqline(test[,17])
# qqnorm(test[,18]);qqline(test[,18])
# # quantile norm.
# qqnorm(test[,19]);qqline(test[,19])
# qqnorm(test[,20]);qqline(test[,20])
# qqnorm(test[,21]);qqline(test[,21])
# qqnorm(test[,22]);qqline(test[,22])
# qqnorm(test[,23]);qqline(test[,23])
# qqnorm(test[,24]);qqline(test[,24])
# qqnorm(test[,25]);qqline(test[,25])
# qqnorm(test[,26]);qqline(test[,26])
# 
# par(mfcol=c(2,4))
# # finalset
# qqnorm(finalset[,3]);qqline(finalset[,3])
# qqnorm(finalset[,4]);qqline(finalset[,4])
# qqnorm(finalset[,5]);qqline(finalset[,5])
# qqnorm(finalset[,6]);qqline(finalset[,6])
# qqnorm(finalset[,7]);qqline(finalset[,7])
# qqnorm(finalset[,8]);qqline(finalset[,8])
# qqnorm(finalset[,9]);qqline(finalset[,9])
# qqnorm(finalset[,10]);qqline(finalset[,10])
# 
# ## Dedcide which transformations to use for each trait
# ## Keep all sets of traits in the following order
# ## Corollawidth,  Corollalength, Petallength, color, nectaramt, nectarconc, Root, fdate,  
# 
# raw <- all.data[ c(1:10)]
# raw.na <- na.exclude(raw)
# natlog <- all.data[ c(1,2,11:18)]
# quantnorm <- all.data[ c(1,2,19:26)]
# 
# 
# finalset <- test[ c(1,2,11,4,13:17,26)]
# ## Exclude NAs for scatterplot matrix and PCA
# ## but use a non-excluded dataset for the Qst analysis
# finalset.na <- na.exclude(finalset)
# 
# ## save files
# write.csv(finalset, "finalset.csv")
# write.csv(finalset.na, "finalset_na.csv")

#########################################################
#### BoxCox transformations

box <- TraitInfo[1:10]
head(box)

# tmp <- boxcox(box$Corollawidth ~ box$dam)
# lam <- tmp$x[max(tmp$y)==tmp$y]
# new <- box$Corollawidth^lam

par(mfrow=c(2,4))
funbox <- function(trait, dam){
  tmp <- boxcox(trait ~ dam)
  lam <- tmp$x[max(tmp$y)==tmp$y]
  return(lam)
}
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
#### Scatterplot matrix, correlations and PCA

### Scatterplot matrix
# pairs.panels(raw.na[3:10], pch=20, cex = 1)
pairs.panels(box.scale.na[4:11], pch=20, cex = 1)

# ### Sort by highest correlation
# # function from http://little-book-of-r-for-multivariate-analysis.readthedocs.org/en/latest/src/multivariateanalysis.html
# mosthighlycorrelated <- function(mydataframe,numtoreport){
#   cormatrix <- cor(mydataframe)
#   diag(cormatrix) <- 0
#   cormatrix[lower.tri(cormatrix)] <- 0
#   fm <- as.data.frame(as.table(cormatrix))
#   names(fm) <- c("First.Variable", "Second.Variable","Correlation")
#   head(fm[order(abs(fm$Correlation),decreasing=T),],n=numtoreport)
# }
# high_cor <- mosthighlycorrelated(box.scale.na[,4:11], 20)
# high_cor
# # write.table(high_cor, file="top_correlations.csv")


### Basic PCA

# means and standard deviations
sapply(box.scale.na[4:11],mean) 
sapply(box.scale.na[4:11],sd) 

# Perforn the PCA on scaled data
x <- box.scale.na[4:11]
pca.x <- princomp(x,scale=TRUE, center=TRUE)
summary(pca.x)
names(pca.x)
# unclass(pca.x) # prints all pca info

screeplot(pca.x, type="barplot")

# Trait loadings on each PC
pca.x$loadings

### Barplots for trait loadings on the first 4 PCs
xname <- c("CW", "CL", "PL", "COL", "NA", "NC", "RT", "FD")
color <- c("darkred", "darkred", "darkred", "darkred", "yellow", "yellow", "darkgreen", "darkgreen")
summary(pca.x)
pdf("pheno_pcs.pdf")
par(mfrow=c(2,2), mar=c(3,3,3,1), oma=c(0,0,3,0))
barplot(pca.x$loadings[,1], names.arg=xname, col=color, main="PC1, 36.5%")
box(which="figure")
barplot(pca.x$loadings[,2], names.arg=xname, col=color, main="PC2, 16.1%")
box(which="figure")
barplot(pca.x$loadings[,3], names.arg=xname, col=color, main="PC3, 13.6%")
box(which="figure")
barplot(pca.x$loadings[,4], names.arg=xname, col=color, main="PC4, 11.6%")
box(which="figure")
title(main="Trait Loadings for Phenotype PCA", cex.main=2, outer=TRUE)
dev.off()
par(mfrow=c(1,1))


# Plot the PCA
pdf("pheno_pca_box.pdf")
# pca.x$scores contains the actual points for plotting
plot(pca.x$scores[,1], pca.x$scores[,2], pch=".", cex=2, xlim=c(-3, 7), ylim=c(-5, 5),
     xlab="Principle Component 1", ylab="Principle Component 2", 
     main="Phenotype PCA: All traits\nMerged Populations 95% CI") 
# Color points by pop
# plot(pca.x$scores[,1], pca.x$scores[,2], col=as.numeric(box.scale.na$Pop), pch=19)
#legend("topright", legend=unique(box.scale.na$Pop), horiz=TRUE, 
#col=unique(as.numeric(box.scale.na$Pop)), pch=19)
# add ellipses and a grid
pops <- c("westchicagocreek","guanellapass","Golden","Smiths","Cuchara","noe",
          "hwy86","72med","72high","caribou_rd","springcreek","Lefthand","MRS",
          "wilkersonpass","TopCaribou","hwy105","Georgetown")
newcol <- c("#C7D34F","#CA53D5","#D65033","#89C9CA","#56335B","#587C3D","#CB98C0",
            "#C6AB8E","#7B4627","#796ECD","#83D598","#CF9A3F","#C84792","#3F4239",
            "#6BD249","#B84B5C","#6286AE")
dataEllipse(pca.x$scores[,1], pca.x$scores[,2], groups=as.factor(box.scale.na$pop), levels=.95, 
            plot.points=FALSE, col=newcol, group.labels=pops)
dev.off()

#biplot of first two principal components
biplot(pca.x)

### Add PCs 1-4 to final dataset
pca.1.4 <- data.frame(pca.x$scores[,c(1:4)])
pca.1.4$id<-rownames(pca.1.4)
box.scale.pca <- merge(box.scale, pca.1.4, by="id",all=T)
write.csv(box.scale.pca, "box.scale.pca.csv")
box.scale.pca <- read.csv("box.scale.pca.csv")
box.scale.pca <- box.scale.pca[,2:16]
head(box.scale.pca)
# 
# ######
# # testing absolute value pcs
# hist(pca.x$scores[,1])
# hist(abs(pca.x$scores[,1]))
# hist(pca.x$scores[,2])
# hist(abs(pca.x$scores[,2]))
# pc1 <- (pca.x$scores[,1])
# pc1abs <- (abs(pca.x$scores[,1]))
# pc2 <- (pca.x$scores[,2])
# pc2abs <- (abs(pca.x$scores[,2]))
# 
# pca.1.2abs <- data.frame(cbind(pc1, pc1abs, pc2, pc2abs))
# pca.1.2abs$id<-rownames(pca.1.2abs)
# box.scale.pca.abs <- merge(box.scale, pca.1.2abs, by="id",all=T)
# write.csv(box.scale.pca.abs, "box.scale.pca.abs.csv")
# box.scale.pca.abs <- read.csv("box.scale.pca.abs.csv")
# box.scale.pca.abs <- box.scale.pca.abs[,2:16]
# head(box.scale.pca.abs)


#########################################################
###########     USE THIS     ############
#########################################################
#### Scatterplot matrix, correlations and PCA for ***floral traits only***

### Scatterplot matrix
# pairs.panels(raw.na[3:10], pch=20, cex = 1)
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
# unclass(pca.x) # prints all pca info

screeplot(pca.x, type="barplot")

# Trait loadings on each PC
pca.x$loadings

### Barplots for trait loadings on the first 4 PCs
xname <- c("CW", "CL", "PL", "COL", "NA", "NC")
color <- c("darkred", "darkred", "darkred", "darkred", "yellow", "yellow")
summary(pca.x)
pdf("pheno_pcs-flower.pdf")
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
dev.off()
par(mfrow=c(1,1))


# Plot the PCA
pdf("pheno_flower_pca_box.pdf")
# pca.x$scores contains the actual points for plotting
plot(pca.x$scores[,1], pca.x$scores[,2], pch="", cex=2, xlim=c(-3, 6), ylim=c(-5.5, 5.5),
     xlab="Principle Component 1", ylab="Principle Component 2", 
     main="Phenotype PCA for Floral Traits") 
# Color points by pop
# plot(pca.x$scores[,1], pca.x$scores[,2], col=as.numeric(box.scale.na$Pop), pch=19)
#legend("topright", legend=unique(box.scale.na$Pop), horiz=TRUE, 
#col=unique(as.numeric(box.scale.na$Pop)), pch=19)
# add ellipses and a grid
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
dev.off()

#biplot of first two principal components
biplot(pca.x)

### Add PCs 1-4 to final dataset
pca.1.4 <- data.frame(pca.x$scores[,c(1:4)])
pca.1.4$id<-rownames(pca.1.4)
box.scale.pca <- merge(box.scale, pca.1.4, by="id",all=T)
write.csv(box.scale.pca, "box.scale.pca.csv")
box.scale.pca <- read.csv("box.scale.pca.csv")
box.scale.pca <- box.scale.pca[,2:16]
head(box.scale.pca)

# ######
# # testing absolute value pcs
# hist(pca.x$scores[,1])
# hist(abs(pca.x$scores[,1]))
# hist(pca.x$scores[,2])
# hist(abs(pca.x$scores[,2]))
# pc1 <- (pca.x$scores[,1])
# pc1abs <- (abs(pca.x$scores[,1]))
# pc2 <- (pca.x$scores[,2])
# pc2abs <- (abs(pca.x$scores[,2]))
# 
# pca.1.2abs <- data.frame(cbind(pc1, pc1abs, pc2, pc2abs))
# pca.1.2abs$id<-rownames(pca.1.2abs)
# box.scale.pca.abs <- merge(box.scale, pca.1.2abs, by="id",all=T)
# write.csv(box.scale.pca.abs, "box.scale.pca.abs.csv")
# box.scale.pca.abs <- read.csv("box.scale.pca.abs.csv")
# box.scale.pca.abs <- box.scale.pca.abs[,2:16]
# head(box.scale.pca.abs)
# 



#########################################################
#### Qst-Fst

# # This script only runs one phenotype at a time.
# # Need to create objects that pair each phenotype with "pop" and "dam" columns
# trait.col <- list()
# for (i in c(4:length(box.scale.pca))){
#   trait.col[[i]] <- box.scale.pca[,c(2,3,i)]
# }
# 
# # Name list items
# # trait.col
# names(trait.col) <- c("ID", "Pop", "Dam", "CorollaWidth",  "CorollaLength", "PetalLength", "Color", 
#                       "NectarAmt",	"NectarConc", "Root", "Fdate",  "PC1",  "PC2",  "PC3","PC4")


# This script only runs one phenotype at a time.
# Need to create objects that pair each phenotype with "pop" and "dam" columns
trait.col <- list()
for (i in c(4:length(box.scale.pca))){
  trait.col[[i]] <- box.scale.pca[,c(2,3,i)]
}

# Name list items
# trait.col
names(trait.col) <- c("ID", "Pop", "Dam", "CorollaWidth",  "CorollaLength", "PetalLength", "Color", 
                      "NectarAmt",  "NectarConc", "Root", "Fdate",  "PC1",  "PC1abs",  "PC2","PCabs")


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
write.csv(results, "results.abs.csv")



#### Plotting Qst

res.name <- c("Corolla Width",  "Corolla Length", "Petal Length", "Color", "Nectar Amount",  
           "Nectar Conc.", "Root Diameter", "Flowering Date",  "PC1",  "PC1abs",  "PC2","PCabs")
col <- c("darkred", "darkred", "darkred", "darkred", "yellow", "yellow", "darkgreen", 
         "darkgreen", "grey", "grey", "grey", "grey")
pdf("qst_bar_withpca.abs.pdf")
par(mfrow=c(1,1), mar=c(9,4,1,2))
barplot(results[,4], names.arg=res.name,las=2, cex.names=1, col=col, main="test")
# error bars
# could not get the confidence intervals out of list format so I created a work around
barx <- barplot(results[,4], names.arg=res.name,las=2, cex.names=1, col=col, main="Qst with 95% CI")
low <- unlist(qst.low, use.names=FALSE)
up<- unlist(qst.up, use.names=FALSE)
arrows(barx,up, barx, low, angle=90, code=3, length=.1)
abline(h=0.03201393)
abline(h=0.023, lty=2)
abline(h=0.042, lty=2)
dev.off()

############################
#### Plotting Qst: no PCs

results.all <- read.csv("results.csv")
results <- results.all[1:8,2:7]
rownames(results) <- results.all$X[1:8]

res.name <- c("Corolla Width",  "Corolla Length", "Petal Length", "Color", "Nectar Volume",  
              "Nectar Conc.", "Root Diameter", "Flowering Date")
col <- c("darkred", "darkred", "darkred", "darkred", "yellow", "yellow", "darkolivegreen3", 
         "darkolivegreen3")
col <- c("darkred", "darkred", "darkred", "darkred", "yellow", "yellow", "darkolivegreen3", 
         "darkolivegreen3")
pdf("qst_bar_nopc.pdf")
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
##############################
# Matching PCA
fpops.code <- c("EBY",  "BLY",  "BLH",  "COI",  "CUC", "ELB",	"FOX",	"GDR",	"LAV",	"MOS",	"NFK",
                "PVR",	"R72",	"RVR", "SPC","SFR",	"LFC",	"MRS",	"SCR","NLV",	"WNR",	"WKP")
fst.col <- c("#404520",  "#BF5BD5",  "#68DA48",	"#DA485A",	"#77A8D0",	"#D1A540",	
             "#623C72",	"#CCC79E",	"#649036",	"#C94B92",	"#C07C7E",	"#6C75CC",	
             "#7ED992",	"#C9DA4E",	"#7AD2CA",	"#506372",	"#D45728",	"#7C2B2C",	
             "#412833",	"#5F8565",	"#D0A9CB",	"#A36E3D")


# TV code
alpha <- "ff"
grey <- rep("#d3d3d3", length(fst.col))
greyq <- rep("#d3d3d3", length(qst.col))
#cex.thor <- c(2,1,0.7)
#cex.main = cex.thor[1]
#col = paste(grey[1], alpha, sep = "")


x <- box.scale.na[4:9]
pca.x <- princomp(x,scale=TRUE, center=TRUE)
qpops.code <- c("WCC","BLY","BLH","COI","CUC","NOE",
                "H86","GDR","R72","CAR","SPC","LFC","MRS",
                "WKP","RVR","H105","GTN")
qst.col <- c("#C7D34F","#CA53D5","#D65033","#89C9CA","#56335B","#587C3D","#CB98C0",
            "#C6AB8E","#7B4627","#796ECD","#83D598","#CF9A3F","#C84792","#3F4239",
            "#6BD249","#B84B5C","#6286AE")
pca1 <- dudi.pca(mydata.na, scannf=FALSE,scale=FALSE)
temp<- as.integer(pop(mydata.na))
#pdf("bothpca.pdf", height=5, width=10)
par(mfrow = c(1,2)) 
par(mar=c(5,5,4,2), oma=c(0,0,0,0))

plot(pca1$li$Axis1, pca1$li$Axis2,col=fst.col, pch="",xlim=c(-2,2), ylim=c(-2, 2),
     xlab="Principle Component 1", ylab="Principle Component 2", 
     main="A) Microsatellite PCA")
dataEllipse(pca1$li$Axis1, pca1$li$Axis2, groups=as.factor(temp), levels=.95, 
            plot.points=FALSE, col=paste(grey,alpha,sep="")) #group.labels=fpops.code

plot(pca.x$scores[,1], pca.x$scores[,2], pch="", cex=2, xlim=c(-3, 6), ylim=c(-5.5, 5.5),
     xlab="Principle Component 1", ylab="Principle Component 2", 
     main="B) Phenotype PCA for Floral Traits") 

dataEllipse(pca.x$scores[,1], pca.x$scores[,2], groups=as.factor(box.scale.na$Pop), levels=.95, 
            plot.points=FALSE, col=qst.col, group.labels=qpops.code)
# dev.off()

micro <- cbind(temp, pca1$li$Axis1, pca1$li$Axis2)
quant <- cbind(box.scale.na$pop,pca.x$scores[,1], pca.x$scores[,2])

quantx <- tapply(quant[,2], quant[,1], mean)
quanty <- tapply(quant[,3], quant[,1], mean)
microx <- tapply(micro[,2], micro[,1], mean)
microy <- tapply(micro[,3], micro[,1], mean)

pdf("bothpca_bw.pdf", height=5, width=10)
par(mfrow = c(1,2)) 
par(mar=c(5,5,4,2), oma=c(0,0,0,0))

plot(microx, microy, xlim=c(-2,2), ylim=c(-2, 2), pch=19,xlab="Principle Component 1", 
     ylab="Principle Component 2", main="Microsatellite PCA")
dataEllipse(pca1$li$Axis1, pca1$li$Axis2, groups=as.factor(temp), levels=.95, 
            plot.points=FALSE, col=rep("grey68",length(fst.col)), 
            group.label=rep("",length(fst.col)),center.pch=FALSE)
points(microx, microy,pch=19)
plot(quantx, quanty, xlim=c(-3, 6), ylim=c(-5.5, 5.5), pch=19, xlab="Principle Component 1",
     ylab="Principle Component 2", main="Phenotype PCA for Floral Traits")
dataEllipse(pca.x$scores[,1], pca.x$scores[,2], groups=as.factor(box.scale.na$pop), 
            levels=.95, plot.points=FALSE, center.pch=FALSE, 
            col=rep("grey68",length(qst.col)), group.label=rep("",length(qst.col)))
points(quantx, quanty,pch=19)
dev.off()

##############################
##############################
#isolation by distance
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
