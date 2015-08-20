## Analysis of matched populations

library(devtools)
library(QstFstComp)
library(psych) # for pairs.panels
library(car) # for dataEllipse
library(MASS) # for boxcox transformation
library(RColorBrewer)
library(adegenet)
library(KernSmooth) # for package neutrality test.r
library(corpcor) # for package neutrality test.r
library(hierfstat)

setwd("~/Dropbox/QSTFST/match/")

MarkerInfo <- read.csv("matched_geno.csv")
TraitInfo <- read.csv("matchpop_unbal.csv")
names(TraitInfo)


## trait distributions

## Look at trait distributions
pdf("match_trait_dist.pdf", height=5, width=10)
for (i in c(3:length(TraitInfo))){
  hist(TraitInfo[,i], main=paste(names(TraitInfo)[[i]]))
}
dev.off()

## QQ visualization plots
par(mfcol=c(2,4))
# raw
qqnorm(TraitInfo[,3]);qqline(TraitInfo[,3])
qqnorm(TraitInfo[,4]);qqline(TraitInfo[,4])
qqnorm(TraitInfo[,5]);qqline(TraitInfo[,5])
qqnorm(TraitInfo[,6]);qqline(TraitInfo[,6])
qqnorm(TraitInfo[,7]);qqline(TraitInfo[,7])
qqnorm(TraitInfo[,8]);qqline(TraitInfo[,8])
qqnorm(TraitInfo[,9]);qqline(TraitInfo[,9])
qqnorm(TraitInfo[,10]);qqline(TraitInfo[,10])

test <-scale(TraitInfo[,3:10])
par(mfcol=c(2,4))
# raw
qqnorm(test[,1]);qqline(test[,1])
qqnorm(test[,2]);qqline(test[,2])
qqnorm(test[,3]);qqline(test[,3])
qqnorm(test[,4]);qqline(test[,4])
qqnorm(test[,5]);qqline(test[,5])
qqnorm(test[,6]);qqline(test[,6])
qqnorm(test[,7]);qqline(test[,7])
qqnorm(test[,8]);qqline(test[,8])



#### BoxCox transformations

box <- TraitInfo
head(box)

funbox <- function(trait, dam){
  tmp <- boxcox(trait ~ dam)
  lam <- tmp$x[max(tmp$y)==tmp$y]
  return(lam)
}
par(mfcol=c(2,4))
lamda <- as.vector(apply(box[3:10], 2, funbox, dam=box$Dam))

head(box)

b1 <- box$CorollaWidth^lamda[1]
b2 <- box$CorollaLength^lamda[2]
b3 <- box$PetalLength^lamda[3]
b4 <- box$Color^lamda[4]
b5 <- box$NectarAmt^lamda[5]
b6 <- box$NectarConc^lamda[6]
b7 <- box$Root^lamda[7]
b8 <- box$Fdate^lamda[8]

all.box <- cbind(box[,1:2], b1, b2, b3, b4, b5, b6, b7, b8)

par(mfrow=c(2,4))
apply(all.box[3:10], 2, hist)

all.box1 <- scale(all.box[3:10])
box.scale <- cbind(all.box[1:2], all.box1)
names(box.scale) <- names(box)

pdf("box_trait_dist.pdf", height=5, width=10)
par(mfrow=c(2,4))
apply(box.scale[3:10], 2, hist)
apply(box.scale[3:10], 2, qqnorm)
dev.off()

box.scale.na <- na.exclude(box.scale)

## save files
write.csv(box.scale, "box.scale.csv")
write.csv(box.scale.na, "box.scale.na.csv")
box.scale <- read.csv("box.scale.csv")
box.scale.na <- read.csv("box.scale.na.csv")
names(box.scale) <- c("ID","Pop","Dam","CorollaWidth","CorolLalength","PetalLength",
                      "Color","NectarAmt","NectarConc","Root","FlowerDate")
names(box.scale.na) <- c("ID","Pop","Dam","CorollaWidth","CorolLalength","PetalLength",
                         "Color","NectarAmt","NectarConc","Root","FlowerDate")


#########################################################
#### Scatterplot matrix, correlations and PCA

### Scatterplot matrix
#pairs.panels(raw.na[3:10], pch=20, cex = 1)
pdf("match_box_corr.pdf", height=5, width=8)
pairs.panels(box.scale.na[4:11], pch=20, cex = 1)
dev.off()

### Sort by highest correlation
# function from http://little-book-of-r-for-multivariate-analysis.readthedocs.org
                                  # /en/latest/src/multivariateanalysis.html
mosthighlycorrelated <- function(mydataframe,numtoreport){
  cormatrix <- cor(mydataframe)
  diag(cormatrix) <- 0
  cormatrix[lower.tri(cormatrix)] <- 0
  fm <- as.data.frame(as.table(cormatrix))
  names(fm) <- c("First.Variable", "Second.Variable","Correlation")
  head(fm[order(abs(fm$Correlation),decreasing=T),],n=numtoreport)
}
high_cor <- mosthighlycorrelated(box.scale.na[,4:11], 20)
high_cor
# write.table(high_cor, file="forage_top_correlations.csv")


### Basic PCA

# means and standard deviations
sapply(box.scale.na[4:11],mean) 
sapply(box.scale.na[4:11],sd) 

# Perforn the PCA on scaled data
x <- box.scale.na[4:11]
pca.x <- princomp(x,scale=TRUE, center=TRUE)
summary(pca.x)

screeplot(pca.x, type="barplot")

pca.x$loadings

### Barplots for trait loadings on the first 4 PCs
xname <- c("CW", "CL", "PL", "COL", "NA", "NC", "RT", "FD")
color <- c("darkred", "darkred", "darkred", "darkred", "yellow", "yellow", "darkgreen", "darkgreen")
summary(pca.x)
pdf("match_box_pheno_pcs.pdf")
par(mfrow=c(2,2), mar=c(3,3,3,1), oma=c(0,0,3,0))
barplot(pca.x$loadings[,1], names.arg=xname, col=color, main="PC1, 33.5%")
box(which="figure")
barplot(pca.x$loadings[,2], names.arg=xname, col=color, main="PC2, 16.8%")
box(which="figure")
barplot(pca.x$loadings[,3], names.arg=xname, col=color, main="PC3, 13.0%")
box(which="figure")
barplot(pca.x$loadings[,4], names.arg=xname, col=color, main="PC4, 11.9%")
box(which="figure")
title(main="Trait Loadings for Phenotype PCA", cex.main=2, outer=TRUE)
dev.off()
par(mfrow=c(1,1))

newcol <- c("#269C49","#8061F0","#CD2B42","#F68E08","#71500F",
            "#9BB710","#841F57","#19542D","#E569E0","#DA9E6C","#658419")

# Plot the PCA
pdf("pheno_pca_matched_box.pdf")
# pca.x$scores contains the actual points for plotting
plot(pca.x$scores[,1], pca.x$scores[,2], pch=".", cex=2, xlim=c(-3, 7), ylim=c(-5, 5),
     xlab="Principle Component 1", ylab="Principle Component 2", 
     main="Phenotype PCA: All traits\nMatched Populations 95% CI") 
# Color points by pop
# plot(pca.x$scores[,1], pca.x$scores[,2], col=as.numeric(box.scale.na$Pop), pch=19)
#legend("topright", legend=unique(box.scale.na$Pop), horiz=TRUE, 
    #col=unique(as.numeric(box.scale.na$Pop)), pch=19)
# add ellipses and a grid
pops <- c("MRS","Lefthand","72med","72high","Smiths","Caribou","Golden",
          "GualellaPass","SpringCreek","WilkersonPass","Cuchara")
dataEllipse(pca.x$scores[,1], pca.x$scores[,2], groups=as.factor(box.scale.na$Pop), levels=.95, 
            plot.points=FALSE, col=newcol, group.labels=pops)
dev.off()


#biplot of first two principal components
biplot(pca.x)

### Add PCs 1-4 to final dataset
pca.1.4 <- data.frame(pca.x$scores[,c(1:4)])
pca.1.4$ID<-rownames(pca.1.4)
box.scale.pca <- merge(box.scale, pca.1.4, by="ID",all=T)



#########################################################
#### Scatterplot matrix, correlations and PCA
## Just for floral traits
### Scatterplot matrix
#pairs.panels(raw.na[3:10], pch=20, cex = 1)
pdf("match_box_corr_floral.pdf", height=5, width=8)
pairs.panels(box.scale.na[4:9], pch=20, cex = 1)
dev.off()

### Sort by highest correlation
# function from http://little-book-of-r-for-multivariate-analysis.readthedocs.org
                                      # /en/latest/src/multivariateanalysis.html
mosthighlycorrelated <- function(mydataframe,numtoreport){
  cormatrix <- cor(mydataframe)
  diag(cormatrix) <- 0
  cormatrix[lower.tri(cormatrix)] <- 0
  fm <- as.data.frame(as.table(cormatrix))
  names(fm) <- c("First.Variable", "Second.Variable","Correlation")
  head(fm[order(abs(fm$Correlation),decreasing=T),],n=numtoreport)
}
high_cor <- mosthighlycorrelated(box.scale.na[,4:9], 20)
high_cor
# write.table(high_cor, file="forage_top_correlations.csv")


### Basic PCA

# means and standard deviations
sapply(box.scale.na[4:9],mean) 
sapply(box.scale.na[4:9],sd) 

# Perforn the PCA on scaled data
x.f <- box.scale.na[4:9]
pca.x.f <- princomp(x.f,scale=TRUE, center=TRUE)
summary(pca.x.f)

screeplot(pca.x.f, type="barplot")

pca.x.f$loadings

### Barplots for trait loadings on the first 4 PCs
xfname <- c("CW", "CL", "PL", "COL", "NA", "NC")
colorf <- c("darkred", "darkred", "darkred", "darkred", "yellow", "yellow")
summary(pca.x.f)
pdf("match_box_pheno_pcs_floral.pdf")
par(mfrow=c(2,2), mar=c(3,3,3,1), oma=c(0,0,3,0))
barplot(pca.x.f$loadings[,1], names.arg=xfname, col=colorf, main="PC1, 41.4%")
box(which="figure")
barplot(pca.x.f$loadings[,2], names.arg=xfname, col=colorf, main="PC2, 17.0%")
box(which="figure")
barplot(pca.x.f$loadings[,3], names.arg=xfname, col=colorf, main="PC3, 15.6%")
box(which="figure")
barplot(pca.x.f$loadings[,4], names.arg=xfname, col=colorf, main="PC4, 11.3%")
box(which="figure")
title(main="Trait Loadings for Phenotype PCA", cex.main=2, outer=TRUE)
dev.off()

newcol <- c("#269C49","#8061F0","#CD2B42","#F68E08","#71500F",
            "#9BB710","#841F57","#19542D","#E569E0","#DA9E6C","#658419")

# Plot the PCA
pdf("pheno_pca_matched_box_floral.pdf")
# pca.x$scores contains the actual points for plotting
plot(pca.x.f$scores[,1], pca.x.f$scores[,2], pch=".", cex=2, xlim=c(-3, 7), ylim=c(-5.2, 5),
     xlab="Principle Component 1", ylab="Principle Component 2", 
     main="Phenotype PCA: Floral Traits\nMatched Populations 95% CI") 
# Color points by pop
# plot(pca.x$scores[,1], pca.x$scores[,2], col=as.numeric(box.scale.na$Pop), pch=19)
#legend("topright", legend=unique(box.scale.na$Pop), horiz=TRUE, 
  #col=unique(as.numeric(box.scale.na$Pop)), pch=19)
# add ellipses and a grid
pops <- c("MRS","Lefthand","72med","72high","Smiths","Caribou","Golden",
          "GualellaPass","SpringCreek","WilkersonPass","Cuchara")
dataEllipse(pca.x.f$scores[,1], pca.x.f$scores[,2], groups=as.factor(box.scale.na$Pop), levels=.95, 
            plot.points=FALSE, col=newcol, group.labels=pops)
dev.off()


#biplot of first two principal components
biplot(pca.x.f)


#########################################################
#### Qst-Fst
## Back to full set of traits

# This script only runs one phenotype at a time.
# Need to create objects that pair each phenotype with "pop" and "dam" columns
trait.col <- list()
for (i in c(4:length(box.scale.pca))){
  trait.col[[i]] <- box.scale.pca[,c(2,3,i)]
}

# Name list items
trait.col
names(trait.col) <- c("ID","Pop","Dam","CorollaWidth","CorolLalength","PetalLength","Color",
                      "NectarAmt","NectarConc","Root","FlowerDate", "PC1", "PC2", "PC3", "PC4")

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

save.image("april_2.RData")

#### Plotting Qst 

res.name <- c("Corolla Width",  "Corolla Length", "Petal Length", "Color", "Nectar Amount",  
              "Nectar Conc.", "Root Diameter", "Flowering Date",  "PC1",  "PC2",  "PC3","PC4")
col <- c("darkred", "darkred", "darkred", "darkred", "yellow", "yellow", "darkgreen", 
         "darkgreen", "grey", "grey", "grey", "grey")
pdf("qst_withpca_box_match.pdf")
par(mfrow=c(1,1), mar=c(8,4,3,2))
#barplot(results[,4], names.arg=res.name,las=2, cex.names=1, col=col, 
main="Qst with 95% CI for boxcox and merged trait set")
# error bars
# could not get the confidence intervals out of list format so I created a work around
barx <- barplot(results[,4], names.arg=res.name,las=2, cex.names=1, col=col, 
                main="Qst with 95% CI for boxcox and merged trait set")
low <- unlist(qst.low, use.names=FALSE)
up<- unlist(qst.up, use.names=FALSE)
arrows(barx,up, barx, low, angle=90, code=3, length=.1)
abline(h=0.03201393)
abline(h=0.023, lty=2)
abline(h=0.042, lty=2)
dev.off()


################################################################
### Neutral markers

mydata <- df2genind(MarkerInfo[,2:12], pop=MarkerInfo$pop)
summary(mydata)

mydata.na <- na.replace(mydata, met=0)
summary(mydata.na)
mydata.na

pca1 <- dudi.pca(mydata.na, scannf=FALSE,scale=FALSE)
temp<- as.integer(pop(mydata.na))
plot(pca1$li, col=temp, pch=19)

### Plot PCA with groups

pops <- c("SawmillRd","SugarLoaf","GrossDamRd","Rt72","CopperdaleInn","Rollinsville",
          "SofBlackhawk","Bailey","SpringCreek","WilkersonPass","Cuchara")
pdf("fst_pca.pdf")
plot(test1, col=newcol, pch=".", cex=2,xlim=c(-2,2), ylim=c(-3, 2),
     xlab="Principle Component 1", ylab="Principle Component 2", 
     main="Microsatellite PCA:\n95% CI for Matched populations")
dataEllipse(test1, groups=as.factor(temp), levels=.95, 
            plot.points=FALSE, col=newcol, group.labels=pops)
dev.off()



hwt <- HWE.test.genind(mydata, res="matrix")
dim(hwt)
colnames(hwt)
idx<- which(hwt<.0001, TRUE)
idx

hwt <- HWE.test.genind(mydata, res="full")
mapply(function(i,j) hwt[[i]][[j]],idx[,2], idx[,1], SIMPLIFY=FALSE)

save.image("matched_april1.RData")


#####################
# DAPC

## only found one cluster based on bayesian information criterion (BIC)
grp <- find.clusters(mydata, max.n.clusters=40)
100
2
head(grp$grp, 10)
grp$size
table(pop(mydata), grp$grp)
dapc1 <- dapc(mydata, grp$grp)
2
2
scatter(dapc1)

#####################
# G and D matricies using MANOVA code from Bessega et al. 2015
# Only floral traits

dam.mean <- aggregate(box.scale.pca[,4:9], list(box.scale.pca$Dam), mean, na.rm=T)
dam.mean
pop.dam <- unique(box.scale.pca[,c(2,3)])
trait.mean <- merge(pop.dam, dam.mean, by.y="Group.1", by.x="Dam")

write.csv(trait.mean, "trait.mean.csv")
tm <- read.csv("trait.mean.csv")
head(tm)
data <- as.matrix(tm[,4:9])
factor <- as.matrix(tm[,3])

#######testing
head(tm)
head(data)

test.data <- cbind(tm$X, tm$Pop, tm$Dam, tm[,4:9])
head(test.data)
names(test.data)
names(test.data) <- c("id", "pop","fam","CorW","CorL","PetL","Col","NecA","NecC")
names(data)

## for data
n=dim(data)[2]-3;nn=table(data$pop);nbpop=length(nn);
traits=names(data[-c(1:3)]);
effect=paste(traits,collapse=",")
form=as.formula(paste("cbind(",effect,")~factor(pop)",sep="")) #automated formula for the manova
manov=summary(manova(form,data=data))                   #do the manova
SSb=manov$SS$"factor(pop)";                             #between pop Sum Squares (SS)
SSw=manov$SS$Residuals;                                 #within pop SS
dfb=manov$stats[1,1];dfw=manov$stats[2,1];              #df between/within pops
nf=mean(nn)-1/nbpop*((mean(nn^2)-mean(nn)^2)/mean(nn)); #compute the equivalent df with unbalanced designs

## for test.data
n=dim(test.data)[2]-3;nn=table(test.data$pop);nbpop=length(nn);
traits=names(test.data[-c(1:3)]);
effect=paste(traits,collapse=",")
form=as.formula(paste("cbind(",effect,")~factor(pop)",sep="")) #automated formula for the manova
manov=summary(manova(form,data=test.data))                   #do the manova
SSb=manov$SS$"factor(pop)";                             #between pop Sum Squares (SS)
SSw=manov$SS$Residuals;                                 #within pop SS
dfb=manov$stats[1,1];dfw=manov$stats[2,1];              #df between/within pops
nf=mean(nn)-1/nbpop*((mean(nn^2)-mean(nn)^2)/mean(nn)); #compute the equivalent df with unbalanced designs

# half-sib test
Gw4=4*MSw;Gb=(MSb-Gw4)/nf;    
DF=c(dfw,dfb);G4=list(Gw4,Gb);MS=list(MSw,MSb); 

# half-sib test2
MSw=(4*SSw)/dfw;MSb=SSb/dfb;                                #compute the corresponding Mean Square matrices
Gw=MSw;Gb=(MSb-Gw)/nf;                                  #compute the corresponding G matrices
DF=c(dfw,dfb);G=list(Gw,Gb);MS=list(MSw,MSb);           #create the inputs for k.prop()

# Correlations
Gwcor<- cov2cor(Gw)
Gw4cor<- cov2cor(Gw4)
Gbcor<- cov2cor(Gb)
Gwcor
Gw4cor
Gbcor
#########end testing

model <- summary(manova(data~factor))
# Perform MANOVA
SSb <- model$SS$factor
#Extract sum of squares between provenances
SSw <- model$SS$Residuals
#Extract sum of squares residuals (between families)
dfb <- model$stats[1,1]
# Degree of freedom between provenances
dfw <- model$stats[2,1]
# Degree of freedom within provenances (between families)
nn <- table(factor)
# Extract number of families per provenance
nbpop<- length(nn) # number of provenances
nf=mean(nn)-1/nbpop*((mean(nn^2)-mean(nn)^2)/mean(nn))
#compute the equivalent df with unbalanced designs
MSw=SSw/dfw
# Matrix of mean squares within
MSb=SSb/dfb
# Matrix of mean squares between
Gw <- MSw
# Matrix of Var-Cov within
Gb=(MSb-Gw)/nf
# Matrix of VarCov between
DF <- c(dfw,dfb)
# Vector with df within and between
G <- list(Gw,Gb)
# List of var-cov matrices
MS <- list(MSw,MSb)
# List of mean square matrices
source("package neutrality test.r")
# Read Martin's script (Check dependencies first: MASS, KernSmooth, corpcor
testG <- k.prop(DF,G)
# Test G
paste("rho_P (population) =",signif(testG$rho[[2]],2))
# you will get rho estimate corresponding to your dataset
#[1] "rho_P (factor) = 0.28"
print("95% CI for rho_P :");signif(testG$CI[[2]],2)
# you will get the CI for your rho estimate
#[1] "95% CI for rho_P :"
#[1] 0.24 0.34
testMS=k.prop(DF,MS)
# Proportionality of matrices
paste("pBartlett=",signif(testMS$pt1,3),"pChi²=",signif(testMS$pX,3))
#[1] "pBartlett= 0.00823 pChi²= 3.13e-05"

save.image("inprocess_april2.RData")

##########################
## bootstrap Fst to calculate confidence intervals for microsat rho.

# genind object
mydata

##########################
##########################

