# G and D matricies using MANOVA code from Martin et al. 2008
# starting with box-cox transformed and scaled traits from merged qst pops.


setwd("~/Dropbox/QSTFST/qstpopmerge/")
# library(expm) # for covariance matrix decomposition
library(KernSmooth) # for package neutrality test.r
library(corpcor) # for package neutrality test.r
library(MASS) # for package neutrality test.r
source("~/Dropbox/QSTFST/qstpopmerge/package neutrality test.r")

box.scale.na <- read.csv("box.scale.na.csv")
box.scale <- read.csv("box.scale.csv")
names(box.scale.na) <- c("Id", "Pop", "Dam", "CorollaWidth", "CorollaLength", "PetalLength",
                        "Color", "NectarAmt", "NectarConc", "Root", "Fdate")
names(box.scale) <- c("Id","Pop","Dam","CorollaWidth","CorollaLength","PetalLength",
                      "Color","NectarAmt","NectarConc","Root","Fdate")


#####################
# G and D matricies using MANOVA code from Bessega et al. 2015
# Only floral traits

dam.mean <- aggregate(box.scale[,4:9], list(box.scale$Dam), mean, na.rm=T)
dam.mean
pop.dam <- unique(box.scale[,c(2,3)])
trait.mean <- merge(pop.dam, dam.mean, by.y="Group.1", by.x="Dam")

write.csv(trait.mean, "trait.mean.csv")
tm <- read.csv("trait.mean.csv")
head(tm)
# data <- as.matrix(tm[,4:9])
# factor <- as.matrix(tm[,3])

#######testing
head(tm)
test.data <- cbind(tm$X, tm$Pop, tm$Dam, tm[,4:9])
head(test.data)
names(test.data)
names(test.data) <- c("id", "pop","fam","CorW","CorL","PetL","Col","NecA","NecC")
head(test.data)

## for data
# n=dim(data)[2]-3;nn=table(data$pop);nbpop=length(nn);
# traits=names(data[-c(1:3)]);
# effect=paste(traits,collapse=",")
# form=as.formula(paste("cbind(",effect,")~factor(pop)",sep="")) #automated formula for the manova
# manov=summary(manova(form,data=data))                   #do the manova
# SSb=manov$SS$"factor(pop)";                             #between pop Sum Squares (SS)
# SSw=manov$SS$Residuals;                                 #within pop SS
# dfb=manov$stats[1,1];dfw=manov$stats[2,1];              #df between/within pops
# nf=mean(nn)-1/nbpop*((mean(nn^2)-mean(nn)^2)/mean(nn)); #compute the equivalent df with unbalanced designs

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

MSw=SSw/dfw;MSb=SSb/dfb;                                #compute the corresponding Mean Square matrices
Gw=MSw;Gb=(MSb-Gw)/nf;                                  #compute the corresponding G matrices
DF=c(dfw,dfb);G=list(Gw,Gb);MS=list(MSw,MSb);           #create the inputs for k.prop()

# half-sib test
Gw4=4*MSw;Gb=(MSb-Gw4)/nf;    
DF=c(dfw,dfb);G4=list(Gw4,Gb);MS=list(MSw,MSb); 

# Correlations
Gwcor<- cov2cor(Gw)
Gw4cor<- cov2cor(Gw4)
Gbcor<- cov2cor(Gb)
Gwcor
Gw4cor
Gbcor

###############
# bending test

# Here is a function that I use for symmetric matrices:
# https://stat.ethz.ch/pipermail/r-help/2003-December/042964.html  

make.positive.definite <-
  function(x, tol=1e-6) {
    eig <- eigen(x, symmetric=TRUE)
    rtol <- tol * eig$values[1]
    if(min(eig$values) < rtol) {
      vals <- eig$values
      vals[vals < rtol] <- rtol
      srev <- eig$vectors %*% (vals * t(eig$vectors))
      dimnames(srev) <- dimnames(x)
      return(srev)
    } else {
      return(x)
    }
  }
nGb <- make.positive.definite(Gb)

# Correlations
Gwcor<- cov2cor(Gw)
Gw4cor<- cov2cor(Gw4)
nGbcor<- cov2cor(nGb)
Gwcor
Gw4cor
nGbcor

######### proportionality testing

#rho_st from test between G matrices
testG=k.prop(DF,G)
paste("rho_P (population) =",signif(testG$rho[[2]],2))
print("95% CI for rho_P :");signif(testG$CI[[2]],2)
#half-sib
testG4=k.prop(DF,G4)
paste("rho_P (population) =",signif(testG4$rho[[2]],2))
print("95% CI for rho_P :");signif(testG4$CI[[2]],2)

#p-values from test between MS matrices
testMS=k.prop(DF,MS);
paste("pBartlett=",signif(testMS$pt1,3),"pChi≤=",signif(testMS$pX,3))


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

# save.image("inprocess_april2.RData")

##########################

# decomposition of Gw, Gb and Gw4

gw.out <- princomp(covmat=Gw)
names(gw.out)
gw.out$loadings
summary(gw.out)
plot(gw.out) #scree plot

gw4.out <- princomp(covmat=Gw4)
names(gw4.out)
gw4.out$loadings
summary(gw4.out)
plot(gw4.out) #scree plot

gb.out <- princomp(covmat=Gb)
names(gb.out)
gb.out$loadings
summary(gb.out)
plot(gb.out) #scree plot

###############
# bending test

# Here is a function that I use for symmetric matrices:
# https://stat.ethz.ch/pipermail/r-help/2003-December/042964.html  

  make.positive.definite <-
  function(x, tol=1e-6) {
    eig <- eigen(x, symmetric=TRUE)
    rtol <- tol * eig$values[1]
    if(min(eig$values) < rtol) {
      vals <- eig$values
      vals[vals < rtol] <- rtol
    srev <- eig$vectors %*% (vals * t(eig$vectors))
    dimnames(srev) <- dimnames(x)
    return(srev)
  } else {
    return(x)
  }
}

nGb <- make.positive.definite(Gb)

gb.out <- princomp(covmat=nGb)
names(gb.out)
gb.out$loadings
summary(gb.out)
plot(gb.out) #scree plot

###############

eigen(Gw)
eigen(Gw)$values/sum(eigen(Gw)$values)
cumsum(eigen(Gw)$values/sum(eigen(Gw)$values))

eigen(Gb)
eigen(Gb)$values/sum(eigen(Gb)$values)
cumsum(eigen(Gb)$values/sum(eigen(Gb)$values))

#Gb

eig <- eigen(nGb)
eigv <- eig$vectors
abs(eigv[,1])/sum(abs(eigv[,1]))
abs(eigv[,1])/max(abs(eigv[,1]))
cumsum(eig$values/sum(eig$values))
eig$values/sum(eig$values)
pdf("Gb_pc.pdf")
par(mfrow=c(2,2), mar=c(2.5,2.5,2.5,2.5),oma = c(0, 0, 2, 0))
barplot(abs(eigv[,1]),names=1:6, main="PC1")
barplot(abs(eigv[,2]),names=1:6, main="PC2")
barplot(abs(eigv[,3]),names=1:6, main="PC3")
barplot(abs(eigv[,4]),names=1:6, main="PC4")
mtext("Trait Loadings on PC 1-4: Gb", outer = TRUE, cex = 1.2)
dev.off()


#Gw
eig <- eigen(Gw)
eigv <- eig$vectors
abs(eigv[,1])/sum(abs(eigv[,1]))
abs(eigv[,1])/max(abs(eigv[,1]))
cumsum(eig$values/sum(eig$values))
eig$values/sum(eig$values)
pdf("Gw_pc.pdf")
par(mfrow=c(2,2), mar=c(2.5,2.5,2.5,2.5),oma = c(0, 0, 2, 0))
barplot(abs(eigv[,1]),names=1:6, main="PC1")
barplot(abs(eigv[,2]),names=1:6, main="PC2")
barplot(abs(eigv[,3]),names=1:6, main="PC3")
barplot(abs(eigv[,4]),names=1:6, main="PC4")
mtext("Trait Loadings on PC 1-4: Gw", outer = TRUE, cex = 1.2)
dev.off()


#Gw
eig <- eigen(Gw4)
eigv <- eig$vectors
abs(eigv[,1])/sum(abs(eigv[,1]))
abs(eigv[,1])/max(abs(eigv[,1]))
cumsum(eig$values/sum(eig$values))
eig$values/sum(eig$values)
pdf("Gw_pc.pdf")
par(mfrow=c(2,2), mar=c(2.5,2.5,2.5,2.5),oma = c(0, 0, 2, 0))
barplot(abs(eigv[,1]),names=1:6, main="PC1")
barplot(abs(eigv[,2]),names=1:6, main="PC2")
barplot(abs(eigv[,3]),names=1:6, main="PC3")
barplot(abs(eigv[,4]),names=1:6, main="PC4")
mtext("Trait Loadings on PC 1-4: Gw", outer = TRUE, cex = 1.2)
dev.off()

##################Combined eigen loadings


eigb <- eigen(nGb)
eigvb <- eigb$vectors
abs(eigvb[,1])/sum(abs(eigvb[,1]))
abs(eigvb[,1])/max(abs(eigvb[,1]))
cumsum(eigb$values/sum(eigb$values))
eigb$values/sum(eigb$values)
#Gw
eigw <- eigen(Gw)
eigvw <- eigw$vectors
abs(eigvw[,1])/sum(abs(eigvw[,1]))
abs(eigvw[,1])/max(abs(eigvw[,1]))
cumsum(eigw$values/sum(eigw$values))
eigw$values/sum(eigw$values)

color <- c("darkred", "darkred", "darkred", "darkred", "yellow", "yellow")
name <- c("CW","CL","PL","Col","NV","NC")

# with absolute value
pdf("G_loadings_abs.pdf", width=8, height=6)
par(mfrow=c(2,2), mar=c(2.5,3,3,3),oma = c(1, 1, 4, 1))
barplot(abs(eigvw[,1]),names=name, main="Gw PC1: 39.3%", col=color, cex.names=.7)
barplot(abs(eigvw[,2]),names=name, main="Gw PC2: 22.2%", col=color, cex.names=.7)
#barplot(abs(eigvw[,3]),names=name, main="Gw PC3: 15.3%", col=color, cex.names=.7)
#barplot(abs(eigvw[,4]),names=name, main="Gw PC4: 10.4%", col=color, cex.names=.7)
barplot(abs(eigvb[,1]),names=name, main="Gb PC1: 90.1%", col=color, cex.names=.7)
barplot(abs(eigvb[,2]),names=name, main="Gb PC2: 9.2%", col=color, cex.names=.7)
mtext("PC trait toadings for G-within and G-between", outer = TRUE, cex = 1.5)
dev.off()
##########################
# without absolute value
## testing to try and recover negative pc components for Gb
pdf("G_loadings_neg.pdf", width=8, height=6)
par(mar=c(2.5,3,3,3),oma = c(1, 1, 4, 1))
layout(matrix(c(1, 2, 3, 4),2,2, byrow=TRUE),heights=c(1, 1.3))
barplot(abs(eigvw[,1]),names=name, main="Gw PC1: 39.3%", col=color, cex.names=.7)
barplot(abs(eigvw[,2]),names=name, main="Gw PC2: 22.2%", col=color, cex.names=.7)
barplot(eigvb[,1],names=name, main="Gb PC1: 90.1%", col=color, cex.names=.7)
barplot(eigvb[,2],names=name, main="Gb PC2: 9.2%", col=color, cex.names=.7)
mtext("PC trait toadings for G-within and G-between", outer = TRUE, cex = 1.5)
dev.off()


##########################
## bootstrap Fst to calculate confidence intervals for microsat rho.

# genind object
# mydata


##################
# compare loadings by correlation
a <- abs(eigvb[,1])/sum(abs(eigvb[,1]))
b <- abs(eigvw[,1])/sum(abs(eigvw[,1]))
cor(a,b)
plot(a,b)
c <- abs(eigvb[,2])/sum(abs(eigvb[,2]))
d <- abs(eigvw[,2])/sum(abs(eigvw[,2]))
cor(c,d)
plot(c,d)

# noabs
a <- (eigvb[,1])/sum((eigvb[,1]))
b <- (eigvw[,1])/sum((eigvw[,1]))
cor(a,b)
plot(a,b)
c <- (eigvb[,2])/sum((eigvb[,2]))
d <- (eigvw[,2])/sum((eigvw[,2]))
cor(c,d)
plot(c,d)

# plot G matrix ellipses
plot(a,c, xlim=c(-6,6), ylim=c(-5,5), col="red", pch=19)
points(b, d, pch=19, col="blue")

dataEllipse(a,c, levels=.95, plot.points=FALSE, col="red")
dataEllipse(b,d, levels=.95, plot.points=FALSE, col="blue")

corgb <- cor(a,c)
test <- cbind(a,c)
dataEllipse(test)
dataEllipse(b, d)
plotcorr(nGbcor)

#####
# need to get populations means from within in the manova to plot onto PCs of G
manov=summary(manova(form,data=test.data))                   #do the manova
SSb=manov$SS$"factor(pop)";                             #between pop Sum Squares (SS)
SSw=manov$SS$Residuals;                                 #within pop SS







##########################
# save.image("rho_fst.RData")
# 
# load("fst_rho/rho_fst.RData")


### sub-script in plot lables
pdf("G_loadings_abs.pdf", width=8, height=6)
par(mfrow=c(2,2), mar=c(2.5,3,3,3),oma = c(1, 1, 4, 1))
barplot(abs(eigvw[,1]),names=name, main=expression("G"[W]*" PC1: 39.3%"), col=color, cex.names=1, cex.main = 1.3)
barplot(abs(eigvw[,2]),names=name, main=expression("G"[W]*" PC2: 22.2%"), col=color, cex.names=1, cex.main = 1.3)
barplot(abs(eigvb[,1]),names=name, main=expression("G"[B]*" PC1: 90.1%"), col=color, cex.names=1, cex.main = 1.3)
barplot(abs(eigvb[,2]),names=name, main=expression("G"[B]*" PC2: 9.2%"), col=color, cex.names=1, cex.main = 1.3)
mtext("PC trait toadings for G-within and G-between", outer = TRUE, cex = 1.5)
dev.off()
