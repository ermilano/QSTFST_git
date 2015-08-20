## driftsel
## R package for detecting signals of natural selection in quantitative traits

## Working with Ipomopsis matched population datasets
## April 2015

setwd("~/Dropbox/QSTFST/match/driftsel/") 

# install.packages("MCMCpack")
# install.packages("SparseM")
# download packates from http://www.helsinki.fi/biosci/egru/software/driftsel.html
# install.packages("~/Downloads/driftsel_2.1.2.tar.gz", repos = NULL, type = "source")
# install.packages("~/Downloads/RAFM_1.2.tar.gz", repos = NULL, type = "source", 
#                  lib="/Library/Frameworks/R.framework/Versions/3.1/Resources/library")
library(driftsel)
library(RAFM)

#################################################################
# This package is not very well documented and there are very specific 
# data format requirements. In brief, there are four files necessary: 
#   genotypes, phenotypes, covariates and pedigree
# The covariate file is necessary even if no covariates are being tested.
# Rquirements for the phenotype, covariate and pedigree files:
#   *id is the first column and in the same order across all files
#   *matrix structure, where dimnames=list(NULL, colnames(data.frame))
#   *only numeric values or NA (ie. sires, dams and populations 
#        must be numeric)
# Populations in the genotype file must equal those in the pedigree 
#        (haven't tested this explicitly but I am making the assumption.)

# Workflow:
# -Load genotype data
# -Run using RAFM::do.all()
#   *test with small number of iterations
#   *figure out how to confirm MCMC convergence 
#   *save the final object because the full MCMC run takes a long time
# -Visualize populations, including the ancestral, with viz.theta()
# 
# -Load phenotype and pedigree files
#   *transform and scale if necessary
#   *create dummy covariate file using the first column of the 
#       traits file (of use real covariates file) 
#   *convert into matrices with the proper dimnames
# -Run MH()
#   *test with small number of iterations
#   *figure out how to confirm MCMC convergence
# -Test for neutrality, neut.test()
#   *0 =/< S < .05 indicates stabilizing selection
#   *S =.5 indicates neutrality
#   *.95 < S </= 1 indicates divergent selection
# -Visualize traits, including the ancestral, with viz.traits()
#   *can only view two at a time
#   *adjust confidence ellipse
#################################################################

#### Genotype Anaylsis ####

geno <- read.csv("geno_ds_match.csv")
head(geno)
afm <- do.all(geno, 20000, 10000, 10)  
# save(afm, file="afm_20000.RData")
load("afm_20000.RData")

names(afm)   
summary(afm$fst)   # range of Fst
par(mfcol=c(1,1))  # do.all messes with panels

# Visualize
viz.theta(afm$theta)   # pattern of relatedness, 'theta'
viz.theta(afm$theta, distance=NA) # without lines
head(afm$theta)
is.array(afm$theta)

## Coancestry matrix, theta
heat.theta(afm$theta)
stats <- heat.theta(afm$theta, mean=FALSE)
round(stats$median, 4)


pdf("theta_coancestry.pdf")
heat.theta(afm$theta)
viz.theta(afm$theta)
viz.theta(afm$theta, distance=NA)
dev.off()

#### Phenotype Anaylsis ####

## Pedigree
ped <- read.csv("ped_ds_match_na.csv")
ped<-ped[,1:5]
head(ped)   
colnames(ped)
pedm <- as.matrix(ped, dimnames=list(NULL,colnames(ped)), quote=FALSE)
head(pedm)
## Traits
traits <- read.csv("box.scale.ds.match.csv")
traits_scale <- scale(traits[2:9])
traits_scale <- cbind(traits[1], traits_scale)
head(traits_scale)
traits <- traits_scale
apply(traits, 2, mean, na.rm=TRUE)
colnames(traits)
traitsm<-as.matrix(traits,dimnames=list(NULL,colnames(traits)), quote=FALSE)
head(traitsm) 
dimnames(traitsm)
## Covariates
covarsm <- as.matrix(traitsm[,1],  dimnames=list(NULL,"id"))
dimnames(covarsm) <- list(NULL, "id")

## Check file formats
is.array(afm$theta)
is.matrix(pedm)
is.matrix(traitsm)
is.matrix(covarsm)

## Metropolis-Hastings algorithm
samp <- MH(afm$theta, pedm, covarsm, traitsm, 5000, 2000, 10, alt=T) 
names(samp) 	
#save.image("samp5000.RData")
load("samp5000.RData")
## Test for signal of selection
neut.test(samp$pop.ef, samp$G, samp$theta,silent=F)  

## Visualize trait space with neutral ellipse expectations
## Try different trait and confidence interaval values
viz.traits(samp$fixed.ef, samp$pop.ef, samp$G, samp$theta, c(1,2), main="Trait space", siz=.95)
viz.traits(samp$fixed.ef, samp$pop.ef, samp$G, samp$theta, c(1,3), main="Trait space", siz=.95)
viz.traits(samp$fixed.ef, samp$pop.ef, samp$G, samp$theta, c(1,4), main="Trait space", siz=.95)
viz.traits(samp$fixed.ef, samp$pop.ef, samp$G, samp$theta, c(1,5), main="Trait space", siz=.95)
viz.traits(samp$fixed.ef, samp$pop.ef, samp$G, samp$theta, c(1,6), main="Trait space", siz=.95)
viz.traits(samp$fixed.ef, samp$pop.ef, samp$G, samp$theta, c(1,7), main="Trait space", siz=.95)
viz.traits(samp$fixed.ef, samp$pop.ef, samp$G, samp$theta, c(1,8), main="Trait space", siz=.95)

viz.traits(samp$fixed.ef, samp$pop.ef, samp$G, samp$theta, c(2,3), main="Trait space", siz=.95)
pdf("driftsel2.4.pdf")
viz.traits(samp$fixed.ef, samp$pop.ef, samp$G, samp$theta, c(2,4), main="Trait space: Corolla Length and Color", siz=.95)
dev.off()
viz.traits(samp$fixed.ef, samp$pop.ef, samp$G, samp$theta, c(2,5), main="Trait space", siz=.95)
viz.traits(samp$fixed.ef, samp$pop.ef, samp$G, samp$theta, c(2,6), main="Trait space", siz=.95)
viz.traits(samp$fixed.ef, samp$pop.ef, samp$G, samp$theta, c(2,7), main="Trait space", siz=.95)
viz.traits(samp$fixed.ef, samp$pop.ef, samp$G, samp$theta, c(2,8), main="Trait space", siz=.95)

viz.traits(samp$fixed.ef, samp$pop.ef, samp$G, samp$theta, c(3,4), main="Trait space", siz=.95)
viz.traits(samp$fixed.ef, samp$pop.ef, samp$G, samp$theta, c(3,5), main="Trait space", siz=.95)
viz.traits(samp$fixed.ef, samp$pop.ef, samp$G, samp$theta, c(3,6), main="Trait space", siz=.95)
viz.traits(samp$fixed.ef, samp$pop.ef, samp$G, samp$theta, c(3,7), main="Trait space", siz=.95)
viz.traits(samp$fixed.ef, samp$pop.ef, samp$G, samp$theta, c(3,8), main="Trait space", siz=.95)

viz.traits(samp$fixed.ef, samp$pop.ef, samp$G, samp$theta, c(4,5), main="Trait space", siz=.95)
viz.traits(samp$fixed.ef, samp$pop.ef, samp$G, samp$theta, c(4,6), main="Trait space", siz=.95)
viz.traits(samp$fixed.ef, samp$pop.ef, samp$G, samp$theta, c(4,7), main="Trait space", siz=.95)
viz.traits(samp$fixed.ef, samp$pop.ef, samp$G, samp$theta, c(4,8), main="Trait space", siz=.95)

viz.traits(samp$fixed.ef, samp$pop.ef, samp$G, samp$theta, c(5,6), main="Trait space", siz=.95)
viz.traits(samp$fixed.ef, samp$pop.ef, samp$G, samp$theta, c(5,7), main="Trait space", siz=.95)
viz.traits(samp$fixed.ef, samp$pop.ef, samp$G, samp$theta, c(5,8), main="Trait space", siz=.95)

viz.traits(samp$fixed.ef, samp$pop.ef, samp$G, samp$theta, c(6,7), main="Trait space", siz=.95)
viz.traits(samp$fixed.ef, samp$pop.ef, samp$G, samp$theta, c(6,8), main="Trait space", siz=.95)

viz.traits(samp$fixed.ef, samp$pop.ef, samp$G, samp$theta, c(7,8), main="Trait space", siz=.95)

# save.image("inprogress_april3.RData")
# load("inprogress_april3.RData")


#################### Just floral traits #######################
## Metropolis-Hastings algorithm

head(traitsm)
trait.f <- traitsm[,1:7]
head(trait.f)

samp.f <- MH(afm$theta, pedm, covarsm, trait.f, 5000, 2000, 10, alt=T) 
names(samp.f)   
#save.image("samp.f5000.RData")
load("samp.f5000.RData")
## Test for signal of selection
neut.test(samp.f$pop.ef, samp.f$G, samp.f$theta,silent=F)  

## Visualize trait space with neutral ellipse expectations
## Try different trait and confidence interaval values
viz.traits(samp.f$fixed.ef, samp.f$pop.ef, samp.f$G, samp.f$theta, c(1,2), main="Trait space", siz=.95)
viz.traits(samp.f$fixed.ef, samp.f$pop.ef, samp.f$G, samp.f$theta, c(1,3), main="Trait space", siz=.95)
viz.traits(samp.f$fixed.ef, samp.f$pop.ef, samp.f$G, samp.f$theta, c(1,4), main="Trait space", siz=.95)
viz.traits(samp.f$fixed.ef, samp.f$pop.ef, samp.f$G, samp.f$theta, c(1,5), main="Trait space", siz=.95)
viz.traits(samp.f$fixed.ef, samp.f$pop.ef, samp.f$G, samp.f$theta, c(1,6), main="Trait space", siz=.95)

viz.traits(samp.f$fixed.ef, samp.f$pop.ef, samp.f$G, samp.f$theta, c(2,3), main="Trait space", siz=.95)
#pdf("driftsel2.4_flower.pdf")
viz.traits(samp.f$fixed.ef, samp.f$pop.ef, samp.f$G, samp.f$theta, c(2,4), main="Trait space: Corolla Length and Color", siz=.95)
#dev.off()
viz.traits(samp.f$fixed.ef, samp.f$pop.ef, samp.f$G, samp.f$theta, c(2,5), main="Trait space", siz=.95)
viz.traits(samp.f$fixed.ef, samp.f$pop.ef, samp.f$G, samp.f$theta, c(2,6), main="Trait space", siz=.95)

viz.traits(samp.f$fixed.ef, samp.f$pop.ef, samp.f$G, samp.f$theta, c(3,4), main="Trait space", siz=.95)
viz.traits(samp.f$fixed.ef, samp.f$pop.ef, samp.f$G, samp.f$theta, c(3,5), main="Trait space", siz=.95)
viz.traits(samp.f$fixed.ef, samp.f$pop.ef, samp.f$G, samp.f$theta, c(3,6), main="Trait space", siz=.95)

viz.traits(samp.f$fixed.ef, samp.f$pop.ef, samp.f$G, samp.f$theta, c(4,5), main="Trait space", siz=.95)
viz.traits(samp.f$fixed.ef, samp.f$pop.ef, samp.f$G, samp.f$theta, c(4,6), main="Trait space", siz=.95)

viz.traits(samp.f$fixed.ef, samp.f$pop.ef, samp.f$G, samp.f$theta, c(5,6), main="Trait space", siz=.95)



### final figure cube
par(mfrow=c(2,2))
viz.theta(afm$theta, distance=NA)
viz.traits(samp.f$fixed.ef, samp.f$pop.ef, samp.f$G, samp.f$theta, c(2,4), main="Trait space: Corolla Width and Color", siz=.95)
viz.traits(samp.f$fixed.ef, samp.f$pop.ef, samp.f$G, samp.f$theta, c(1,5), main="Trait space: Corolla Length and Nectar Volume", siz=.95)
viz.traits(samp.f$fixed.ef, samp.f$pop.ef, samp.f$G, samp.f$theta, c(3,6), main="Trait space: Petal Length and Nectar Concentration", siz=.95)
dev.off()
#newviz(samp.f$fixed.ef, samp.f$pop.ef, samp.f$G, samp.f$theta, c(3,6), main="Trait space: Petal Length and Nectar Concentration", siz=.95, xlab="Petal Length", ylab="Nectar Concentration", xlim=c(-1.3,1.3), ylim=c(-1,1))

par(mfrow=c(2,2))
viz.theta(afm$theta, distance=NA)
viz.traits(samp.f$fixed.ef, samp.f$pop.ef, samp.f$G, samp.f$theta, c(2,4), main="Trait space: Corolla Width and Color", siz=.5)
viz.traits(samp.f$fixed.ef, samp.f$pop.ef, samp.f$G, samp.f$theta, c(1,5), main="Trait space: Corolla Length and Nectar Volume", siz=.5)
viz.traits(samp.f$fixed.ef, samp.f$pop.ef, samp.f$G, samp.f$theta, c(3,6), main="Trait space: Petal Length and Nectar Concentration", siz=.5)
dev.off()

pdf("driftsel_95_quad.pdf")
par(mfrow=c(2,2))
viz.theta(afm$theta)
viz.traits(samp.f$fixed.ef, samp.f$pop.ef, samp.f$G, samp.f$theta, c(2,4), main="Trait space: Corolla Width and Color", siz=.95)
viz.traits(samp.f$fixed.ef, samp.f$pop.ef, samp.f$G, samp.f$theta, c(1,2), main="Trait space: Corolla Length and Corolla Width", siz=.95)
viz.traits(samp.f$fixed.ef, samp.f$pop.ef, samp.f$G, samp.f$theta, c(3,6), main="Trait space: Petal Length and Nectar Concentration", siz=.95)
dev.off()

source('~/Dropbox/QSTFST/match/driftsel/driftsel_modified_figures.R')
pdf("driftsel_50_quad.pdf")
par(mfrow=c(2,2))
viz.theta2(afm$theta, distance=NA)
viz.traits2(samp.f$fixed.ef, samp.f$pop.ef, samp.f$G, samp.f$theta, c(2,4), xlab="Corolla Width", ylab="Color",  siz=.5)
viz.traits2(samp.f$fixed.ef, samp.f$pop.ef, samp.f$G, samp.f$theta, c(1,2), xlab="Corolla Length", ylab="Corolla Width", siz=.5)
viz.traits2(samp.f$fixed.ef, samp.f$pop.ef, samp.f$G, samp.f$theta, c(3,6), xlab="Petal Length", ylab="Nectar Concentration",siz=.5)
dev.off()