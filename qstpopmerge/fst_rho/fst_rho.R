setwd("~/Dropbox/QSTFST/qstpopmerge/fst_rho/")
load("rho_fst.RData")

rho<-(2*fst.out)/(1-fst.out)
rho.05<-rho[rank(rho)==length(rho)*.05]
rho.95<-rho[rank(rho)==length(rho)*.95]
hist(rho, breaks=100)
abline(v=rho.05, col="red")
abline(v=rho.95, col="blue")

ave(rho)
mean(rho)
