QstFstComp_fst.out<-function (fst.dat, qst.dat, numpops, nsim = 1000, AFLP = FALSE, 
          breeding.design, dam.offspring.relatedness = 0.25, output = "concise") 
{
  if (missing(fst.dat)) 
    stop("Genotypic data must be provided.")
  if (missing(qst.dat)) 
    stop("Phenotypic data must be provided.")
  if (missing(numpops)) 
    stop("Number of populations must be defined.")
  if (missing(breeding.design)) 
    stop("Breeding design must be defined.")
  if (breeding.design == "half.sib.sire") {
    qst.MS <- MeanSq.unbalanced.sire(qst.dat)
    qst.temp <- QSTfromSireModel(qst.MS$MSpops, qst.MS$MSsires, 
                                 qst.MS$MSdams, qst.MS$MSwithin, qst.MS$n0primeprime, 
                                 qst.MS$n0prime, qst.MS$n0, qst.MS$nc0prime, qst.MS$nc0, 
                                 qst.MS$ncb0)
    qst.obs <- qst.temp$Qst
    mean.trait.value <- mean(qst.dat[, 4], na.rm = TRUE)
    Va <- qst.temp$Va
    if (Va < 0) {
      Va <- 0
    }
    CVa <- sqrt(Va)/mean.trait.value * 100
  }
  if (breeding.design == "half.sib.dam") {
    qst.MS <- MeanSq.unbalanced.dam(qst.dat)
    qst.obs <- QSTfromDamModel(qst.MS$MSpops, qst.MS$MSdams, 
                               qst.MS$MSwithin, qst.MS$n0prime, qst.MS$n0, qst.MS$nb0, 
                               dam.offspring.relatedness)
    mean.trait.value <- mean(qst.dat[, 3], na.rm = TRUE)
    Va <- 1/dam.offspring.relatedness * (qst.MS$MSdams - 
                                           qst.MS$MSwithin)/qst.MS$n0
    if (Va < 0) {
      Va <- 0
    }
    CVa <- sqrt(Va)/mean.trait.value * 100
  }
  if (AFLP == FALSE) {
    nloci <- ncol(fst.dat) - 1
    abc.mat <- wc.calc(fst.dat, diploid = TRUE)
    nalleles <- nrow(abc.mat)
    fst.obs <- sum(abc.mat[, 1])/sum(abc.mat[, 1] + abc.mat[, 
                                                            2] + abc.mat[, 3])
  }
  if (AFLP == TRUE) {
    full.fst.dat <- read.fst.input(fst.dat, num.pops = numpops)
    fst.data <- as.matrix(full.fst.dat[[1]])
    var.dat <- as.matrix(full.fst.dat[[2]])
    nloci <- nrow(fst.data)
    fst.obs <- mean.fst(dat = fst.data, vardat = var.dat, 
                        num.loci = nloci, num.pops = numpops)
  }
  sim.est <- vector(length = nsim)
  qst.neut <- vector(length = nsim)
  fst.est <- vector(length = nsim)
  qstForCI.est <- vector(length = nsim)
  Va.est <- vector(length = nsim)
  CVa.est <- vector(length = nsim)
  for (i in 1:nsim) {
    if (AFLP == FALSE) {
      fst.repl <- fst.sample(abc.mat, nalleles)
    }
    if (AFLP == TRUE) {
      fst.repl <- fst.sample.aflp(fst.data, var.dat, nloci)
    }
    if (breeding.design == "half.sib.sire") 
      qst.repl <- qst.parboot.siremodel(qst.MS, fst.obs)[[1]]
    if (breeding.design == "half.sib.dam") 
      qst.repl <- qst.parboot.dammodel(qst.MS, fst.obs, 
                                       dam.offspring.relatedness)
    sim.est[i] <- qst.repl - fst.repl
    fst.est[i] <- fst.repl
    qst.neut[i] <- qst.repl
    if (breeding.design == "half.sib.sire") {
      temp <- qstVa.parbootForCI.siremodel(qst.MS)
    }
    if (breeding.design == "half.sib.dam") {
      temp <- qstVa.parbootForCI.dammodel(qst.MS, dam.offspring.relatedness)
    }
    qstForCI.est[i] <- temp[1]
    Va.est[i] <- temp[2]
    if (Va.est[i] < 0) {
      Va.est[i] <- 0
    }
    CVa.est[i] <- sqrt(Va.est[i])/mean.trait.value * 100
  }
  fst.out<<-fst.est
  time.stamp <- Sys.time()
  formatted.time.stamp <- gsub(" ", "_", time.stamp)
  formatted.time.stamp <- gsub(":", "-", formatted.time.stamp)
  writeLines(as.character(sim.est), paste("QminusFvalues_", 
                                          formatted.time.stamp, ".txt", sep = ""))
  diff.repl <- qst.neut - fst.est
  Q.obsMinusF.obs <- qst.obs - fst.obs
  right.one.tailed.p <- sum(Q.obsMinusF.obs < diff.repl)/length(diff.repl)
  left.one.tailed.p <- sum(Q.obsMinusF.obs > diff.repl)/length(diff.repl)
  two.tailed.p <- 2 * min(right.one.tailed.p, left.one.tailed.p)
  if (output == "concise") {
    return(list(QminusF <- c(`Calculated Qst-Fst` = Q.obsMinusF.obs, 
                             `Lower Bound crit. value` = quantile(sim.est, 0.025, 
                                                                  na.rm = TRUE), `Upper bound crit. value` = quantile(sim.est, 
                                                                                                                      0.975, na.rm = TRUE)), Distribution <- paste("Qst-Fst values output to file:", 
                                                                                                                                                                   paste(getwd(), "/QminusFvalues_", formatted.time.stamp, 
                                                                                                                                                                         ".txt", sep = "")), QminusF.p.values <- c(`Lower one-tailed p value` = left.one.tailed.p, 
                                                                                                                                                                                                                   `Upper one-tailed p value` = right.one.tailed.p, 
                                                                                                                                                                                                                   `Two-tailed p value` = two.tailed.p), Fst <- c(`Estimated Fst` = fst.obs, 
                                                                                                                                                                                                                                                                  `Lower Bound CI` = quantile(fst.est, 0.025, na.rm = TRUE), 
                                                                                                                                                                                                                                                                  `Upper bound CI` = quantile(fst.est, 0.975, na.rm = TRUE)), 
                Qst <- c(`Estimated Qst` = qst.obs, `Lower Bound CI` = quantile(qstForCI.est, 
                                                                                0.025, na.rm = TRUE), `Upper bound CI` = quantile(qstForCI.est, 
                                                                                                                                  0.975, na.rm = TRUE)), Va <- c(Va = Va, `Lower bound CI` = quantile(Va.est, 
                                                                                                                                                                                                      0.025, na.rm = TRUE), `Upper bound CI` = quantile(Va.est, 
                                                                                                                                                                                                                                                        0.975, na.rm = TRUE))))
  }
  if (output == "full") {
    if (breeding.design == "half.sib.sire") {
      return(list(QminusF <- c(`Calculated Qst-Fst` = Q.obsMinusF.obs, 
                               `Lower Bound crit. value` = quantile(sim.est, 
                                                                    0.025, na.rm = TRUE), `Upper bound crit. value` = quantile(sim.est, 
                                                                                                                               0.975, na.rm = TRUE)), Distribution <- paste("Qst-Fst values output to file:", 
                                                                                                                                                                            paste(getwd(), "/QminusFvalues_", formatted.time.stamp, 
                                                                                                                                                                                  ".txt", sep = "")), QminusF.p.values <- c(`Lower one-tailed p value` = left.one.tailed.p, 
                                                                                                                                                                                                                            `Upper one-tailed p value` = right.one.tailed.p, 
                                                                                                                                                                                                                            `Two-tailed p value` = two.tailed.p), Fst <- c(`Estimated Fst` = fst.obs, 
                                                                                                                                                                                                                                                                           `Lower Bound CI` = quantile(fst.est, 0.025, na.rm = TRUE), 
                                                                                                                                                                                                                                                                           `Upper bound CI` = quantile(fst.est, 0.975, na.rm = TRUE)), 
                  Fst.resampled <- c(`Fst Resampled` = mean(fst.est, 
                                                            na.rm = TRUE), `Fst std. dev.` = sd(fst.est, 
                                                                                                na.rm = TRUE), `Fst 95% CI` = quantile(fst.est, 
                                                                                                                                       c(0.025, 0.975), na.rm = TRUE), `Fst 99% CI` = quantile(fst.est, 
                                                                                                                                                                                               c(0.005, 0.995), na.rm = TRUE)), Qst <- c(`Estimated Qst` = qst.obs, 
                                                                                                                                                                                                                                         `Lower Bound CI` = quantile(qstForCI.est, 0.025, 
                                                                                                                                                                                                                                                                     na.rm = TRUE), `Upper bound CI` = quantile(qstForCI.est, 
                                                                                                                                                                                                                                                                                                                0.975, na.rm = TRUE)), Qst.neutral <- c(`Qst Resampled` = mean(qst.neut, 
                                                                                                                                                                                                                                                                                                                                                                               na.rm = TRUE), `Qst std. dev.` = sd(qst.neut, 
                                                                                                                                                                                                                                                                                                                                                                                                                   na.rm = TRUE), `Qst 95% crit. values` = quantile(qst.neut, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                    c(0.025, 0.975), na.rm = TRUE), `Qst 99% crit. values` = quantile(qst.neut, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      c(0.005, 0.995), na.rm = TRUE)), ANOVA.table <- c(`MS Pop` = qst.MS[[1]], 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        `MS Sire` = qst.MS[[2]], `MS Dam` = qst.MS[[3]], 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        `MS Within` = qst.MS[[4]], n0primeprime = qst.MS[[5]], 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        n0prime = qst.MS[[6]], n0 = qst.MS[[7]], nc0prime = qst.MS[[8]], 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        nc0 = qst.MS[[9]], ncb0 = qst.MS[[10]], `df Pop` = qst.MS[[11]], 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        `df Sire` = qst.MS[[12]], `df Dam` = qst.MS[[13]], 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        `df Within` = qst.MS[[14]], `Sigma^2 sires` = qst.MS[[15]]), 
                  Va <- c(Va = Va, `Lower bound CI` = quantile(Va.est, 
                                                               0.025, na.rm = TRUE), `Upper bound CI` = quantile(Va.est, 
                                                                                                                 0.975, na.rm = TRUE)), CVa <- c("CVa" <- CVa, 
                                                                                                                                                 `Lower bound CI` = quantile(CVa.est, 0.025, 
                                                                                                                                                                             na.rm = TRUE), `Upper bound CI` = quantile(CVa.est, 
                                                                                                                                                                                                                        0.975, na.rm = TRUE))))
    }
    if (breeding.design == "half.sib.dam") {
      return(list(QminusF <- c(`Calculated Qst-Fst` = Q.obsMinusF.obs, 
                               `Lower Bound crit. value` = quantile(sim.est, 
                                                                    0.025, na.rm = TRUE), `Upper bound crit. value` = quantile(sim.est, 
                                                                                                                               0.975, na.rm = TRUE)), Distribution <- paste("Qst-Fst values output to file:", 
                                                                                                                                                                            paste(getwd(), "/QminusFvalues_", formatted.time.stamp, 
                                                                                                                                                                                  ".txt", sep = "")), QminusF.p.values <- c(`Lower one-tailed p value` = left.one.tailed.p, 
                                                                                                                                                                                                                            `Upper one-tailed p value` = right.one.tailed.p, 
                                                                                                                                                                                                                            `Two-tailed p value` = two.tailed.p), Fst <- c(`Estimated Fst` = fst.obs, 
                                                                                                                                                                                                                                                                           `Lower Bound CI` = quantile(fst.est, 0.025, na.rm = TRUE), 
                                                                                                                                                                                                                                                                           `Upper bound CI` = quantile(fst.est, 0.975, na.rm = TRUE)), 
                  Fst.resampled <- c(`Fst Resampled` = mean(fst.est, 
                                                            na.rm = TRUE), `Fst std. dev.` = sd(fst.est, 
                                                                                                na.rm = TRUE), `Fst 95% CI` = quantile(fst.est, 
                                                                                                                                       c(0.025, 0.975), na.rm = TRUE), `Fst 99% CI` = quantile(fst.est, 
                                                                                                                                                                                               c(0.005, 0.995), na.rm = TRUE)), Qst <- c(`Estimated Qst` = qst.obs, 
                                                                                                                                                                                                                                         `Lower Bound CI` = quantile(qstForCI.est, 0.025, 
                                                                                                                                                                                                                                                                     na.rm = TRUE), `Upper bound CI` = quantile(qstForCI.est, 
                                                                                                                                                                                                                                                                                                                0.975, na.rm = TRUE)), Qst.neutral <- c(`Qst Resampled` = mean(qst.neut, 
                                                                                                                                                                                                                                                                                                                                                                               na.rm = TRUE), `Qst std. dev.` = sd(qst.neut, 
                                                                                                                                                                                                                                                                                                                                                                                                                   na.rm = TRUE), `Qst 95% crit. values` = quantile(qst.neut, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                    c(0.025, 0.975), na.rm = TRUE), `Qst 99% crit. values` = quantile(qst.neut, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      c(0.005, 0.995), na.rm = TRUE)), ANOVA.table <- c(`MS Pop` = qst.MS[[1]], 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        `MS Dam` = qst.MS[[2]], `MS Within` = qst.MS[[3]], 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        n0prime = qst.MS[[4]], n0 = qst.MS[[5]], nb0 = qst.MS[[6]], 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        `df Pop` = qst.MS[[7]], `df Dam` = qst.MS[[8]], 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        `df Within` = qst.MS[[9]], `Sigma^2 dams` = qst.MS[[10]]), 
                  Va <- c(Va = Va, `Lower bound CI` = quantile(Va.est, 
                                                               0.025, na.rm = TRUE), `Upper bound CI` = quantile(Va.est, 
                                                                                                                 0.975, na.rm = TRUE)), CVa <- c(CVa = CVa, 
                                                                                                                                                 `Lower bound CI` = quantile(CVa.est, 0.025, 
                                                                                                                                                                             na.rm = TRUE), `Upper bound CI` = quantile(CVa.est, 
                                                                                                                                                                                                                        0.975, na.rm = TRUE))))
    }
  }
}