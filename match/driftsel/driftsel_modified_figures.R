ellipsis <-
  function(mu, Sig, prob.mass, lwd, col){
    radius <- sqrt(qchisq(prob.mass, 2))
    theta <- 2*pi*(0:360) / 360
    unit.circle <- cbind(cos(theta), sin(theta))
    Q <- base::chol(Sig, pivot=TRUE)
    order <- order(attr(Q, "pivot"))
    ellipse <- t(mu + radius * t(unit.circle %*% Q[, order]))
    lines(ellipse[,1], ellipse[,2], lwd=lwd, col=col)
  } 


viz.traits2 <- function (fixedpost, popefpost, Gpost, THpost, traits, siz = 0.5, 
          main = NA, xlab = NA, ylab = NA, plot.id = NA) 
{
  fixedpost = fixedpost[1, traits, ]
  popefpost = popefpost[, traits, ]
  Gpost = Gpost[traits, traits, ]
  popef = matrix(NA, nrow = dim(popefpost)[1], ncol = dim(popefpost)[2])
  G = matrix(NA, nrow = dim(Gpost)[1], ncol = dim(Gpost)[2])
  TH = matrix(NA, nrow = dim(THpost)[1], ncol = dim(THpost)[2])
  mu = rep(NA, dim(Gpost)[1])
  for (i in 1:length(traits)) {
    mu[i] = mean(fixedpost[i, ])
    for (j in 1:length(traits)) {
      G[i, j] = mean(Gpost[i, j, ])
    }
    for (j in 1:dim(popefpost)[1]) {
      popef[j, i] = mean(popefpost[j, i, ])
    }
  }
  for (i in 1:nrow(TH)) {
    for (j in 1:ncol(TH)) {
      TH[i, j] = mean(THpost[i, j, ])
    }
  }
  npop = ncol(TH)
  ntr = length(traits)
  npop = nrow(popef)
  colvec = 1:9
  colvec[7] = "orange"
  colvec[9] = "chartreuse"
  i = 1
  j = 2
  ei = mu[i]
  ej = mu[j]
  Gthis = matrix(c(G[i, i], G[i, j], G[j, i], G[j, j]), ncol = 2)
#   xlab = paste("trait", traits[i])
#   ylab = paste("trait", traits[j])
  M = max(c(sqrt(2 * G[i, i] * diag(TH)), sqrt(2 * G[j, j] * 
                                                 diag(TH)), popef))
  plot(ei, ej, pch = 16, cex = 0, xlim = c(ei - M, ei + M), 
       ylim = c(ej - M, ej + M), xlab = xlab, ylab = ylab, main = main)
  for (k in 1:npop) {
    mu = c(ei, ej)
    Sigma = 2 * TH[k, k] * Gthis
    if (k > 9) {
      k = 1 + k%%8
    }
    ellipsis(mu, Sigma, siz, 1, colvec[k])
  }
  text(ei + popef[, i], ej + popef[, j], 1:npop, cex = 1, col = colvec)
  text(ei, ej, "A", cex=1.5)
#   mtext(plot.num, side=3, line=5 )
}

viz.theta2 <- function (thetapost, distance = T, center = F, main = NA) 
{
  npop = dim(thetapost)[1]
  theta = D = matrix(NA, npop, npop)
  for (i in 1:npop) {
    for (j in 1:i) {
      theta[i, j] = theta[j, i] = mean(thetapost[i, j, 
                                                 ])
    }
  }
  for (i in 1:npop) {
    for (j in 1:i) {
      D2 = 2 * (theta[i, i] + theta[j, j] - 2 * theta[i, 
                                                      j])
      D[i, j] = D[j, i] = sqrt(D2)
    }
  }
  D = rbind(D, sqrt(2 * diag(theta)))
  D = cbind(D, c(sqrt(2 * diag(theta)), 0))
  D = sqrt(2/pi) * D
  fit = cmdscale(D, k = 2)
  if (center) {
    fit[nrow(fit), ] = c(0, 0)
  }
  if (is.na(main)) {
#     main = "Expected drift distances (units of ancestral SD)"
    if (!is.na(distance)) {
      if (!distance) {
        main = "Population-level coancestry"
      }
    }
  }
  eps = max(abs(fit)) * 0.1
  M = max(abs(fit)) + eps
  cx1 = 1 - 1 * is.na(distance)
  plot(fit[, 1], fit[, 2], ylab = "MDS 2", main = main, pch = 16, 
       xlim = c(-M, M), ylim = c(-M, M), axes = T, xlab = "MDS 1", 
       cex = cx1)
  eps = max(abs(fit)) * 0.1
  colvec = 1:9
  colvec[7] = "orange"
  colvec[9] = "chartreuse"
  eps = (1 - 1 * is.na(distance)) * eps
  text(fit[1:npop, 1] + eps, fit[1:npop, 2], 1:npop, col = colvec)
  text(fit[npop + 1, 1] + eps, fit[npop + 1, 2], "A", cex=1.5)
  if (!is.na(distance)) {
    for (i in 1:(npop + 1)) {
      for (j in 1:(npop + 1)) {
        lines(c(fit[i, 1], fit[j, 1]), c(fit[i, 2], fit[j, 
                                                        2]))
      }
    }
    eps = 0.8 * eps
    if (distance) {
      dd = round(D, 2)
    }
    else {
      self = diag(theta)
      dd = rbind(theta, self)
      dd = cbind(dd, c(self, 0))
      dd = round(dd, 2)
    }
    for (i in 1:nrow(fit)) {
      if (i > 1) {
        for (j in 1:(i - 1)) {
          text(0.5 * (fit[i, 1] + fit[j, 1]), 0.5 * (fit[i, 
                                                         2] + fit[j, 2]) + eps, dd[i, j], cex = 0.75)
        }
      }
    }
  }
}