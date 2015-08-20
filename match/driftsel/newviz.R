newviz <- function (fixedpost, popefpost, Gpost, THpost, traits, siz = 0.5, 
          main = NA, xlab, ylab, xlim, ylim) 
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
 # xlab = paste("trait", traits[i])
 # ylab = paste("trait", traits[j])
  M = max(c(sqrt(2 * G[i, i] * diag(TH)), sqrt(2 * G[j, j] * 
                                                 diag(TH)), popef))
  plot(ei, ej, pch = 16, cex = 0, xlim = xlim, 
       ylim = ylim, xlab = xlab, ylab = ylab, main = main)
  for (k in 1:npop) {
    mu = c(ei, ej)
    Sigma = 2 * TH[k, k] * Gthis
    if (k > 9) {
      k = 1 + k%%8
    }
    ellipsis(mu, Sigma, siz, 1, colvec[k])
  }
  text(ei + popef[, i], ej + popef[, j], 1:npop, cex = 1, col = colvec)
  text(ei, ej, "A")
}
<environment: namespace:driftsel>