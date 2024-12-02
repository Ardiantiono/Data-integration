elpd <- function(samples.fit, samples.pred, access.pred, tree.pred, 
                 canopy.pred, my.iter, n.iter, I, type) {
  
  param.fit.names <- attr(samples.fit, 'dimnames')[[2]]
  curr.indx <- which(substr(param.fit.names, 1, 7) %in% c('beta.0.')) #7 = #letters
  beta.0.fit <- samples.fit[my.iter, curr.indx] #error during logit
  #perhaps skip logit (cause beta.0 still log like beta.1)
  #dim(samples.fit)
  curr.indx <- which(substr(param.fit.names, 1, 7) %in% c('beta.1.')) 
  beta.1.fit <- samples.fit[my.iter, curr.indx]
  curr.indx <- which(substr(param.fit.names, 1, 7) %in% c('beta.2.'))
  beta.2.fit <- samples.fit[my.iter, curr.indx]
  curr.indx <- which(substr(param.fit.names, 1, 7) %in% c('beta.3.'))
  beta.3.fit <- samples.fit[my.iter, curr.indx]
  curr.indx <- which(substr(param.fit.names, 1, 7) %in% c('beta.4.'))
  beta.4.fit <- samples.fit[my.iter, curr.indx]
  curr.indx <- which(substr(param.fit.names, 1, 7) %in% c('beta.5.'))
  beta.5.fit <- samples.fit[my.iter, curr.indx]
  
  # Predicted values
  param.pred.names <- attr(samples.pred, 'dimnames')[[2]]
  val <- ifelse(type == 'sm', 'z.sm', ifelse(type == 'ct', 'z.ct', 'z.tr'))
  #val <- 'z.ct'
  curr.indx <- which(substr(param.pred.names, 1, 4) %in% val)
  z.pred.samples <- samples.pred[my.iter, curr.indx]
  J.pred <- ncol(z.pred.samples) / I 
  z.pred.samples <- array(z.pred.samples, dim = c(n.iter, I, J.pred))
  z.fit.samples <- array(NA, dim(z.pred.samples))
  pred.dens.hbef <- array(NA, dim (z.pred.samples))
  # Composition sampling algorithm to get cross validation metrics ----------
  for (a in 1:n.iter) {
    print(paste("Currently on iteration ", a, " out of ", n.iter, sep = ''))
    for (i in 1:I) {
      for (j in 1:J.pred) {
            psi <- logit.inv(beta.0.fit[a, i] + 
                               beta.1.fit[a, i] * access.pred[j] + 
                               beta.2.fit[a, i] * access.pred[j]^2 + 
                               beta.3.fit[a, i] * tree.pred[j] + 
                               beta.4.fit[a, i] * tree.pred[j]^2 + 
                               beta.5.fit[a, i] * canopy.pred[j])
            z.fit.samples[a, i, j] <- rbinom(1, 1, psi)
            pred.dens.hbef[a, i, j] <- dbinom(z.pred.samples[a, i, j], 1, psi)
      } # j
    } # i
  }
  # Compute ELPD estimate
  return(apply(pred.dens.hbef, c(2, 3), mean))
}
