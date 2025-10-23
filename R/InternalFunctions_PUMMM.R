post_puMmm <-
  function(X,
           Y,
           pp,
           lambda,
           mu,
           Sigma) {
    G <- dim(mu)[1]
    w <- matrix(as.double(NA), dim(X)[1], G)
    for (g in 1:G) {
      w[,g] = log(pp[g]) + X%*%lambda[,g] + mclust::dmvnorm(Y[,,g], mu[g,], as.matrix(Sigma[[g]]), log = TRUE)
    }
    wnorm <- apply(w, 1, max)
    w <- exp(sweep(w, 1, wnorm, "-"))
    w <- w / rowSums(w)
    w[which(w < sqrt(.Machine$double.eps), arr.ind = TRUE)] <- sqrt(.Machine$double.eps)
    return(w)
  }

component_param_PUMMM <-
  function(Y,
           w) {
    n.obs <- dim(Y)[1]
    p <- dim(Y)[2]
    G <- dim(w)[2]
    Sigma.new <- list()
    
    mu.new = array(0, dim=c(G,p))
    for (g in 1:G){
      mu.new[g,] = colSums(w[,g]*Y[,,g])/colSums(w)[g]
    }
    
    r <- sqrt(w)
    for (g in 1:G) {
      Yo <- sweep(Y[,,g], 2, mu.new[g,], "-", check.margin = FALSE)
      Yo <- sweep(Yo, 1, r[, g], "*", check.margin = FALSE)
      Sigma.new[[g]] <-  (t(Yo) %*% Yo) / colSums(w)[g]
    }
    
    return(list(mu = mu.new,
                Sigma = Sigma.new))
  }

lambda_update <-
  function(X,
           w,
           lambda,
           mu,
           sigma) {
    G = dim(w)[2]
    n = dim(X)[1]
    p = dim(X)[2]
    
    index = matrix(rep(0, p*G), nrow=p, ncol=G)
    for (g in 1:G){
      for (i in 1:p){
        if (abs(lambda[i,g]) < 1e-12){
          index[i,g] = 0
          lambda[i,g] = 0
        } else {
          index[i,g] = 1
          lambda[i,g] = lambda[i,g]
        }
      }
    }
    
    lambda.new = matrix(0, nrow = p, ncol = G)
    for (g in 1:G){
      lambda.new[,g] <- simplex(Q,
                                n1 = n,
                                p = p,
                                index_vec = index[,g],
                                X_mat = as.matrix(X),
                                z_vec = w[,g],
                                mu_vec = mu[g,],
                                sigma_mat = sigma[[g]],
                                start_vec = lambda[,g],
                                EPSILON = 1e-06,
                                scale = 0.1)
    }
    return(lambda.new)
  }

manly_transform_mat <- function(X, lambda){
  #lambda is a matrix
  n = dim(X)[1]
  p = dim(X)[2]
  G = dim(lambda)[2]
  Y = array(NA, dim=c(n,p,G))
  for (g in 1:G){
    for (i in 1:p){
      if (abs(lambda[i,g]) < 1e-12){
        Y[,i,g] = X[,i]
      } else {
        Y[,i,g] = (exp(X[,i]*lambda[i,g]) - 1)/lambda[i,g]
      }
    }
  }
  return(Y)
}

manly_transform <-
  function(X,
           lambda) {
    #lambda is a vector
    n = dim(X)[1]
    p = dim(X)[2]
    Y = array(NA, dim=c(n,p))
    for (j in 1:p){
      for (i in 1:n){
        if (abs(lambda[j]) < 1e-12){
          Y[i,j] = X[i,j]
        } else {
          Y[i,j] = (exp(X[i,j]*lambda[j]) - 1)/lambda[j]
        }
      }
    }
    return(Y)
  }


Q = function(lambda, index, X, w, mu, sigma){
  return(-sum(w*(mclust::dmvnorm(manly_transform(X,lambda), mu, sigma, log = TRUE) + X%*%lambda)))
}

loglikH_EI_G_puMmm <-
  function(X,
           Y,
           pp,
           lambda,
           mu,
           Sigma,
           Sv,
           Sw,
           Sb,
           V,
           w) {
    n <- dim(X)[1]
    p <- dim(X)[2]
    G <- length(pp)
    lf <- matrix(as.double(NA), n, G)
    for (g in 1:G) {
      lf[, g] <- mclust::dmvnorm(Y[,,g], mu[g, , drop = FALSE], Sigma[[g]], log = TRUE) + X%*%lambda[,g]  + log(pp[g])
    }
    loglik = sum(lf * w - w * log(w))
    return(loglik)
  }

loglikH_F_puMmm <-
  function(X,
           Y,
           g,
           pp,
           lambda,
           mu,
           Sigma.current,
           Sigma.other,
           Sv.current,
           Sv.other,
           Sw.current,
           Sw.other,
           Sb.current,
           Sb.other,
           V.current,
           V.other,
           w) {
    n.obs <- dim(X)[1]
    G <- length(pp)
    lf <- matrix(as.double(NA), n.obs, G)
    for (k in 1:G) {
      if (k == g) {
        lf[, k] <-
          as.matrix(mclust::dmvnorm(Y[,,k], mu[k, , drop = FALSE], Sigma.current, log = TRUE)) + X%*%lambda[,k] + log(pp[k])
      } else {
        lf[, k] <-
          as.matrix(mclust::dmvnorm(Y[,,k], mu[k, , drop = FALSE], Sigma.other[[k]], log = TRUE)) + X%*%lambda[,k] + log(pp[k])
      }
    }
    loglik = sum(lf * w - w * log(w))
    return(loglik)
  }

number_param_puMmm <-
  function(p,
           G,
           m,
           model.name = c("EUUU", "EUUE", "EUEE", "EEEU", "EEEE", "EEEF", "EEFF", "EFFF", "FIII", "FIIF", "FIFF", "FFFI", "FFFF")) {
    model.name <- match.arg(model.name)
    switch(
      model.name,
      "EUUU" = {pm <- G + (G + 1) * p + 2 + G*p},
      "EUUE" = {pm <- G + (G + 1) * p + m + G*p},
      "EUEE" = {pm <- G + (G + 1) * p + 2* m - 1 + G*p},
      "EEEU" = {pm <- G + ((G + 1) * p) + 2 * m + G*p},
      "EEEE" = {pm <- G + (G + 1) * p + 3 * m - 2 + G*p},
      "EEEF" = {pm <- (G + 1) * p + (G + 2) * m - 1 + G*p},
      "EEFF" = {pm <- (G + 1) * p + (2 * G  + 1) * m - 1 + G*p},
      "EFFF" = {pm <- (G + 1) * p + 3 * G * m - 1 + G*p},
      "FIII" = {pm <- 2 * G * (p + 2) - 1 + G*p},
      "FIIF" = {pm <- G * (2 * p + m + 2) - 1 + G*p},
      "FIFF" = {pm <- G * (2 * p + 2* m + 1) - 1 + G*p},
      "FFFI" = {pm <- 2 * G * (p + m + 1) - 1 + G*p},
      "FFFF" = {pm <- G * (2 * p + 3 * m) - 1 + G*p}
    )
    return(pm)
  }
