tr <-
  function(x)
    sum(diag(as.matrix(x)))


norm <-
  function(x,
           type = c("none", "standard", "center", "range", "SVD")) {
    type <- match.arg(type,
                      choices = eval(formals(norm)$type),
                      several.ok = FALSE)
    x <- as.matrix(x)
    switch(
      type,
      "none" = x,
      "standard" = scale(x, center = TRUE, scale = TRUE),
      "center"   = scale(x, center = TRUE, scale = FALSE),
      "range"    = apply(x, 2, function(x)
        (x - min(x)) / (max(x) - min(x))),
      "SVD"      = {
        x <- scale(x, center = TRUE, scale = TRUE)
        p <- ncol(x)
        SVD <- svd(x, nu = 0)
        x %*% SVD$v %*% diag(1 / sqrt(SVD$d), p, p)
      }
    )
  }


#' Random partition of objects into classes
#'
#' @description Performs a random partition of objects into classes.
#' @param n.obs Number of objects
#' @param G Number of classes
#'
#' @return A binary and row-stochastic matrix with \eqn{n.obs} rows and \eqn{G} columns.
#' @details No empty classes can occur.
#'
#' @examples
#' rand.member(10, 3)
#' @export
#' rand.member
rand.member <-
  function(n.obs,
           G) {
    if (G > n.obs)
      stop('The number of groups is larger than the number of observations')
    if (G == n.obs) {
      U <- diag(n.obs)
    } else if (G == 1) {
      U <- c(rep(1, n.obs))
    }
    else {
      U <- matrix(0, n.obs, G)
      U[1:G,] = diag(G)
      U[(G + 1):n.obs, 1] <- 1
      for (i in (G + 1):n.obs) {
        U[i,] <- U[i, sample(G)]
      }
      U <- U[sample(n.obs),]
    }
    return(as.matrix(U))
  }


UCMSigmaV <-
  function(X,
           S,
           m,
           rndst,
           maxiter = 100,
           eps = sqrt(.Machine$double.eps),
           model.name = c("EUUU", "EUUE", "EUEE", "EEEU", "EEEE", "EEEF", "EEFF", "EFFF", "FIII", "FIIF", "FIFF", "FFFI", "FFFF")) {
    model.name <- match.arg(model.name)
    p <- dim(S)[1]
    Im <- diag(m)
    if (m == p) {
      Vopt <- diag(p)
    } else if (m == 1) {
      Vopt <- as.matrix(rep(1, p))
    } else {
      for (loop in 1:rndst) {
        it <- 0
        if (m > p - 2){
          V <- mclust::unmap(stats::kmeans(t(X), centers = m, iter.max = 100)$cluster)
        } else {
        V <- mclust::unmap(ClusterR::KMeans_rcpp(t(X), clusters = m)$clusters)}
        Sv <- estimate_Sv(S, V, model.name)
        Sw <- estimate_Sw(S, Sv, V, model.name)
        Sb <- estimate_Sb(S, V, model.name)
        if (m > 1) {
          Sw <- check_constraint_SwSb(Sw, Sb, model.name)$Sw
        }
        if (m < p) {
          Sv <- check_constraint_SvSw(Sv, Sw, model.name)$Sv
        }
        St.psd <- psd_sigma(V %*% (Sw + Sb) %*% t(V) - diag(diag(V %*% Sw %*% t(V))) + diag(diag(V %*% Sv %*% t(V))))
        St <- St.psd$S
        if (St.psd$constr.psd != 0) {
          Sv <- estimate_Sv(St, V, model.name)
          Sw <- estimate_Sw(St, Sv, V, model.name)
          Sb <- estimate_Sb(St, V, model.name)
          if (m > 1) {
            Sw <- check_constraint_SwSb(Sw, Sb, model.name)$Sw
            if (m == p) {
              Sv <- Sw
            }
          }
          if (m < p) {
            Sv <- check_constraint_SvSw(Sv, Sw, model.name)$Sv
          }
          St <- V %*% (Sw + Sb) %*% t(V) - diag(diag(V %*% Sw %*% t(V))) + diag(diag(V %*% Sv %*% t(V)))
        }
        while (min(eigen(St, symmetric = TRUE, only.values = TRUE)$values) <= 0 || det(St) <= sqrt(.Machine$double.eps)) {
          a <- abs(min(eigen(St, symmetric = TRUE, only.values = TRUE)$values)) + sqrt(.Machine$double.eps)
          St <- St + a * diag(p)
          Sv <- Sv + a * diag(m)
        }
        fo <- tr(MASS::ginv(St,  tol = .Machine$double.xmin) %*% S) + log(det(St))
        fmin <- fo
        if (loop == 1) {
          Vopt <- V
          fopt <- fmin
        }
        fdif <- fmin
        while (fdif > eps && it <= maxiter) {
          it <- it + 1
          fo <- fmin
          for (i in 1:p) {
            posmin <- which(V[i, ] == 1)
            V[i,] <- Im[posmin, ]
            for (j in 1:m) {
              V[i,] <- Im[j, ]
              if (sum(V[, posmin]) > 0) {
                Sv <- estimate_Sv(S, V, model.name)
                Sw <- estimate_Sw(S, Sv, V, model.name)
                Sb <- estimate_Sb(S, V, model.name)
                if (m > 1) {
                  Sw <- check_constraint_SwSb(Sw, Sb, model.name)$Sw
                }
                if (m < p) {
                  Sv <- check_constraint_SvSw(Sv, Sw, model.name)$Sv
                }
                St.psd <- psd_sigma(V %*% (Sw + Sb) %*% t(V) - diag(diag(V %*% Sw %*% t(V))) + diag(diag(V %*% Sv %*% t(V))))
                St <- St.psd$S
                if (St.psd$constr.psd != 0) {
                  Sv <- estimate_Sv(S, V, model.name)
                  Sw <- estimate_Sw(S, Sv, V, model.name)
                  Sb <- estimate_Sb(S, V, model.name)
                  if (m > 1) {
                    Sw <- check_constraint_SwSb(Sw, Sb, model.name)$Sw
                    if (m == p) {
                      Sv <- Sw
                    }
                  }
                  if (m < p) {
                    Sv <- check_constraint_SvSw(Sv, Sw, model.name)$Sv
                  }
                  St <- V %*% (Sw + Sb) %*% t(V) - diag(diag(V %*% Sw %*% t(V))) + diag(diag(V %*% Sv %*% t(V)))
                }
                while (min(eigen(St, symmetric = TRUE, only.values = TRUE)$values) <= 0 || det(St) <= sqrt(.Machine$double.eps)) {
                  a <- abs(min(eigen(St, symmetric = TRUE, only.values = TRUE)$values)) + sqrt(.Machine$double.eps)
                  St <- St + a * diag(p)
                  Sv <- Sv + a * diag(m)
                }
                ff <- tr(MASS::ginv(St,  tol = .Machine$double.xmin) %*% S) + log(det(St))
                if (ff < fmin) {
                  fmin <- ff
                  posmin <- j
                }
              }
            }
            V[i,] <- Im[posmin, ]
          }
          fdif <- fo - fmin
          if (fmin < fopt) {
            Vopt <- V
            fopt <- fmin
          }
        }
      }
    }
    return(Vopt)
  }


ultrcov <-
  function(Sb,
           V) {
    m <- nrow(Sb)
    A <- matrix(0, m - 1, 1)
    B <- matrix(0, m - 1, 1)
    levfus <- matrix(0, m - 1, 1)
    PP <- matrix(0, m, 1)
    P <- matrix(0, m, 1)
    class <- matrix(0, m, 1)
    SU <- matrix(0, m, m)
    for (q in 1:m) {
      PP[q, 1] <- q
      class[q, 1] <- sum(V[, q])
    }
    for (ustep in 1:(m - 1)) {
      rmax <- (-.Machine$double.xmax)
      for (i in 1:(m - 1)) {
        if (PP[i, 1] == i) {
          for (j in (i + 1):m) {
            if (PP[j, 1] == j) {
              if (Sb[i, j] >= rmax) {
                ic <- i
                jc <- j
                rmax <- Sb[i, j]
              }
            }
          }
        }
      }
      for (j in 1:m) {
        if (PP[j, 1] == jc) {
          PP[j, 1] <- ic
        }
      }
      A[ustep, 1] <- ic
      B[ustep, 1] <- jc
      levfus[ustep, 1] <- rmax
      for (i in 1:m) {
        if (i != ic && PP[i, 1] == i) {
          rs <-
            (class[ic, 1] * Sb[ic, i] + class[jc, 1] * Sb[jc, i]) / (class[ic, 1] + class[jc, 1])
          Sb[ic, i] <- rs
          Sb[i, ic] <- rs
        }
      }
      class[ic, 1] <- class[ic, 1] + class[jc, 1]
    }
    for (i in 1:m) {
      P[i, 1] <- i
    }
    for (k in 1:(m - 1)) {
      for (i in A[k, 1]:m) {
        if (P[i, 1] == A[k, 1]) {
          for (j in B[k, 1]:m) {
            if (P[j, 1] == B[k, 1]) {
              SU[i, j] <- levfus[k, 1]
              SU[j, i] <- levfus[k, 1]
            }
          }
        }
      }
      for (i in A[k, 1]:m) {
        if (P[i, 1] == B[k, 1]) {
          P[i, 1] <- A[k, 1]
        }
      }
    }
    return(SU)
  }

prior <-
  function(w) {
    pp <- colSums(w) / dim(w)[1]
    return(pp)
  }

post <-
  function(X,
           pp,
           mu,
           Sigma) {
    G <- dim(mu)[1]
    w <- matrix(as.double(NA), dim(X)[1], G)
    for (g in 1:G) {
      w[,g] = log(pp[g]) + mclust::dmvnorm(X, mu[g,], as.matrix(Sigma[[g]]), log = TRUE)
    }
    wnorm <- apply(w, 1, max)
    w <- exp(sweep(w, 1, wnorm, "-"))
    w <- w / rowSums(w)
    w[which(w < sqrt(.Machine$double.eps), arr.ind = TRUE)] <- sqrt(.Machine$double.eps)
    return(w)
  }


component_param <-
  function(X,
           w) {
    n.obs <- dim(X)[1]
    p <- dim(X)[2]
    G <- dim(w)[2]
    Sigma.new <- list()

    mu.new <- sweep(t(w) %*% X, 1, 1 / colSums(w), "*")

    r <- sqrt(w)
    for (g in 1:G) {
      Xo <- sweep(X, 2, mu.new[g,], "-", check.margin = FALSE)
      Xo <- sweep(Xo, 1, r[, g], "*", check.margin = FALSE)
      Sigma.new[[g]] <-  (t(Xo) %*% Xo) / colSums(w)[g]
    }
    return(list(mu = mu.new,
                Sigma = Sigma.new))
  }


loglikH_EI_G_pugmm <-
  function(X,
           pp,
           mu,
           Sigma,
           Sv,
           Sw,
           Sb,
           V,
           w,
           gaussian) {
    n <- dim(X)[1]
    p <- dim(X)[2]
    G <- length(pp)
    llh <- rep(as.double(NA), G)
    lf <- matrix(as.double(NA), n, G)
    for (g in 1:G) {
      nv <- colSums(V[[g]])
      if (any(nv == 1)) {
        gaussian = "mclust"
      }
      if (gaussian == "mclust") {
        lf[, g] <-
          mclust::dmvnorm(X, mu[g, , drop = FALSE], Sigma[[g]], log = TRUE) + log(pp[g])
      } else {
        lf[, g] <-
          canonical_norm(X, mu[g, , drop = FALSE], Sv[[g]], Sw[[g]], Sb[[g]], V[[g]]) + log(pp[g])
      }
    }
    llh <- lf * w - w * log(w)
    loglik <- sum(llh)
    return(loglik)
  }

loglikH_F_pugmm <-
  function(X,
           g,
           pp,
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
           w,
           gaussian) {
    n.obs <- dim(X)[1]
    G <- length(pp)
    llh <- rep(as.double(NA), G)
    lf <- matrix(as.double(NA), n.obs, G)
    nv0 <- colSums(V.current)
    for (k in 1:G) {
      nv <- colSums(V.other[[k]])
      if (any(nv0 == 1) || any(nv == 1)) {
        gaussian = "mclust"
      }
      if (gaussian == "mclust") {
        if (k == g) {
          lf[, k] <-
            as.matrix(mclust::dmvnorm(X, mu[k, , drop = FALSE], Sigma.current, log = TRUE)) + log(pp[k])
        } else {
          lf[, k] <-
            as.matrix(mclust::dmvnorm(X, mu[k, , drop = FALSE], Sigma.other[[k]], log = TRUE)) + log(pp[k])
        }
      } else {
        for (k in 1:G) {
          if (k == g) {
            lf[, k] <-
              canonical_norm(X, mu[k, , drop = FALSE], Sv.current, Sw.current, Sb.current, V.current) + log(pp[k])
          } else {
            lf[, k] <-
              canonical_norm(X, mu[k, , drop = FALSE], Sv.other[[k]], Sw.other[[k]], Sb.other[[k]], V.other[[k]]) + log(pp[k])
          }
        }
      }
    }
    llh <- lf * w - w * log(w)
    loglik <- sum(llh)
    return(loglik)
  }

#' PUGMM Model Names
#'
#' @description
#' Description  of the model names used in the \emph{PUGMM} package.
#'
#' @return Available models in PUGMM, i.e., the thirteen extended ultrametric covariance structures of PUGMM.
#'
#' @details
#' The PUGMM model names in the \emph{PUGMM} package are characterized by four letters:
#' \itemize{
#' \item First letter: it refers to the variable-group membership matrix \eqn{V}, which can be equal (E) or free to vary (F) across components.
#' \item Second, third, fourth letters: they refer to the matrices of the group variances \eqn{\Sigma_{V}}, the within-group covariances \eqn{\Sigma_{W}} and the between-group covariances \eqn{\Sigma_{B}}, respectively, by indicating if they are unique (U, i.e., equal within and across components), isotropic (I, i.e., equal within components), equal (E, i.e., equal across components) or free to vary across components (F).
#' }
#' @references Cavicchia, C., Vichi, M., Zaccaria, G. (2024) Parsimonious ultrametric Gaussian mixture models. \emph{Statistics and Computing}, 34, 108.
#' @seealso [pugmm()]
#'
#' @examples
#' pugmm_available_models()
#'
#' @export pugmm_available_models
#' @export
pugmm_available_models <-
  function() {
    available.models <- c("EUUU",
                          "EUUE",
                          "EUEE",
                          "EEEU",
                          "EEEE",
                          "EEEF",
                          "EEFF",
                          "EFFF",
                          "FIII",
                          "FIIF",
                          "FIFF",
                          "FFFI",
                          "FFFF")
    return(available.models)
  }


number_param <-
  function(p,
           G,
           m,
           model.name = c("EUUU", "EUUE", "EUEE", "EEEU", "EEEE", "EEEF", "EEFF", "EFFF", "FIII", "FIIF", "FIFF", "FFFI", "FFFF")) {
    model.name <- match.arg(model.name)
    switch(
      model.name,
      "EUUU" = {pm <- G + (G + 1) * p + 2},
      "EUUE" = {pm <- G + (G + 1) * p + m},
      "EUEE" = {pm <- G + (G + 1) * p + 2* m - 1},
      "EEEU" = {pm <- G + ((G + 1) * p) + 2 * m},
      "EEEE" = {pm <- G + (G + 1) * p + 3 * m - 2},
      "EEEF" = {pm <- (G + 1) * p + (G + 2) * m - 1},
      "EEFF" = {pm <- (G + 1) * p + (2 * G  + 1) * m - 1},
      "EFFF" = {pm <- (G + 1) * p + 3 * G * m - 1},
      "FIII" = {pm <- 2 * G * (p + 2) - 1},
      "FIIF" = {pm <- G * (2 * p + m + 2) - 1},
      "FIFF" = {pm <- G * (2 * p + 2* m + 1) - 1},
      "FFFI" = {pm <- 2 * G * (p + m + 1) - 1},
      "FFFF" = {pm <- G * (2 * p + 3 * m) - 1}
    )
    return(pm)
  }

ver_ultrametric <-
  function(U) {
    n.obs <- dim(U)[1]
    value <- 0
    for (i in 1:(n.obs - 2)) {
      for (j in (i + 1):(n.obs - 1)) {
        for (h in (j + 1):n.obs) {
          if (U[i, j] >= U[i, h] &&
              U[i, j] >= U[j, h]) {
            value <- value + ((U[i, h] - U[j, h]) ^ 2)
          }
          else if (U[i, h] >= U[i, j] &&
                   U[i, h] >= U[j, h]) {
            value <- value + (U[i, j] - U[j, h]) ^ 2
          }
          else {
            value <- value + (U[i, h] - U[i, j]) ^ 2
          }
        }
      }
    }
    return(value)
  }

canonical_norm <-
  function(X,
           mu,
           Sv,
           Sw,
           Sb,
           V) {
    p <- dim(V)[1]
    m <- dim(V)[2]
    if (m == 1) {
      Xo <- X
      muo <- mu
      nv <- p
      A <- Sv + (nv - 1) * Sw
    } else {
    cmax <- apply(V, 2, which.max)
    omax <- order(cmax)
    Vo <- V[, omax]
    v <- mclust::map(Vo)
    vs <- sort(v)
    ov <- order(v)
    Sv <- Sv[omax, omax]
    Sw <- Sw[omax, omax]
    Sb <- Sb[omax, omax]
    Vo <- mclust::unmap(vs)
    nv <- colSums(Vo)
    A <- Sv + diag(nv - 1) * Sw + Sb * sqrt(nv %*% t(nv))
    Xo <- X[, ov]
    muo <- mu[ov]
    }
    Xs <- sweep(Xo , 2, muo, "-", check.margin = FALSE)
    D <- A
    Qv <- c(rep(1 / sqrt(nv[1]), nv[1]))
    Qvt <- mcompanion::null_complement(Qv)
    for (q in 1:m) {
      DD <- (Sv[q, q] - Sw[q, q]) * diag(nv[q] - 1)
      D <- Matrix::bdiag(D, DD)
      if (q < m) {
        QQ <- c(rep(1 / sqrt(nv[q + 1]), nv[q + 1]))
        Qv <- Matrix::bdiag(Qv, QQ)
        Qvt <- Matrix::bdiag(Qvt, mcompanion::null_complement(QQ))
      }
    }
    Q <- as.matrix(cbind(Qv, Qvt))
    d <- matrix(D[(m + 1):p, (m + 1):p])
    detD <- det(A) * prod(d[which(d != 0)])
    invD <- Matrix::bdiag(MASS::ginv(A,  tol = .Machine$double.xmin), diag(1 / d[which(d != 0)]))
    lf <-
      diag(as.matrix((-1 / 2) * (
        p * log(2 * pi) + log(detD) + Xs %*% Q %*% invD %*% t(Q) %*% t(Xs)
      )))
    return(lf)
  }



