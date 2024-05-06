#' Ultrametric Correlation Matrix
#' @description
#' Fit an ultrametric correlation matrix on a nonnegative correlation one.
#' @param R (\eqn{p \times p}) nonnegative correlation matrix.
#' @param m Integer specifying the number of variable groups.
#' @param rndstart Integer value specifying the number of random starts.
#' @param maxiter Integer value specifying the maximum number of iterations of the EM algorithm (default: maxiter = 100).
#' @param eps Numeric value specifying the tolerance for the convergence criterion used in the coordinate descent algorithm (default: eps = 1e-6).
#'
#' @return A list with the following elements: \cr
#' @return `call` Matched call.
#' @return `V` Optimal binary and row-stochastic (\eqn{p \times m}) variable-group membership matrix.
#' @return `Rt` Optimal (\eqn{p \times p}) ultrametric correlation matrix.
#' @return `Rw` Optimal (\eqn{m \times m}) within-concept consistency (diagonal) matrix.
#' @return `Rb` Optimal (\eqn{m \times m}) between-concept correlation matrix.
#' @return `of` Objective function corresponding to the optimal solution.
#' @return `loop` Random start corresponding to the optimal solution.
#' @return `iter` Number of iterations needed to obtain the optimal solution.
#' @references Cavicchia, C., Vichi, M., Zaccaria, G. (2020) The ultrametric correlation matrix for modelling hierarchical latent concepts. \emph{Advances in Data Analysis and Classification}, 14(4), 837-853.
#' @examples
#' data(penguins)
#' R <- cor(penguins[, 2:5])
#' UCM(R, 4, 1)
#' @export
#' UCM
UCM <- function(R,
                m,
                rndstart,
                maxiter = 100,
                eps = 1e-6){
  call <- match.call()
  p <- dim(R)[1]
  pp <- (1:p)
  Im <- diag(m)
  ts <- sum(R*R)
  fmm <- matrix(0, rndstart, 1)
  fOtt <- Inf
  if (m == 1 || m == p) {
    rndstart <- 1
    maxiter <- -1
  }
  for (loop in 1:rndstart){
    it <- 0
    V <- rand.member(p, m)
    if (m < p) {
      Rw <- t(V) %*% (R-diag(p)) %*% V %*% MASS::ginv((t(V) %*% V)^2 - (t(V) %*% V),  tol = .Machine$double.xmin)
      if (m > 1) {
        Rw <- diag(diag(Rw))
      }
      if (any(colSums(V) == 1)){
        Rw[which(colSums(V) == 1), which(colSums(V) == 1)] <- diag(length(which(colSums(V) == 1)))
      }
    } else if (m == p) {
      Rw <- diag(m)
    }
    Rw_wo <- Rw
    if (m > 1) {
      su <- t(colSums(V))
      Den <- t(su)%*%su
      Db <- (t(V)%*%R%*%V)/Den
      Db <- Db-diag(diag(Db))+diag(m)-Matrix::tril(Db)+Matrix::tril(t(Db))
      out <- ultrcorr(R, Db, V)
      Rb <- out[[1]]
      A <- out[[2]]
      B <- out[[3]]
      levfus <- out[[4]]
    } else if (m == 1) {
      Rb <- diag(m)
    }
    if (m > 1 && m < p) {
      dRw <- diag(Rw)
      minRw <- min(dRw)
      maxRb <- max(Rb-diag(m))
      dm <- maxRb-minRw
      if (dm > 0) {
        idRw <- which(dRw<maxRb)
        dRw[idRw] <- maxRb
        Rw <- diag(dRw)
      }
    }
    Rt <- V %*% (Rb-diag(m)) %*% t(V) + V %*% Rw %*% t(V) - diag(diag(V%*%Rw%*%t(V))) + diag(p)
    fo <- tr(t(R-Rt)%*%(R-Rt))/ts
    fmin <- fo
    fdif <- fmin
    while (fdif > eps && it <=  maxiter) {
      it <- it + 1
      fo <- fmin
      for (i in 1:p){
        posmin <- pp[V[i, ] == 1]
        posmin <- posmin[1]
        V[i, ] <- Im[posmin, ]
        for (z in 1:m){
          V[i, ] <- Im[z, ]
          if (sum(V[, posmin]) > 0){
            Rw <- t(V) %*% (R-diag(p)) %*% V %*% MASS::ginv((t(V) %*% V)^2 - (t(V) %*% V),  tol = .Machine$double.xmin)
            if (m > 1) {
              Rw <- diag(diag(Rw))
            }
            if (any(colSums(V) == 1)){
              Rw[which(colSums(V) == 1), which(colSums(V) == 1)] <- diag(length(which(colSums(V) == 1)))
            }
            Rw_wo <- Rw
            su <- t(colSums(V))
            ss <- (su==0)
            su <- su+ss
            Den <- t(su)%*%su
            Db <- (t(V)%*%R%*%V)/Den
            Db <- Db-diag(diag(Db))+diag(m)-Matrix::tril(Db)+Matrix::tril(t(Db))

            out <- ultrcorr(R, Db, V)
            Rb <- out[[1]]
            A <- out[[2]]
            B <- out[[3]]
            levfus <- out[[4]]

            dRw <- diag(Rw)
            minRw <- min(dRw)
            maxRb <- max(Rb-diag(m))
            dm <- maxRb-minRw
            if (dm > 0){
              idRw <- which(dRw<maxRb)
              dRw[idRw] <- maxRb
              Rw <- diag(dRw)
            }

            Rt <- V %*% (Rb-diag(m)) %*% t(V) + V %*% Rw %*% t(V) - diag(diag(V%*%Rw%*%t(V))) + diag(p)
            ff <- tr(t(R-Rt)%*%(R-Rt))/ts
            if (ff <= fmin){
              fmin <- ff
              posmin <- z
              Rwit <- Rw
              Rbit <- Rb
              Rtit <- Rt
              Rw_woit <- Rw_wo
              Ait <- A
              Bit <- B
              levfusit <- levfus
            }
          }
        }
        V[i, ] <- Im[posmin, ]
      }
      fdif <- fo-fmin
    }

    fmm[loop, 1] <- fmin

    # Case: m == 1
    if (m == 1){
      VOtt <- V
      RwOtt <- Rw
      Rw_woOtt <- Rw_wo
      RbOtt <- Rb
      RtOtt <- Rt
      fOtt <- fo
      loopOtt <- 1
      iterOtt <- 1
      AOtt <- 0
      BOtt <- 0
      levfusOtt <- Rw
      break
    }

    # Case: m == p
    if  (m == p){
      VOtt <- V
      RwOtt <- Rw
      Rw_woOtt <- Rw_wo
      RbOtt <- Rb
      RtOtt <- Rt
      fOtt <- fo
      loopOtt <- 1
      iterOtt <- 1
      AOtt <- A
      BOtt <- B
      levfusOtt <- levfus
    }

    # Optimal solution
    if (fmin < fOtt){
      VOtt <- V
      RwOtt <- Rwit
      Rw_woOtt <- Rw_woit
      RbOtt <- Rbit
      RtOtt <- Rtit
      fOtt <- fmin
      loopOtt <- loop
      iterOtt <- it
      AOtt <- Ait
      BOtt <- Bit
      levfusOtt <- levfusit
    }
  }

  nfOtt <- length(which(round(fmm, 6) == round(fOtt, 6)))

  # Graphical Representation
  Dt <- matrix(1, p, p) - RtOtt

  if (m == 1) {
    m <- p
    Rbp <- as.numeric(Rw)*(matrix(1, m, m) - diag(m))
    Vp <- diag(m)
    ultrcov.sim <- ultrcov_sim(Rbp, Vp)
    Ap <- ultrcov.sim$A
    Bp <- ultrcov.sim$B
    levfusp <- ultrcov.sim$levfus
    Rwp <- diag(m)
    m <- dim(R)[1]
  } else if (m > 1) {
    Vp <- VOtt
    Ap <- AOtt
    Bp <- BOtt
    levfusp <- levfusOtt
    Rbp <- RbOtt
    Rwp <- RwOtt
  }

  list <- list()
  for (i in 1:m) {
    list[[i]] <- list(i)
  }
  for (i in 1:length(Ap)) {
    for (j in 1:length(list)) {
      if (Ap[i] %in% list[[j]]) {
        z <- j
      }
    }
    for (j in 1:length(list)) {
      if (Bp[i] %in% list[[j]]) {
        y <- j
      }
    }
    if (Bp[i] == Ap[1]) {
      w <- c(list[[y]], list[[z]])
    } else{
      w <- c(list[[z]], list[[y]])
    }
    list[[z]] <- w
    list[[y]] <- NULL
  }

  outV <- unlist(list)
  outV <- as.matrix(outV)
  outV <- t(outV)

  PathDiagram(p, m, Vp, outV, Ap, Bp, levfusp, Rbp - diag(m), Rwp, diag(m), 1, cluster_names = "")

  return(list(call = call,
                V = VOtt,
                Rt = RtOtt,
                Rw = RwOtt,
                Rb = RbOtt,
                of = fOtt,
                loop = loopOtt,
                iter = iterOtt))
}


ultrcorr <- function(R,
                     Rb,
                     V){

  m <- dim(Rb)[1]
  p <- dim(V)[1]
  A <- matrix(0, m-1, 1)
  B <- matrix(0, m-1, 1)
  levfus <- matrix(0, m-1, 1)
  PP <- matrix(0, m, 1)
  P <- matrix(0, m, 1)

  for (q in 1:m){
    PP[q] <- q
  }
  class <- t(colSums(V))

  for (ustep in 1:(m-1)){
    rmax <- -Inf
    for (i in 1:(m-1)){
      if (PP[i] == i){
        for (j in (i+1):m){
          if (PP[j] == j){
            if (Rb[i,j] >= rmax){
              ic <- i
              jc <- j
              rmax <- Rb[i, j]
            }
          }
        }
      }
    }

    for (j in 1:m){
      if (PP[j] == jc){
        PP[j] <- ic
      }
    }
    A[ustep] <- ic
    B[ustep] <- jc
    levfus[ustep] <- rmax

    for (i in 1:m){
      if (i != ic && PP[i] == i){
        rs <- (class[ic]*Rb[ic, i]+class[jc]*Rb[jc, i])/(class[ic]+class[jc])
        Rb[ic,i] <- rs
        Rb[i,ic] <- rs
      }
    }
    class[ic] <- class[ic]+class[jc]
  }

  for (i in 1:m){
    P[i] <- i
  }

  RU <- matrix(0, m, m)

  for (k in 1:(m-1)){
    for (i in A[k]:m){
      if (P[i] == A[k]){
        for (j in B[k]:m){
          if (P[j] == B[k]){
            RU[i,j] <- levfus[k]
            RU[j,i] <- levfus[k]
          }
        }
      }
    }
    for (i in A[k]:m){
      if (P[i] == B[k]){
        P[i] <- A[k]
      }
    }
  }
  RU <- RU + diag(m)

  return(list(RU = RU,
              A = A,
              B = B,
              levfus = levfus))
}












