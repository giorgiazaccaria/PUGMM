ultrcov_sim <-
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
    return(list(Su = SU,
                A = A,
                B = B,
                levfus = levfus))
  }






flipudr <- function (a) {
  a <- a[nrow(a):1, ]
  return(a)
}







