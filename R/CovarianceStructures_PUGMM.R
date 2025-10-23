estimate_Sv <-
  function(Sigma,
           V,
           model.name = c("EUUU", "EUUE", "EUEE", "EEEU", "EEEE", "EEEF", "EEFF", "EFFF", "FIII", "FIIF", "FIFF", "FFFI", "FFFF")) {
    model.name <- match.arg(model.name)

    p <- dim(V)[1]
    m <- dim(V)[2]
    Sv <- MASS::ginv(V,  tol = .Machine$double.xmin) %*% diag(diag(Sigma)) %*% V
    if (m > 1) {
      switch(
        model.name,
        "EUUU" = ,
        "EUUE" = ,
        "EUEE" = ,
        "FIII" = ,
        "FIIF" = ,
        "FIFF" = {
          Sv <- (sum(diag(Sv)) / m) * diag(m)
        },
        "EEEU" = ,
        "EEEE" = ,
        "EEEF" = ,
        "EEFF" = ,
        "EFFF" = ,
        "FFFI" = ,
        "FFFF" = {
          if (m < p) {
            Sv <- diag(diag(Sv))
          }
        }
      )
    }
    return(Sv)
  }


estimate_Sw <-
  function(Sigma,
           Sv,
           V,
           model.name = c("EUUU", "EUUE", "EUEE", "EEEU", "EEEE", "EEEF", "EEFF", "EFFF", "FIII", "FIIF", "FIFF", "FFFI", "FFFF")) {
    model.name <- match.arg(model.name)

    p <- dim(V)[1]
    m <- dim(V)[2]
    if (m == p) {
      Sw <- Sv
    } else {
      Sw <- t(V) %*% (Sigma - diag(diag(V %*% Sv %*% t(V)))) %*% V %*% MASS::ginv((t(V) %*% V)^2 - (t(V) %*% V),  tol = .Machine$double.xmin)
      if (m > 1) {
        switch(
          model.name,
          "EUUU" = ,
          "EUUE" = ,
          "FIII" = ,
          "FIIF" = {
            Sw <-(sum(diag(Sw)) / m) * diag(m)
          },
          "EUEE" = ,
          "EEEU" = ,
          "EEEE" = ,
          "EEEF" = ,
          "EEFF" = ,
          "EFFF" = ,
          "FIFF" = ,
          "FFFI" = ,
          "FFFF" = {
            Sw <- diag(diag(Sw))
          }
        )
      }
    }
    return(Sw)
  }



estimate_Sb <-
  function(Sigma,
           V,
           model.name = c("EUUU", "EUUE", "EUEE", "EEEU", "EEEE", "EEEF", "EEFF", "EFFF", "FIII", "FIIF", "FIFF", "FFFI", "FFFF")) {
    model.name <- match.arg(model.name)

    p <- dim(V)[1]
    m <- dim(V)[2]
    if (m == 1) {
      Sb <- as.matrix(0)
    } else {
    switch(
      model.name,
      "EUUU" = ,
      "EEEU" = ,
      "FIII" = ,
      "FFFI" = {
        Sb <- ((MASS::ginv(V,  tol = .Machine$double.xmin) %*% Sigma %*% MASS::ginv(t(V),  tol = .Machine$double.xmin)) * (matrix(1, m, m) - diag(m))) %*% MASS::ginv(matrix(1, m, m) - diag(m),  tol = .Machine$double.xmin)
        Sb <- (sum(diag(Sb)) / m) * (matrix(1, m, m) - diag(m))
      },
      "EUUE" = ,
      "EUEE" = ,
      "EEEE" = ,
      "EEEF" = ,
      "EEFF" = ,
      "EFFF" = ,
      "FIIF" = ,
      "FIFF" = ,
      "FFFF" = {
        Db <- MASS::ginv(V) %*% Sigma %*% MASS::ginv(t(V),  tol = .Machine$double.xmin)
        Sb <- ultrcov(Db, V)
        }
    )
    }
   return(Sb)
  }


check_constraint_SwSb <-
  function(Sw,
           Sb,
           model.name = c("EUUU", "EUUE", "EUEE", "EEEU", "EEEE", "EEEF", "EEFF", "EFFF", "FIII", "FIIF", "FIFF", "FFFI", "FFFF")) {
    model.name <- match.arg(model.name)
    m <- dim(Sw)[1]
    constr.SwSb <- 0
    switch(
      model.name,
      "EUUU" = ,
      "FIII" = {
        if (Sb[1, 2] - Sw[1, 1] > 0) {
          constr.SwSb <- 1
          Sw <- Sb[1, 2] * diag(m)
        }
      },
      "EUUE" = ,
      "FIIF" = {
        maxSb <- max(max(Sb + (min(min(Sb)) * diag(m))))
        if (maxSb - Sw[1, 1] > 0) {
          constr.SwSb <- 1
          Sw <- maxSb * diag(m)
        }
      },
      "EUEE" = ,
      "EEEU" = ,
      "EEEE" = ,
      "EEEF" = ,
      "EEFF" = ,
      "EFFF" = ,
      "FIFF" = ,
      "FFFI" = ,
      "FFFF" = {
        dSw <- diag(Sw)
        minSw <- min(dSw)
        maxSb <- max(max(Sb + (min(min(Sb)) * diag(m))))
        if (maxSb - minSw > 0) {
          constr.SwSb <- length(dSw[dSw < maxSb])
          dSw[dSw < maxSb] <- maxSb;
          Sw <- diag(dSw)
        }
      }
    )
    return(list(Sw = Sw,
                constr.SwSb = constr.SwSb))
  }


check_constraint_SvSw <-
  function(Sv,
           Sw,
           model.name = c("EUUU", "EUUE", "EUEE", "EEEU", "EEEE", "EEEF", "EEFF", "EFFF", "FIII", "FIIF", "FIFF", "FFFI", "FFFF")) {
    model.name <- match.arg(model.name)
    m <- dim(Sv)[1]
    G <- length(Sw)
    constr.SvSw <- 0
    switch(
      model.name,
      "EUUU" = ,
      "EUUE" = ,
      "FIII" = ,
      "FIIF" = {
        if (Sv[1, 1] < abs(Sw[1, 1])) {
          constr.SvSw <- 1
          Sv <- as.matrix(abs(Sw[1, 1]) * diag(m) + sqrt(.Machine$double.eps))
        }
      },
      "EUEE" = ,
      "FIFF" = {
        maxSw <- max(abs(diag(Sw)))
        if (Sv[1, 1] < maxSw) {
          constr.SvSw <- 1
          Sv <- as.matrix(maxSw * diag(m) + sqrt(.Machine$double.eps))
        }
      },
      "EEEU" = ,
      "EEEE" = ,
      "EEEF" = ,
      "EEFF" = ,
      "EFFF" = ,
      "FFFI" = ,
      "FFFF" = {
        if (m == 1) {
          if (Sv < abs(Sw)) {
            constr.SvSw <- 1
            Sv <- Sw
          }
        } else if (m > 1) {
          dSw <- abs(diag(Sw))
          dSv <- diag(Sv)
          constr.SvSw <- length(dSv[dSv < dSw])
          dSv[dSv < dSw] <- dSw[dSv < dSw] + sqrt(.Machine$double.eps)
          Sv <- diag(dSv)
        }
      }
    )
    return(list(Sv = Sv,
                constr.SvSw = constr.SvSw))
  }


psd_sigma <-
  function(S) {
    p <- dim(S)[1]
    constr.psd <- 0
    if (min(eigen(S, symmetric = TRUE)$values) <= 0) {
      constr.psd <- 1
      S.sd <- eigen(S, symmetric = TRUE)
      H <- S.sd$vectors %*% abs(diag(S.sd$values)) %*% t(S.sd$vectors)
      S <- (S + H)/2
    }
    return(list(S = S,
                constr.psd = constr.psd))
  }

