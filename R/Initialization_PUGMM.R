init_PUGMM <-
  function(X,
           G,
           m,
           initG,
           initm,
           model.name = c("EUUU", "EUUE", "EUEE", "EEEU", "EEEE", "EEEF", "EEFF", "EFFF", "FIII", "FIIF", "FIFF", "FFFI", "FFFF"),
           pd) {
    model.name <- match.arg(model.name)
    n.obs <- dim(X)[1]
    p <- dim(X)[2]
    Sigma <- vector(mode = "list", length = G)
    model <- list()
    if (initG == "kmeans") {
      if (G > n.obs - 2){
        w <- mclust::unmap(stats::kmeans(X, centers = G, iter.max = 100)$cluster)
      } else {
      w <- mclust::unmap(ClusterR::KMeans_rcpp(X, clusters = G)$clusters)}
    } else if (initG == "kmeansf") {
      w <- unname(ppclust::fcm(X, G, iter.max = 100, con.val = sqrt(.Machine$double.eps))$u)
    } else if (initG == "random"){
      w <- rand.member(n.obs, G)
    }
    w[which(w < .Machine$double.xmin, arr.ind = TRUE)] <-
      .Machine$double.xmin
    label <- apply(w, 1, which.max)
    V <- vector(mode = "list", length = G)
    if (m == 1) {
      V[seq(1:G)] <- list(as.matrix(rep(1, p)))
    } else if (m == p) {
      V[seq(1:G)] <- list(diag(p))
    } else {
      if (initm == "ucms") {
        switch(
          model.name,
          "EUUU" = ,
          "EUUE" = ,
          "EUEE" = ,
          "EEEU" = ,
          "EEEE" = ,
          "EEEF" = ,
          "EEFF" = ,
          "EFFF" = {
            S <- matrix(0, p, p)
            for (g in 1:G) {
              if (nrow(X[label == g, , drop = FALSE]) != 1) {
                S <- S + (stats::cov(X[label == g, ]) * colSums(w)[g])
              }
            }
            S <- S / n.obs
            V[seq(1:G)] <- list(UCMSigmaV(X, S, m, rndst = 1, model.name = model.name))
          },
          "FIII" = ,
          "FIIF" = ,
          "FIFF" = ,
          "FFFI" = ,
          "FFFF" = {
            for (g in 1:G) {
              if (nrow(X[label == g, , drop = FALSE]) == 1) {
                V[[g]] <- rand.member(p, m)
              } else {
                S <- matrix(0, p, p)
                S <- stats::cov(X[label == g, ])
                V[[g]] <- UCMSigmaV(X, S, m, rndst = 1, model.name = model.name)
              }
            }
          }
        )
      } else if (initm == "random") {
        switch(
          model.name,
          "EUUU" = ,
          "EUUE" = ,
          "EUEE" = ,
          "EEEU" = ,
          "EEEE" = ,
          "EEEF" = ,
          "EEFF" = ,
          "EFFF" = {
            V[seq(1:G)] <- list(rand.member(p, m))
          },
          "FIII" = ,
          "FIIF" = ,
          "FIFF" = ,
          "FFFI" = ,
          "FFFF" = {
            V[seq(1:G)] <- replicate(G, list(rand.member(p, m)))
          }
        )
      }
    }
    model$pp <- prior(w)
    comp.param <- component_param(X, w)
    model$mu <- comp.param$mu
    Sigma <- comp.param$Sigma
    Sbar <- apply(sweep(array(unlist(Sigma), dim = c(p, p, G)), 3, model$pp, "*"), MARGIN = c(1, 2), sum)
    if (model.name == "EFFF" | model.name == "FIII" | model.name == "FIIF" | model.name == "FIFF" | model.name == "FFFI" | model.name == "FFFF") {
      for (g in 1:G) {
        Sv <- estimate_Sv(Sigma[[g]], V[[g]], model.name = model.name)
        Sw <- estimate_Sw(Sigma[[g]], Sv, V[[g]], model.name = model.name)
        Sb <- estimate_Sb(Sigma[[g]], V[[g]], model.name = model.name)
        if (m > 1) {
          Sw <- check_constraint_SwSb(Sw, Sb, model.name = model.name)$Sw
          if (m == p) {
            Sv <- Sw
          }
        }
        if (m < p) {
          Sv <- check_constraint_SvSw(Sv, Sw, model.name = model.name)$Sv
        }
        St.psd <- psd_sigma(V[[g]] %*% (Sw + Sb) %*% t(V[[g]]) - diag(diag(V[[g]] %*% Sw %*% t(V[[g]]))) + diag(diag(V[[g]] %*% Sv %*% t(V[[g]]))))
        St <- St.psd$S
        if (St.psd$constr.psd != 0) {
          Sv <- estimate_Sv(St, V[[g]], model.name = model.name)
          Sw <- estimate_Sw(St, Sv, V[[g]], model.name = model.name)
          Sb <- estimate_Sb(St, V[[g]], model.name = model.name)
          if (m > 1) {
            Sw <- check_constraint_SwSb(Sw, Sb, model.name = model.name)$Sw
            if (m == p) {
              Sv <- Sw
            }
          }
          if (m < p) {
            Sv <- check_constraint_SvSw(Sv, Sw, model.name = model.name)$Sv
          }
          St <- V[[g]] %*% (Sw + Sb) %*% t(V[[g]]) - diag(diag(V[[g]] %*% Sw %*% t(V[[g]]))) + diag(diag(V[[g]] %*% Sv %*% t(V[[g]])))
        }
        while (min(eigen(St, symmetric = TRUE, only.values = TRUE)$values) <= 0 || det(St) <= sqrt(.Machine$double.eps)) {
          a <- abs(min(eigen(St, symmetric = TRUE, only.values = TRUE)$values)) + sqrt(.Machine$double.eps)
          St <- St + a * diag(p)
          Sv <- Sv + a * diag(m)
        }
        if (any(apply(V[[g]], 2, sum) == 1)) {
          Sw[which(apply(V[[g]], 2, sum) == 1), which(apply(V[[g]], 2, sum) == 1)] <- Sv[which(apply(V[[g]], 2, sum) == 1), which(apply(V[[g]], 2, sum) == 1)]
        }
        model$Sigma[[g]] <- St
        model$V[[g]] <- V[[g]]
        model$Sv[[g]] <- Sv
        model$Sw[[g]] <- Sw
        model$Sb[[g]] <- Sb
      }
    } else {
      if (model.name == "EUUU" | model.name == "EUUE" | model.name == "EUEE" | model.name == "EEEU" | model.name == "EEEE") {
        Sv <- estimate_Sv(Sbar, V[[1]], model.name = model.name)
        Sw <- estimate_Sw(Sbar, Sv, V[[1]], model.name = model.name)
        Sb <- estimate_Sb(Sbar, V[[1]], model.name = model.name)
        if (m > 1) {
          Sw <- check_constraint_SwSb(Sw, Sb, model.name = model.name)$Sw
          if (m == p) {
            Sv <- Sw
          }
        }
        if (m < p) {
          Sv <- check_constraint_SvSw(Sv, Sw, model.name = model.name)$Sv
        }
        St.psd <- psd_sigma(V[[1]] %*% (Sw + Sb) %*% t(V[[1]]) - diag(diag(V[[1]] %*% Sw %*% t(V[[1]]))) + diag(diag(V[[1]] %*% Sv %*% t(V[[1]]))))
        St <- St.psd$S
        if (St.psd$constr.psd != 0) {
          Sv <- estimate_Sv(St, V[[1]], model.name = model.name)
          Sw <- estimate_Sw(St, Sv, V[[1]], model.name = model.name)
          Sb <- estimate_Sb(St, V[[1]], model.name = model.name)
          if (m > 1) {
            Sw <- check_constraint_SwSb(Sw, Sb, model.name = model.name)$Sw
            if (m == p) {
              Sv <- Sw
            }
          }
          if (m < p) {
            Sv <- check_constraint_SvSw(Sv, Sw, model.name = model.name)$Sv
          }
          St <- V[[1]] %*% (Sw + Sb) %*% t(V[[1]]) - diag(diag(V[[1]] %*% Sw %*% t(V[[1]]))) + diag(diag(V[[1]] %*% Sv %*% t(V[[1]])))
        }
        while (min(eigen(St, symmetric = TRUE, only.values = TRUE)$values) <= 0 || det(St) <= sqrt(.Machine$double.eps)) {
          a <- abs(min(eigen(St, symmetric = TRUE, only.values = TRUE)$values)) + sqrt(.Machine$double.eps)
          St <- St + a * diag(p)
          Sv <- Sv + a * diag(m)
        }
        if (any(apply(V[[1]], 2, sum) == 1)) {
          Sw[which(apply(V[[1]], 2, sum) == 1), which(apply(V[[1]], 2, sum) == 1)] <- Sv[which(apply(V[[1]], 2, sum) == 1), which(apply(V[[1]], 2, sum) == 1)]
        }
        model$Sigma[seq(1:G)] <- list(St)
        model$V[seq(1:G)] <- list(V[[1]])
        model$Sv[seq(1:G)] <- list(Sv)
        model$Sw[seq(1:G)] <- list(Sw)
        model$Sb[seq(1:G)] <- list(Sb)
      } else if (model.name == "EEEF") {
        Sv <- estimate_Sv(Sbar, V[[1]], model.name = model.name)
        Sw <- estimate_Sw(Sbar, Sv, V[[1]], model.name = model.name)
        check.maxSb <- (-.Machine$double.xmax)
        for (g in 1:G) {
          Sb <- estimate_Sb(Sigma[[g]], V[[1]], model.name = model.name)
          maxSb <- max(max(Sb + (min(min(Sb)) * diag(m))))
          if (check.maxSb < maxSb) {
            check.maxSb <- maxSb
            which.maxSb <- g
          }
          model$Sb[[g]] <- Sb
        }
        if (m > 1) {
          Sw <- check_constraint_SwSb(Sw, model$Sb[[which.maxSb]], model.name = model.name)$Sw
          if (m == p) {
            Sv <- Sw
          }
        }
        if (m < p) {
          Sv <- check_constraint_SvSw(Sv, Sw, model.name = model.name)$Sv
        }
        constr.psd <- 0
        for (g in 1:G) {
          St.psd <- psd_sigma(V[[1]] %*% (Sw + model$Sb[[g]]) %*% t(V[[1]]) - diag(diag(V[[1]] %*% Sw %*% t(V[[1]]))) + diag(diag(V[[1]] %*% Sv %*% t(V[[1]]))))
          model$Sigma[[g]] <- St.psd$S
          if (St.psd$constr.psd != 0) {
            model$Sb[[g]] <- estimate_Sb(model$Sigma[[g]], V[[1]], model.name = model.name)
            constr.psd <- constr.psd + 1
          }
        }
        if (constr.psd > 0) {
          Sbar <- apply(sweep(array(unlist(model$Sigma), dim = c(p, p, G)), 3, model$pp, "*"), MARGIN = c(1, 2), sum)
          Sv <- estimate_Sv(Sbar, V[[1]], model.name = model.name)
          Sw <- estimate_Sw(Sbar, Sv, V[[1]], model.name = model.name)
          which.maxSb <- which.max(sapply(model$Sb, function(x){max(max(x + (min(min(x)) * diag(m))))}))
          if (m > 1) {
            Sw <- check_constraint_SwSb(Sw, model$Sb[[which.maxSb]], model.name = model.name)$Sw
            if (m == p) {
              Sv <- Sw
            }
          }
          if (m < p) {
            Sv <- check_constraint_SvSw(Sv, Sw, model.name = model.name)$Sv
          }
          for (g in 1:G) {
            model$Sigma[[g]] <- V[[1]] %*% (Sw + model$Sb[[g]]) %*% t(V[[1]]) - diag(diag(V[[1]] %*% Sw %*% t(V[[1]]))) + diag(diag(V[[1]] %*% Sv %*% t(V[[1]])))
          }
        }
        max.a <- 0
        for (g in 1:G) {
          if (min(eigen(model$Sigma[[g]], symmetric = TRUE, only.values = TRUE)$values) <= 0 || det(model$Sigma[[g]]) <= sqrt(.Machine$double.eps)) {
            eig.St <- abs(min(eigen(model$Sigma[[g]], symmetric = TRUE, only.values = TRUE)$values)) + sqrt(.Machine$double.eps)
            if (max.a < eig.St) {
              max.a <- eig.St
            }
          }
        }
        Sv <- Sv + max.a * diag(m)
        if (any(apply(V[[1]], 2, sum) == 1)) {
          Sw[which(apply(V[[1]], 2, sum) == 1), which(apply(V[[1]], 2, sum) == 1)] <- Sv[which(apply(V[[1]], 2, sum) == 1), which(apply(V[[1]], 2, sum) == 1)]
        }
        model$Sigma[seq(1:G)] <- lapply(model$Sigma, function (x) return(x + max.a * diag(p)))
        model$V[seq(1:G)] <- list(V[[1]])
        model$Sv[seq(1:G)] <- list(Sv)
        model$Sw[seq(1:G)] <- list(Sw)
      } else if (model.name == "EEFF") {
        Sv <- estimate_Sv(Sbar, V[[1]], model.name = model.name)
        for (g in 1:G) {
          Sw <- estimate_Sw(Sigma[[g]], Sv, V[[1]], model.name = model.name)
          Sb <- estimate_Sb(Sigma[[g]], V[[1]], model.name = model.name)
          if (m > 1) {
            Sw <- check_constraint_SwSb(Sw, Sb, model.name = model.name)$Sw
          }
          model$Sw[[g]] <- Sw
          model$Sb[[g]] <- Sb
        }
        Sw.check <- diag(apply(matrix(unlist(lapply(model$Sw, diag)), m, G), 1, function(x) max(abs(x))), m, m)
        if (m < p) {
          Sv <- check_constraint_SvSw(Sv, Sw.check, model.name = model.name)$Sv
        }
        constr.psd <- 0
        for (g in 1:G) {
          St.psd <- psd_sigma(V[[1]] %*% (model$Sw[[g]] + model$Sb[[g]]) %*% t(V[[1]]) - diag(diag(V[[1]] %*% model$Sw[[g]] %*% t(V[[1]]))) + diag(diag(V[[1]] %*% Sv %*% t(V[[1]]))))
          model$Sigma[[g]] <- St.psd$S
          if (St.psd$constr.psd != 0) {
            model$Sb[[g]] <- estimate_Sb(model$Sigma[[g]], V[[1]], model.name)
            constr.psd <- constr.psd + 1
          }
        }
        if (constr.psd > 0) {
          Sbar <- apply(sweep(array(unlist(model$Sigma[[g]]), dim = c(p, p, G)), 3, model$pp, "*"), MARGIN = c(1, 2), sum)
          count.SwSb <- 0
          count.SvSw <- 0
          Sv <- estimate_Sv(Sbar, V[[1]], model.name)
          for (g in 1:G) {
            Sw <- estimate_Sw(model$Sigma[[g]], Sv, V[[1]], model.name)
            if (m > 1) {
              Sw.constr <- check_constraint_SwSb(Sw, model$Sb[[g]], model.name)
              Sw <- Sw.constr$Sw
              count.SwSb <- count.SwSb + Sw.constr$constr.SwSb
            }
            model$Sw[[g]] <- Sw
          }
          Sw.check <- diag(apply(matrix(unlist(lapply(model$Sw, diag)), m, G), 1, function(x) max(abs(x))), m, m)
          if (m < p) {
            Sv.constr <- check_constraint_SvSw(Sv, Sw.check, model.name)
            Sv <- Sv.constr$Sv
            count.SvSw <- Sv.constr$constr.SvSw
          }
          for (g in 1:G) {
            model$Sigma[[g]] <- V[[1]] %*% (model$Sw[[g]] + model$Sb[[g]]) %*% t(V[[1]]) - diag(diag(V[[1]] %*% model$Sw[[g]] %*% t(V[[1]]))) + diag(diag(V[[1]] %*% Sv %*% t(V[[1]])))
          }
        }
        max.a <- 0
        for (g in 1:G) {
          if (min(eigen(model$Sigma[[g]], symmetric = TRUE, only.values = TRUE)$values) <= 0 || det(model$Sigma[[g]]) <= sqrt(.Machine$double.eps)) {
            eig.St <- abs(min(eigen(model$Sigma[[g]], symmetric = TRUE, only.values = TRUE)$values)) + sqrt(.Machine$double.eps)
            if (max.a < eig.St) {
              max.a <- eig.St
            }
          }
        }
        Sv <- Sv + max.a * diag(m)
        if (any(apply(V[[1]], 2, sum) == 1)) {
          model$Sw <- lapply(model$Sw, function(x){x[which(apply(V[[1]], 2, sum) == 1), which(apply(V[[1]], 2, sum) == 1)] <- Sv[which(apply(V[[1]], 2, sum) == 1), which(apply(V[[1]], 2, sum) == 1)]
          return(x)})
        }
        model$Sigma[seq(1:G)] <- lapply(model$Sigma, function (x) return(x + max.a * diag(p)))
        model$V[seq(1:G)] <- list(V[[1]])
        model$Sv[seq(1:G)] <- list(Sv)
      }
    }
    model$post <- post(X, model$pp, model$mu, model$Sigma)
    return(list(
      label = label,
      model = model,
      Sg = Sigma
    ))
  }
