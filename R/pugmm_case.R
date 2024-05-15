pugmm_case <-
  function(X,
           G,
           m,
           model.name,
           model,
           maxiter,
           tol,
           stop,
           rndstart,
           initG,
           initm,
           gaussian) {
    call <- match.call()
    n.obs <- dim(X)[1]
    p <- dim(X)[2]
   change_model <- 0
    if(m == p && (model.name == "EUEE" | model.name == "EEFF" | model.name == "FIFF")) {
      if ((model.name == "EUEE" && "EUUE" %in% model) |
          (model.name == "EEFF" && "EEEF" %in% model) |
          (model.name == "FIFF" && "FIIF" %in% model) ){
      ans <- c(list(call = call,
                    X = X,
                    G = G,
                    m = m,
                    pm.free = +Inf,
                    bic = NA,
                    model.name = model.name))
      orderedNames <- c("call","X","G","m","bic","model.name")
      structure(ans[orderedNames], class = "PUGMM_big")
      return()
      }
      if (model.name == "EUEE" && !("EUUE" %in% model)){
        model.name <- "EUUE"
        change_model <- 1
      }
      if (model.name == "EEFF" && !("EEEF" %in% model)){
        model.name <- "EEEF"
        change_model <- 1
      }
      if (model.name == "FIFF" && !("FIIF" %in% model)){
        model.name <- "FIIF"
        change_model <- 1
      }
    }
   IQ <- diag(m)
    wpc <- G * (G + 1) / 2
    for (loop in 1:rndstart) {
      count <- 1
      conv <- 1
      wp <- 0
      init <- init_PUGMM(X, G, m, model.name = model.name, initG = initG, initm = initm)
      model <- init$model
      label <- init$label
      Sigma <- init$Sg
      loglik <-
        loglikH_EI_G_pugmm(
          X,
          model$pp,
          model$mu,
          model$Sigma,
          model$Sv,
          model$Sw,
          model$Sb,
          model$V,
          model$post,
          gaussian
        )
      while (conv > tol && count <= maxiter) {
        if (count > 1) {
          model$post <- post(X, model$pp, model$mu, model$Sigma)
          label <- apply(model$post, 1, which.max)
        }
        fmax <- loglik
        if (sum(unique(label)) != wpc) {
          wp <- 1
          loglik <- -.Machine$double.xmax
          break
        }
        model$pp <- prior(model$post)
        comp.param <- component_param(X, model$post)
        model$mu <- comp.param$mu
        Sigma <- comp.param$Sigma
        Sbar <- apply(sweep(array(unlist(Sigma), dim = c(p, p, G)), 3, model$pp, "*"), MARGIN = c(1, 2), sum)
        if (model.name == "FIII" | model.name == "FIIF" | model.name == "FIFF" | model.name == "FFFI" | model.name == "FFFF") {
          constr.SwSb <- 0
          constr.SvSw <- 0
          for (g in 1:G) {
            Stit <- model$Sigma[[g]]
            Svit <- model$Sv[[g]]
            Swit <- model$Sw[[g]]
            Sbit <- model$Sb[[g]]
            constr.SwSbit <- 0
            constr.SvSwit <- 0
            V0 <- model$V[[g]]
            for (i in 1:p) {
              posmax <- which(V0[i, ] == 1)
              V0[i, ] <- IQ[posmax, ]
              for (h in 1:m) {
                V0[i, ] <- IQ[h, ]
                if (sum(V0[, posmax]) > 0) {
                  count.SwSb <- 0
                  count.SvSw <- 0
                  Sv <- estimate_Sv(Sigma[[g]], V0, model.name)
                  Sw <- estimate_Sw(Sigma[[g]], Sv, V0, model.name)
                  Sb <- estimate_Sb(Sigma[[g]], V0, model.name)
                  if (m > 1) {
                    Sw.constr <- check_constraint_SwSb(Sw, Sb, model.name)
                    Sw <- Sw.constr$Sw
                    count.SwSb <- Sw.constr$constr.SwSb
                    if (m == p) {
                      Sv <- Sw
                    }
                  }
                  if (m < p) {
                    Sv.constr <- check_constraint_SvSw(Sv, Sw, model.name)
                    Sv <- Sv.constr$Sv
                    count.SvSw <- Sv.constr$constr.SvSw
                  }
                  St.psd <- psd_sigma(V0 %*% (Sw + Sb) %*% t(V0) - diag(diag(V0 %*% Sw %*% t(V0))) + diag(diag(V0 %*% Sv %*% t(V0))))
                  St <- St.psd$S
                  if (St.psd$constr.psd != 0) {
                    count.SwSb <- 0
                    count.SvSw <- 0
                    Sv <- estimate_Sv(St, V0, model.name)
                    Sw <- estimate_Sw(St, Sv, V0, model.name)
                    Sb <- estimate_Sb(St, V0, model.name)
                    if (m > 1) {
                      Sw.constr <- check_constraint_SwSb(Sw, Sb, model.name)
                      Sw <- Sw.constr$Sw
                      count.SwSb <- Sw.constr$constr.SwSb
                      if (m == p) {
                        Sv <- Sw
                      }
                    }
                    if (m < p) {
                      Sv.constr <- check_constraint_SvSw(Sv, Sw, model.name)
                      Sv <- Sv.constr$Sv
                      count.SvSw <- Sv.constr$constr.SvSw
                    }
                    St <- V0 %*% (Sw + Sb) %*% t(V0) - diag(diag(V0 %*% Sw %*% t(V0))) + diag(diag(V0 %*% Sv %*% t(V0)))
                  }
                  while (min(eigen(St, symmetric = TRUE, only.values = TRUE)$values) <= 0 || det(St) <= .Machine$double.xmin) {
                    a <- abs(min(eigen(St, symmetric = TRUE, only.values = TRUE)$values)) + sqrt(.Machine$double.eps)
                    St <- St + a * diag(p)
                    Sv <- Sv + a * diag(m)
                  }
                  ff <-
                    loglikH_F_pugmm(
                      X,
                      g,
                      model$pp,
                      model$mu,
                      St,
                      model$Sigma,
                      Sv,
                      model$Sv,
                      Sw,
                      model$Sw,
                      Sb,
                      model$Sb,
                      V0,
                      model$V,
                      model$post,
                      gaussian
                    )
                  if (ff >= fmax) {
                    fmax <- ff
                    posmax <- h
                    Stit <- St
                    Svit <- Sv
                    Swit <- Sw
                    Sbit <- Sb
                    constr.SwSbit <- count.SwSb
                    constr.SvSwit <- count.SvSw
                  }
                }
              }
              V0[i, ] <- IQ[posmax, ]
            }
            if (any(apply(V0, 2, sum) == 1)) {
              Swit[which(apply(V0, 2, sum) == 1), which(apply(V0, 2, sum) == 1)] <- Svit[which(apply(V0, 2, sum) == 1), which(apply(V0, 2, sum) == 1)]
            }
            model$Sigma[[g]] <- Stit
            model$V[[g]] <- V0
            model$Sv[[g]] <- Svit
            model$Sw[[g]] <- Swit
            model$Sb[[g]] <- Sbit
            constr.SwSb <- constr.SwSb + constr.SwSbit
            constr.SvSw <- constr.SvSw + constr.SvSwit
          }
        } else if (model.name == "EUUU" | model.name == "EUUE" | model.name == "EUEE" | model.name == "EEEU" | model.name == "EEEE") {
          Stit <- model$Sigma[[1]]
          Svit <- model$Sv[[1]]
          Swit <- model$Sw[[1]]
          Sbit <- model$Sb[[1]]
          constr.SwSb <- 0
          constr.SvSw <- 0
          V0 <- model$V[[1]]
          for (i in 1:p) {
            posmax <- which(V0[i, ] == 1)
            V0[i, ] <- IQ[posmax, ]
            for (h in 1:m) {
              V0[i, ] <- IQ[h, ]
              if (sum(V0[, posmax]) > 0) {
                count.SwSb <- 0
                count.SvSw <- 0
                Sv <- estimate_Sv(Sbar, V0, model.name)
                Sw <- estimate_Sw(Sbar, Sv, V0, model.name)
                Sb <- estimate_Sb(Sbar, V0, model.name)
                if (m > 1) {
                  Sw.constr <- check_constraint_SwSb(Sw, Sb, model.name)
                  Sw <- Sw.constr$Sw
                  count.SwSb <- Sw.constr$constr.SwSb
                  if (m == p) {
                    Sv <- Sw
                  }
                }
                if (m < p) {
                  Sv.constr <- check_constraint_SvSw(Sv, Sw, model.name)
                  Sv <- Sv.constr$Sv
                  count.SvSw <- Sv.constr$constr.SvSw
                }
                St.psd <- psd_sigma(V0 %*% (Sw + Sb) %*% t(V0) - diag(diag(V0 %*% Sw %*% t(V0))) + diag(diag(V0 %*% Sv %*% t(V0))))
                St <- St.psd$S
                if (St.psd$constr.psd != 0) {
                  count.SwSb <- 0
                  count.SvSw <- 0
                  Sv <- estimate_Sv(St, V0, model.name)
                  Sw <- estimate_Sw(St, Sv, V0, model.name)
                  Sb <- estimate_Sb(St, V0, model.name)
                  if (m > 1) {
                    Sw.constr <- check_constraint_SwSb(Sw, Sb, model.name)
                    Sw <- Sw.constr$Sw
                    count.SwSb <- Sw.constr$constr.SwSb
                    if (m == p) {
                      Sv <- Sw
                    }
                  }
                  if (m < p) {
                    Sv.constr <- check_constraint_SvSw(Sv, Sw, model.name)
                    Sv <- Sv.constr$Sv
                    count.SvSw <- Sv.constr$constr.SvSw
                  }
                  St.psd <- psd_sigma(V0 %*% (Sw + Sb) %*% t(V0) - diag(diag(V0 %*% Sw %*% t(V0))) + diag(diag(V0 %*% Sv %*% t(V0))))
                  St <- St.psd$S
                }
                while (min(eigen(St, symmetric = TRUE, only.values = TRUE)$values) <= 0 || det(St) <= .Machine$double.xmin) {
                  a <- abs(min(eigen(St, symmetric = TRUE, only.values = TRUE)$values)) + sqrt(.Machine$double.eps)
                  St <- St + a * diag(p)
                  Sv <- Sv + a * diag(m)
                }
                ff <-
                  loglikH_EI_G_pugmm(
                    X,
                    model$pp,
                    model$mu,
                    replicate(G, list(St)),
                    replicate(G, list(Sv)),
                    replicate(G, list(Sw)),
                    replicate(G, list(Sb)),
                    replicate(G, list(V0)),
                    model$post,
                    gaussian
                  )
                if (ff >= fmax) {
                  fmax <- ff
                  posmax <- h
                  Stit <- St
                  Svit <- Sv
                  Swit <- Sw
                  Sbit <- Sb
                  constr.SwSb <- count.SwSb
                  constr.SvSw <- count.SvSw
                }
              }
            }
            V0[i, ] <- IQ[posmax, ]
          }
          if (any(apply(V0, 2, sum) == 1)) {
            Swit[which(apply(V0, 2, sum) == 1), which(apply(V0, 2, sum) == 1)] <- Svit[which(apply(V0, 2, sum) == 1), which(apply(V0, 2, sum) == 1)]
          }
          model$Sigma[seq(1:G)] <- list(Stit)
          model$V[seq(1:G)] <- list(V0)
          model$Sv[seq(1:G)] <- list(Svit)
          model$Sw[seq(1:G)] <- list(Swit)
          model$Sb[seq(1:G)] <- list(Sbit)
        } else if (model.name == "EEEF") {
          Stit <- model$Sigma
          Svit <- model$Sv[[1]]
          Swit <- model$Sw[[1]]
          Sbit <- model$Sb
          constr.SwSb <- 0
          constr.SvSw <- 0
          V0 <- model$V[[1]]
          Sb.list <- vector(mode = "list", length = G)
          St.list <- vector(mode = "list", length = G)
          for (i in 1:p) {
            posmax <- which(V0[i, ] == 1)
            V0[i, ] <- IQ[posmax, ]
            for (h in 1:m) {
              V0[i, ] <- IQ[h, ]
              if (sum(V0[, posmax]) > 0) {
                count.SwSb <- 0
                count.SvSw <- 0
                Sv <- estimate_Sv(Sbar, V0, model.name)
                Sw <- estimate_Sw(Sbar, Sv, V0, model.name)
                check.maxSb <- (-.Machine$double.xmax)
                for (g in 1:G) {
                  Sb <- estimate_Sb(Sigma[[g]], V0, model.name)
                  maxSb <- max(max(Sb + (min(min(Sb)) * diag(m))))
                  if (check.maxSb < maxSb) {
                    check.maxSb <- maxSb
                    which.maxSb <- g
                  }
                  Sb.list[[g]] <- Sb
                }
                if (m > 1) {
                  Sw.constr <- check_constraint_SwSb(Sw, Sb.list[[which.maxSb]], model.name)
                  Sw <- Sw.constr$Sw
                  count.SwSb <- Sw.constr$constr.SwSb
                  if (m == p) {
                    Sv <- Sw
                  }
                }
                if (m < p) {
                  Sv.constr <- check_constraint_SvSw(Sv, Sw, model.name)
                  Sv <- Sv.constr$Sv
                  count.SvSw <- Sv.constr$constr.SvSw
                }
                constr.psd <- 0
                for (g in 1:G) {
                  St.psd <- psd_sigma(V0 %*% (Sw + Sb.list[[g]]) %*% t(V0) - diag(diag(V0 %*% Sw %*% t(V0))) + diag(diag(V0 %*% Sv %*% t(V0))))
                  St.list[[g]] <- St.psd$S
                  if (St.psd$constr.psd != 0) {
                    Sb.list[[g]] <- estimate_Sb(St.list[[g]], V0, model.name)
                    constr.psd <- constr.psd + 1
                  }
                }
                if (constr.psd > 0) {
                  Sbar <- apply(sweep(array(unlist(St.list), dim = c(p, p, G)), 3, model$pp, "*"), MARGIN = c(1, 2), sum)
                  count.SwSb <- 0
                  count.SvSw <- 0
                  Sv <- estimate_Sv(Sbar, V0, model.name)
                  Sw <- estimate_Sw(Sbar, Sv, V0, model.name)
                  which.maxSb <- which.max(sapply(Sb.list, function(x){max(max(x + (min(min(x)) * diag(m))))}))
                  if (m > 1) {
                    Sw.constr <- check_constraint_SwSb(Sw, Sb.list[[which.maxSb]], model.name)
                    Sw <- Sw.constr$Sw
                    count.SwSb <- Sw.constr$constr.SwSb
                    if (m == p) {
                      Sv <- Sw
                    }
                  }
                  if (m < p) {
                    Sv.constr <- check_constraint_SvSw(Sv, Sw, model.name)
                    Sv <- Sv.constr$Sv
                    count.SvSw <- Sv.constr$constr.SvSw
                  }
                  for (g in 1:G) {
                    St.list[[g]] <- V0 %*% (Sw + Sb.list[[g]]) %*% t(V0) - diag(diag(V0 %*% Sw %*% t(V0))) + diag(diag(V0 %*% Sv %*% t(V0)))
                  }
                }
                max.a <- 0
                for (g in 1:G) {
                  if (min(eigen(St.list[[g]], symmetric = TRUE, only.values = TRUE)$values) <= 0 || det(St.list[[g]]) <= .Machine$double.xmin) {
                    eig.St <- abs(min(eigen(St.list[[g]], symmetric = TRUE, only.values = TRUE)$values)) + sqrt(.Machine$double.eps)
                    if (max.a < eig.St) {
                      max.a <- eig.St
                    }
                  }
                }
                Sv <- Sv + max.a * diag(m)
                St.list[seq(1:G)] <- lapply(St.list, function (x) return(x + max.a * diag(p)))
                ff <-
                  loglikH_EI_G_pugmm(
                    X,
                    model$pp,
                    model$mu,
                    St.list,
                    replicate(G, list(Sv)),
                    replicate(G, list(Sw)),
                    Sb.list,
                    replicate(G, list(V0)),
                    model$post,
                    gaussian
                  )
                if (ff >= fmax) {
                  fmax <- ff
                  posmax <- h
                  Stit <- St.list
                  Svit <- Sv
                  Swit <- Sw
                  Sbit <- Sb.list
                  constr.SwSb <- count.SwSb
                  constr.SvSw <- count.SvSw
                }
              }
            }
            V0[i, ] <- IQ[posmax, ]
          }
          if (any(apply(V0, 2, sum) == 1)) {
            Swit[which(apply(V0, 2, sum) == 1), which(apply(V0, 2, sum) == 1)] <- Svit[which(apply(V0, 2, sum) == 1), which(apply(V0, 2, sum) == 1)]
          }
          model$Sigma <- Stit
          model$V[seq(1:G)] <- list(V0)
          model$Sv[seq(1:G)] <- list(Svit)
          model$Sw[seq(1:G)] <- list(Swit)
          model$Sb <- Sbit
        } else if (model.name == "EEFF") {
          Stit <- model$Sigma
          Svit <- model$Sv[[1]]
          Swit <- model$Sw
          Sbit <- model$Sb
          constr.SwSb <- 0
          constr.SvSw <- 0
          V0 <- model$V[[1]]
          Sw.list <- vector(mode = "list", length = G)
          Sb.list <- vector(mode = "list", length = G)
          St.list <- vector(mode = "list", length = G)
          for (i in 1:p) {
            posmax <- which(V0[i, ] == 1)
            V0[i, ] <- IQ[posmax, ]
            for (h in 1:m) {
              V0[i, ] <- IQ[h, ]
              if (sum(V0[, posmax]) > 0) {
                count.SwSb <- 0
                count.SvSw <- 0
                Sv <- estimate_Sv(Sbar, V0, model.name)
                for (g in 1:G) {
                  Sw <- estimate_Sw(Sigma[[g]], Sv, V0, model.name)
                  Sb <- estimate_Sb(Sigma[[g]], V0, model.name)
                  if (m > 1) {
                    Sw.constr <- check_constraint_SwSb(Sw, Sb, model.name)
                    Sw <- Sw.constr$Sw
                    count.SwSb <- count.SwSb + Sw.constr$constr.SwSb
                  }
                  Sw.list[[g]] <- Sw
                  Sb.list[[g]] <- Sb
                }
                Sw.check <- diag(apply(matrix(unlist(lapply(Sw.list, diag)), m, G), 1, function(x) max(abs(x))), m, m)
                if (m < p) {
                  Sv.constr <- check_constraint_SvSw(Sv, Sw.check, model.name)
                  Sv <- Sv.constr$Sv
                  count.SvSw <- Sv.constr$constr.SvSw
                }
                constr.psd <- 0
                for (g in 1:G) {
                  St.psd <- psd_sigma(V0 %*% (Sw.list[[g]] + Sb.list[[g]]) %*% t(V0) - diag(diag(V0 %*% Sw.list[[g]] %*% t(V0))) + diag(diag(V0 %*% Sv %*% t(V0))))
                  St.list[[g]] <- St.psd$S
                  if (St.psd$constr.psd != 0) {
                    Sb.list[[g]] <- estimate_Sb(St.list[[g]], V0, model.name)
                    constr.psd <- constr.psd + 1
                  }
                }
                if (constr.psd > 0) {
                  Sbar <- apply(sweep(array(unlist(St.list), dim = c(p, p, G)), 3, model$pp, "*"), MARGIN = c(1, 2), sum)
                  count.SwSb <- 0
                  count.SvSw <- 0
                  Sv <- estimate_Sv(Sbar, V0, model.name)
                  for (g in 1:G) {
                    Sw <- estimate_Sw(St.list[[g]], Sv, V0, model.name)
                    if (m > 1) {
                      Sw.constr <- check_constraint_SwSb(Sw, Sb.list[[g]], model.name)
                      Sw <- Sw.constr$Sw
                      count.SwSb <- count.SwSb + Sw.constr$constr.SwSb
                    }
                    Sw.list[[g]] <- Sw
                  }
                  Sw.check <- diag(apply(matrix(unlist(lapply(Sw.list, diag)), m, G), 1, function(x) max(abs(x))), m, m)
                  if (m < p) {
                    Sv.constr <- check_constraint_SvSw(Sv, Sw.check, model.name)
                    Sv <- Sv.constr$Sv
                    count.SvSw <- Sv.constr$constr.SvSw
                  }
                  for (g in 1:G) {
                    St.list[[g]] <- V0 %*% (Sw.list[[g]] + Sb.list[[g]]) %*% t(V0) - diag(diag(V0 %*% Sw.list[[g]] %*% t(V0))) + diag(diag(V0 %*% Sv %*% t(V0)))
                  }
                }
                max.a <- 0
                for (g in 1:G) {
                  if (min(eigen(St.list[[g]], symmetric = TRUE, only.values = TRUE)$values) <= 0 || det(St.list[[g]]) <= .Machine$double.xmin) {
                    eig.St <- abs(min(eigen(St.list[[g]], symmetric = TRUE, only.values = TRUE)$values)) + sqrt(.Machine$double.eps)
                    if (max.a < eig.St) {
                      max.a <- eig.St
                    }
                  }
                }
                Sv <- Sv + max.a * diag(m)
                St.list[seq(1:G)] <- lapply(St.list, function (x) return(x + max.a * diag(p)))
                ff <-
                  loglikH_EI_G_pugmm(
                    X,
                    model$pp,
                    model$mu,
                    St.list,
                    replicate(G, list(Sv)),
                    Sw.list,
                    Sb.list,
                    replicate(G, list(V0)),
                    model$post,
                    gaussian
                  )
                if (ff >= fmax) {
                  fmax <- ff
                  posmax <- h
                  Stit <- St.list
                  Svit <- Sv
                  Swit <- Sw.list
                  Sbit <- Sb.list
                  constr.SwSb <- count.SwSb
                  constr.SvSw <- count.SvSw
                }
              }
            }
            V0[i, ] <- IQ[posmax, ]
          }
          if (any(apply(V0, 2, sum) == 1)) {
            Swit <- lapply(Swit, function(x){x[which(apply(V0, 2, sum) == 1), which(apply(V0, 2, sum) == 1)] <- Svit[which(apply(V0, 2, sum) == 1), which(apply(V0, 2, sum) == 1)]
            return(x)})
          }
          model$Sigma <- Stit
          model$V[seq(1:G)] <- list(V0)
          model$Sv[seq(1:G)] <- list(Svit)
          model$Sw <- Swit
          model$Sb <- Sbit
        } else if (model.name == "EFFF") {
          Stit <- model$Sigma
          Svit <- model$Sv
          Swit <- model$Sw
          Sbit <- model$Sb
          constr.SwSb <- 0
          constr.SvSw <- 0
          V0 <- model$V[[1]]
          Sv.list <- vector(mode = "list", length = G)
          Sw.list <- vector(mode = "list", length = G)
          Sb.list <- vector(mode = "list", length = G)
          St.list <- vector(mode = "list", length = G)
          for (i in 1:p) {
            posmax <- which(V0[i, ] == 1)
            V0[i, ] <- IQ[posmax, ]
            for (h in 1:m) {
              V0[i, ] <- IQ[h, ]
              if (sum(V0[, posmax]) > 0) {
                count.SwSb <- 0
                count.SvSw <- 0
                for (g in 1:G) {
                  Sv <- estimate_Sv(Sigma[[g]], V0, model.name)
                  Sw <- estimate_Sw(Sigma[[g]], Sv, V0, model.name)
                  Sb <- estimate_Sb(Sigma[[g]], V0, model.name)
                  if (m > 1) {
                    Sw.constr <- check_constraint_SwSb(Sw, Sb, model.name)
                    Sw <- Sw.constr$Sw
                    count.SwSb <- count.SwSb + Sw.constr$constr.SwSb
                    if (m == p) {
                      Sv <- Sw
                    }
                  }
                  if (m < p) {
                    Sv.constr <- check_constraint_SvSw(Sv, Sw, model.name)
                    Sv <- Sv.constr$Sv
                    count.SvSw <- count.SvSw + Sv.constr$constr.SvSw
                  }
                  St.psd <- psd_sigma(V0 %*% (Sw + Sb) %*% t(V0) - diag(diag(V0 %*% Sw %*% t(V0))) + diag(diag(V0 %*% Sv %*% t(V0))))
                  St <- St.psd$S
                  if (St.psd$constr.psd != 0) {
                    if (m > 1) {
                    count.SwSb <- count.SwSb - Sw.constr$constr.SwSb
                    }
                    if (m < p) {
                    count.SvSw <- count.SvSw - Sv.constr$constr.SvSw
                    }
                    Sv <- estimate_Sv(St, V0, model.name)
                    Sw <- estimate_Sw(St, Sv, V0, model.name)
                    Sb <- estimate_Sb(St, V0, model.name)
                    if (m > 1) {
                      Sw.constr <- check_constraint_SwSb(Sw, Sb, model.name)
                      Sw <- Sw.constr$Sw
                      count.SwSb <- count.SwSb + Sw.constr$constr.SwSb
                      if (m == p) {
                        Sv <- Sw
                      }
                    }
                    if (m < p) {
                      Sv.constr <- check_constraint_SvSw(Sv, Sw, model.name)
                      Sv <- Sv.constr$Sv
                      count.SvSw <- count.SvSw + Sv.constr$constr.SvSw
                    }
                    St <- V0 %*% (Sw + Sb) %*% t(V0) - diag(diag(V0 %*% Sw %*% t(V0))) + diag(diag(V0 %*% Sv %*% t(V0)))
                  }
                  while (min(eigen(St, symmetric = TRUE, only.values = TRUE)$values) <= 0 || det(St) <= .Machine$double.xmin) {
                    a <- abs(min(eigen(St, symmetric = TRUE, only.values = TRUE)$values)) + sqrt(.Machine$double.eps)
                    St <- St + a * diag(p)
                    Sv <- Sv + a * diag(m)
                  }
                  St.list[[g]] <- St
                  Sv.list[[g]] <- Sv
                  Sw.list[[g]] <- Sw
                  Sb.list[[g]] <- Sb
                }
                ff <-
                  loglikH_EI_G_pugmm(
                    X,
                    model$pp,
                    model$mu,
                    St.list,
                    Sv.list,
                    Sw.list,
                    Sb.list,
                    replicate(G, list(V0)),
                    model$post,
                    gaussian
                  )
                if (ff >= fmax) {
                  fmax <- ff
                  posmax <- h
                  Stit <- St.list
                  Svit <- Sv.list
                  Swit <- Sw.list
                  Sbit <- Sb.list
                  constr.SwSb <- count.SwSb
                  constr.SvSw <- count.SvSw
                }
              }
            }
            V0[i, ] <- IQ[posmax, ]
          }
          if (any(apply(V0, 2, sum) == 1)) {
            Swit <- mapply(function(x, y){x[which(apply(V0, 2, sum) == 1), which(apply(V0, 2, sum) == 1)] <- y[which(apply(V0, 2, sum) == 1), which(apply(V0, 2, sum) == 1)]
            return(x)}, Swit, Svit, SIMPLIFY = FALSE)
          }
          model$Sigma <- Stit
          model$V[seq(1:G)] <- list(V0)
          model$Sv <- Svit
          model$Sw <- Swit
          model$Sb <- Sbit
        }
        loglik.c <-
          loglikH_EI_G_pugmm(
            X,
            model$pp,
            model$mu,
            model$Sigma,
            model$Sv,
            model$Sw,
            model$Sb,
            model$V,
            model$post,
            gaussian
          )
        if (stop == "aitken") {
          if (count == 1) {
            conv <- (loglik.c - loglik)
            loglik.p <- loglik
          } else {
            if (loglik.c > loglik) {
              ait.c <- (loglik.c - loglik) / (loglik - loglik.p)
            } else {
              ait.c <- 0
            }
            loglik.inf <- loglik  + (loglik.c - loglik)/(1 - ait.c)
            conv <- loglik.inf - loglik
            loglik.p <- loglik
          }
        } else if (stop == "relative") {
          conv <- abs(loglik.c - loglik) / abs(loglik)
        }
        count <- count + 1
        loglik <- loglik.c
      }
      if (loop == 1) {
        label.best <- label
        model.best <- model
        loglik.best <- loglik
        constr.SwSb.best <- constr.SwSb
        constr.SvSw.best <- constr.SvSw
        loop.best <- loop
        iter.best <- count - 1
      }
      if (wp == 1) {
        iter.best <- count
      }
      if (loglik > loglik.best) {
        label.best <- label
        model.best <- model
        constr.SwSb.best <- constr.SwSb
        constr.SvSw.best <- constr.SvSw
        loglik.best <- loglik
        loop.best <- loop
        iter.best <- count - 1
      }
    }
    pm.best <- number_param(p, G, m, model.name)
    if (model.name == "FIII" | model.name == "FIIF" | model.name == "FIFF" | model.name == "FFFI" | model.name == "FFFF") {
      free.param <- pm.best - (G * m + constr.SwSb.best + constr.SvSw.best)
    } else {
      free.param <- pm.best - (m + constr.SwSb.best + constr.SvSw.best)
    }
    bic.best <- 2 * loglik.best - (free.param * log(n.obs))

        if (change_model == 1){
      if (model.name == "EUUE"){
        model.name = "EUEE"
      }
      if (model.name == "EEEF"){
        model.name = "EEFF"
      }
      if (model.name == "FIIF"){
        model.name = "FIFF"
      }
    }
    ans <- c(list(call = call,
                  X = X,
                  G = G,
                  m = m,
                  label = label.best,
                  pp = model.best$pp,
                  mu = model.best$mu,
                  sigma = model.best$Sigma,
                  V = model.best$V,
                  Sv = model.best$Sv,
                  Sw = model.best$Sw,
                  Sb = model.best$Sb,
                  post = model.best$post,
                  pm = pm.best,
                  pm.cov = pm.best - (G - 1 + G * p),
                  pm.free = free.param,
                  count.constr.SwSb = constr.SwSb.best,
                  count.constr.SvSw = constr.SvSw.best,
                  bic = bic.best,
                  loglik = loglik.best,
                  loop = loop.best,
                  iter = iter.best,
                  model.name = model.name))

    orderedNames <- c("call","X","G","m","label","pp","mu","sigma","V","Sv","Sw","Sb",
                      "post","pm","pm.cov","pm.free","count.constr.SwSb","count.constr.SvSw",
                      "bic","loglik","loop","iter","model.name")

    structure(ans[orderedNames], class = "PUGMM_big")
  }
