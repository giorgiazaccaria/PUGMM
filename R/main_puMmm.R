main_puMmm <-
  function(X,
           G = NULL,
           m = NULL,
           lambda = NULL,
           normalization = NULL,
           model = NULL,
           maxiter = 500,
           tol = 1e-4,
           stop = "aitken",
           rndstart = 1,
           initG = "ManlyMix",
           initm = "ucms",
           parallel = FALSE
  ) {
    X <- as.matrix(X)
    if (!is.numeric(X)) {
      stop("X must be numeric.")
    }
    n <- dim(X)[1]
    p <- dim(X)[2]
    if (is.null(G)) {
      if (n >= 5) {
        G <- 1:5
      }
      else{
        G <- 1:n
      }
    } else{
      if (!(is.numeric(G) && is.null(dim(G)) && all.equal(G, as.integer(G)) && max(G) < n)) {
        stop("G must be integer values << n.")
      } else{
        G <- sort(unique(G))
      }
    }
    if (is.null(m)) {
      if (p >= 5) {
        m <- 1:5
      }
      else{
        m <- 1:p
      }
    }
    else{
      if (!(is.numeric(m) &&
            is.null(dim(m)) &&
            all.equal(m, as.integer(m)) && max(m) <= p)) {
        stop("m must be integer values included in [1,p].")
      }
      else{
        m <- sort(unique(m))
      }
    }
    if (!is.null(lambda) && (dim(lambda)[1] != dim(X)[2] || dim(lambda)[2] != G)){
      stop("lambda must be of p x G dimension.")
    }
    if (!is.null(normalization)) {
      if (!is.character(normalization) &&
          length(normalization) != 1 &&
          normalization %in%  c("none", "standard", "center", "range", "SVD")) {
        stop("noramlization must be 'none' or 'standard' or 'center' or 'range' or 'SVD'.")
      }
      X <- norm(X, normalization)
    }
    if (is.null(model)) {
      model <-
        c(
          "EUUU",
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
          "FFFF"
        )
    }
    else{
      if (!all(
        model %in% c(
          "EUUU",
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
          "FFFF"
        )
      )) {
        stop(
          "model must be 'EUUU', 'EUUE', 'EUEE', 'EEEU', 'EEEE', 'EEEF', 'EEFF', 'EFFF', 'FIII', 'FIIF', 'FIFF', 'FFFI', 'FFFF'."
        )
      }
      model <- unique(model)
    }
    if (!(is.numeric(maxiter) &&
          length(maxiter) == 1 && all.equal(maxiter, as.integer(maxiter)))) {
      stop(
        "maxiter must be an integer value. We suggest to use maxiter not higher than 500 for computational reason."
      )
    }
    if (!(is.numeric(tol) && length(tol) == 1)) {
      stop("tol must be a numeric value.")
    }
    if (!(stop %in% c("aitken", "relative"))) {
      stop("stop must be 'aitken' or 'relative'.")
    }
    if (!(is.numeric(rndstart) &&
          length(rndstart) == 1 &&
          all.equal(rndstart, as.integer(rndstart)))) {
      stop("rndstart must be an integer value.")
    }
    if (!(initG %in% c("random", "kmeansf", "kmeans", "ManlyMix"))) {
      stop("stop must be 'random' or 'kmeansf' or 'kmeans' or 'ManlyMix'.")
    }
    if (!(initm %in% c("random", "ucms"))) {
      stop("stop must be 'random' or 'ucms'.")
    }
    if (!is.logical(parallel)) {
      stop("parallel must be TRUE or FALSE.")
    }

    call <- match.call()
    num_mod <- length(model)

    all_triples = expand.grid(m, G, 1:num_mod)

    w <- nrow(all_triples)
    gm <- length(G) * length(m)

    # Parallelization
    if (parallel == TRUE) {
      cl <- parallel::makeCluster(parallel::detectCores() - 2)
      doParallel::registerDoParallel(cl)

      res <- list()

      '%dopar%' <- foreach::'%dopar%'

      res <-
        foreach::foreach(
          i = 1:w,
          .inorder = TRUE,
          .export = ".GlobalEnv",
          .packages = c("PUGMM")
        ) %dopar% {
          puMmm_case(
            X,
            G = all_triples[i, 2],
            m = all_triples[i, 1],
            lambda = lambda,
            model.name = model[all_triples[i, 3]],
            model = model,
            maxiter = maxiter,
            tol = tol,
            stop = stop,
            rndstart = rndstart,
            initG = initG,
            initm = initm
          )
        }
      parallel::stopCluster(cl)
    } else {
      res <- list()
      for (i in 1:w) {
        res[[i]] <- puMmm_case(
          X,
          G = all_triples[i, 2],
          m = all_triples[i, 1],
          lambda = lambda,
          model.name = model[all_triples[i, 3]],
          model = model,
          maxiter = maxiter,
          tol = tol,
          stop = stop,
          rndstart = rndstart,
          initG = initG,
          initm = initm
        )
      }
    }


    # Warnings
    one_message <- 0
    if (p %in% m) {
      if ("EUEE" %in% model) {
        one_message <- 1
      }
      if ("EEFF" %in% model) {
        one_message <- 1
      }
      if ("FIFF" %in% model) {
        one_message <- 1
      }
    }
    if (one_message == 0){
      messages <- NULL
    }else{
      messages <- "When the number of variable groups is equal to the number of variables (m = p), the models EUEE, EEFF, FIFF are changed into EUUE, EEEF, FIFF, respectively."
    }

    # BIC
    pars <- list()
    for (i in G) {
      for (j in m) {
        pars[[length(pars) + 1]] <- paste0("(", i, ",", j, ")")
      }
    }

    BIC <- matrix(0, gm, num_mod, dimnames = list(pars, model))
    for (i in 1:num_mod) {
      for (j in 1:gm) {
        if (is.null(res[[(i - 1) * gm + j]]$bic)){
          BIC[j, i] <- NA
        }
        else{
          BIC[j, i] <- res[[(i - 1) * gm + j]]$bic
        }
      }
    }

    # Best index. If more models have similar bic,
    # then we chose the one with lowest pm

    indbicbest <- which(BIC >=  max(BIC, na.rm = TRUE) - tol, arr.ind = TRUE)
    #indbicbest <- which(BIC <= min(BIC, na.rm = TRUE) + tol, arr.ind = TRUE)

    indbest <- matrix(0, dim(indbicbest)[1], 1)
    for (i in 1: dim(indbicbest)[1]) {
      indbest[i] <- (indbicbest[i,2] - 1) * gm + indbicbest[i,1]
    }

    pmbest <- +Inf

    for (i in 1: dim(indbicbest)[1]) {
      if (res[[indbest[i]]]$pm.free < pmbest){
        pmbest <- res[[indbest[i]]]$pm.free
        indbestbest <- indbest[i]
        bic <- res[[indbestbest]]$bic
      }
    }
    indbest <- indbestbest

    p <- dim(X)[2]
    resbest <- res[[indbest]]

    if (bic == -Inf) {
      warning("The solution has a number of clusters < G. The user can run the model with a reduced number of clusters.")
    }

    ans <- c(
      list(
        call = call,
        X = resbest$X,
        Y = resbest$Y,
        G = resbest$G,
        m = resbest$m,
        label = resbest$label,
        pp = resbest$pp,
        lambda = resbest$lambda,
        mu = resbest$mu,
        sigma = resbest$sigma,
        V = resbest$V,
        Sv = resbest$Sv,
        Sw = resbest$Sw,
        Sb = resbest$Sb,
        post = resbest$post,
        pm = resbest$pm,
        pm.cov = resbest$pm.cov,
        pm.free = resbest$pm.free,
        count.constr.SwSb = resbest$count.constr.SwSb,
        count.constr.SvSw = resbest$count.constr.SvSw,
        BIC = BIC,
        bic = bic,
        loglik = resbest$loglik,
        loop = resbest$loop,
        iter = resbest$iter,
        model.name = resbest$model.name,
        messages = messages
      )
    )

    orderedNames <-
      c(
        "call",
        "X",
        "Y",
        "G",
        "m",
        "label",
        "pp",
        "lambda",
        "mu",
        "sigma",
        "V",
        "Sv",
        "Sw",
        "Sb",
        "post",
        "pm",
        "pm.cov",
        "pm.free",
        "count.constr.SwSb",
        "count.constr.SvSw",
        "BIC",
        "bic",
        "loglik",
        "loop",
        "iter",
        "model.name",
        "messages"
      )

    structure(ans[orderedNames], class = "puMmm")
  }
