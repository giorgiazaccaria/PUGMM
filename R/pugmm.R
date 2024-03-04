#' Parsimonious Ultrametric Gaussian Mixture Models
#' @description
#' Model-based clustering via Parsimonious Ultrametric Gaussian Mixture Models. Hierarchical relationships among variables within and between clusters are inspected. The grouped coordinate ascent algorithm is used for the parameter estimation. The optimal model is selected according to BIC.
#' @param X (\eqn{n \times p}) numeric matrix or data frame, where \eqn{n} and \eqn{p} represent the number of units and variables, respectively. Categorical variables are not allowed.
#' @param G Integer (vector) specifying the number of mixture components (default: G = 1:5).
#' @param m Integer (vector) specifying the number of variable groups (default: m = 1:5).
#' @param normalization Character string specifying the data transformation. If "NULL", no transformation is applied to the data matrix (default). Other options are: "standard" for the standardization; "center" for centering the data; "range" for the MinMax transformation; "SVD" for the Singular Value Decomposition transformation.
#' @param model Vector of character strings indicating the model names to be fitted. If "NULL", all the possible models are fitted (default). See the possible models using `available_models()`.
#' @param maxiter Integer value specifying the maximum number of iterations of the EM algorithm (default: maxiter = 500).
#' @param tol Numeric value specifying the tolerance for the convergence criteria used in the EM algorithm (default: tol = sqrt(.Machine$double.eps)).
#' @param stop Character string specifying the convergence criteria. If "aitken", the Aitken acceleration-based stopping rule is used (default); if "relative", the relative log-likelihood in two sequential iterations is evaluated.
#' @param rndstart Integer value specifying the number of random starts (default: rndstart = 1).
#' @param initG Character string specifying the method for the initialization of the unit-component membership. If "kmeans", hard k-means via RcppArmadillo is used (default). Other options are: "random" for random assignment; "kmeansf" for fuzzy c-means.
#' @param initm Character string specifying the method for the initialization of the variable-group membership. If "ucms", the Ultrametric Covariance Model is used (default); if "random", random assignment is performed.
#' @param gaussian Character string specifying the way to compute the log-likelihood. If "mclust", `dmvnorm` of `mclust` is used (default); if "canonical", the log-likelihood computation is based upon the canonical representation of an extended ultrametric covariance matrix.
#' @param showprogress A logical value indicating whether the fitting progress should be displayed (default: showprogress = TRUE).
#'
#' @return An object of class `pugmm` containing the results of the optimal - according to BIC - Parsimonious Ultrametric Gaussian Mixture Model estimation. \cr
#' @return `call` Matched call.
#' @return `X` Input data matrix.
#' @return `G` Number of components of the best model.
#' @return `m` Number of variable groups of the best model.
#' @return `label` Integer vector of dimension \eqn{n}, taking values in \eqn{\{1, \ldots, G\}}. It identifies the unit classification according to the maximum a posteriori of the best model.
#' @return `pp` Numeric vector of dimension \eqn{G} containing the prior probabilities for the best model.
#' @return `mu` (\eqn{G \times p}) numeric matrix containing the component mean vectors (by row) for the best model.
#' @return `sigma` List of dimension \eqn{G} containing the (\eqn{p \times p}) numeric component extended ultrametric covariance matrices for the best model.
#' @return `V` List of dimension \eqn{G} containing the (\eqn{p \times m}) binary variable-group membership matrices for the best model.
#' @return `Sv` List of dimension \eqn{G} containing the (\eqn{m \times m}) numeric diagonal matrices of the group variances for the best model.
#' @return `Sw` List of dimension \eqn{G} containing the (\eqn{m \times m}) numeric diagonal matrices of the within-group covariances for the best model.
#' @return `Sb` List of dimension \eqn{G} containing the (\eqn{m \times m}) numeric hallow matrices of the between-group covariances for the best model.
#' @return `post` (\eqn{n \times G}) numeric matrix containing the posterior probabilities for the best model.
#' @return `pm` Number of parameters of the best model.
#' @return `pm.cov` Number of covariance parameters of the best model.
#' @return `pm.free` Number of free parameters of the best model (`pm` - (constraints on \eqn{V} + `count.constr.SwSb` + `count.constr.SvSw`)).
#' @return `count.constr.SwSb` Number of times the constraint between `Sw` and `Sb` has been turned on for the best model.
#' @return `count.constr.SvSw` Number of times the constraint between `Sv` and `Sw` has been turned on for the best model.
#' @return `BIC` BIC values for all the fitted models.
#' @return `bic` BIC value of the best model.
#' @return `loglik` Log-likelihood of the best model.
#' @return `loop` Random start corresponding to the selected solution of the best model.
#' @return `iter` Number of iterations needed to estimate the best model.
#' @return `model.name` Character string denoting the PUGMM model name of the best model among the ones fitted.
#' @details
#' The grouped coordinate ascent algorithm used for the estimation of PUGMMs parameters was demonstrated to be equivalent to an Expectation-Maximization (EM) algorithm in the GMM framework (Hathaway, 1986).
#' @references Hathaway, R. (1986) Another interpretation of the EM algorithm for mixture distributions. \emph{Statistics and Probability Letters}, 4(2), 53-56.
#' @seealso [available_models()], [plot.pugmm()]
#' @examples
#' data(wine, package = "HDclassif")
#' x <- scale(wine[, -1])
#' pugmm.wine <- pugmm(x, 3, 5)
#' table(wine[, 1], pugmm.wine$label)
#' @export
#' pugmm
pugmm <-
  function(X,
           G = NULL,
           m = NULL,
           normalization = NULL,
           model = NULL,
           maxiter = 500,
           tol = sqrt(.Machine$double.eps),
           stop = "aitken",
           rndstart = 1,
           initG = "kmeans",
           initm = "ucms",
           gaussian = "mclust",
           showprogress = TRUE) {

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
    if (!(initG %in% c("random", "kmeansf", "kmeans"))) {
      stop("stop must be 'random' or 'kmeansf' or 'kmeans'.")
    }
    if (!(initm %in% c("random", "ucms"))) {
      stop("stop must be 'random' or 'ucms'.")
    }
    if (!(gaussian %in% c("mclust", "canonical"))) {
      stop("stop must be 'mclust' or 'canonical'.")
    }
    if (!(is.numeric(rndstart) && length(rndstart) == 1)) {
      stop("showprogress must be TRUE or FALSE.")
    }

    call <- match.call()
    num_mod <- length(model)

    all_triples = expand.grid(m, G, 1:num_mod)
    w <- nrow(all_triples)
    gm <- length(G) * length(m)

    # Parallelization
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
        pugmm_case(X,
                   G = all_triples[i, 2],
                   m = all_triples[i, 1],
                   model.name = model[all_triples[i, 3]],
                   model = model,
                   maxiter = maxiter,
                   tol = sqrt(.Machine$double.eps),
                   stop = stop,
                   rndstart = rndstart,
                   initG = initG,
                   initm = initm,
                   gaussian = gaussian,
                   showprogress = showprogress)
      }
    parallel::stopCluster(cl)


    # Warnings
    if (p %in% m) {
      if ("EUEE" %in% model) {
        warning(
          "When the number of variable groups is equal to the number of variables (m = p), the model EUEE is changed into EUUE."
        )
      }
      if ("EEFF" %in% model) {
        warning(
          "When the number of variable groups is equal to the number of variables (m = p), the model EEFF is changed into EEEF."
        )
      }
      if ("FIFF" %in% model) {
        warning(
          "When the number of variable groups is equal to the number of variables (m = p), the model FIFF is changed into FIIF."
        )
      }
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

    ans <- c(
      list(
        call = call,
        X = resbest$X,
        G = resbest$G,
        m = resbest$m,
        label = resbest$label,
        pp = resbest$pp,
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
        model.name = resbest$model.name
      )
    )

    orderedNames <-
      c(
        "call",
        "X",
        "G",
        "m",
        "label",
        "pp",
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
        "model.name"
      )

    structure(ans[orderedNames], class = "pugmm")
  }

