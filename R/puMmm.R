#' @useDynLib PUGMM, .registration = TRUE
#' @import Rcpp
NULL

#' Parsimonious Ultrametric Manly Mixture Models
#' @description
#' Model-based clustering via Parsimonious Ultrametric Manly Mixture Models. Hierarchical relationships among variables within and between clusters are inspected. The grouped coordinate ascent algorithm is used for the parameter estimation. The optimal model is selected according to BIC.
#' @param X (\eqn{n \times p}) numeric matrix or data frame, where \eqn{n} and \eqn{p} represent the number of units and variables, respectively. Categorical variables are not allowed.
#' @param G Integer (vector) specifying the number of mixture components (default: `G = 1:5`).
#' @param m Integer (vector) specifying the number of variable groups (default: `m = 1:5`).
#' @param lambda (\eqn{G \times p}) numeric matrix containing the initial transformation parameters for a single specified G.
#' @param normalization Character string specifying the data transformation. If `NULL`, no transformation is applied to the data matrix (default). Other options are: "standard" for the standardization; "center" for centering the data; "range" for the MinMax transformation; "SVD" for the Singular Value Decomposition transformation.
#' @param model Vector of character strings indicating the model names to be fitted. If `NULL`, all the possible models are fitted (default). See the possible models using `pugmm_available_models()`.
#' @param modelselect Character string indicating the model selection method to be used. If "BIC", the best model is selected according to the BIC (default); if "two-step", the best model is selected according to the two-step model selection method.
#' @param maxiter Integer value specifying the maximum number of iterations of the EM algorithm (default: `maxiter = 500`).
#' @param tol Numeric value specifying the tolerance for the convergence criteria used in the EM algorithm (default: `tol = 1e-6`).
#' @param stop Character string specifying the convergence criteria. If "aitken", the Aitken acceleration-based stopping rule is used (default); if "relative", the relative log-likelihood in two sequential iterations is evaluated.
#' @param rndstart Integer value specifying the number of random starts (default: `rndstart = 1`).
#' @param initG Character string specifying the method for the initialization of the unit-component membership. If "ManlyMix", the `Manly.model()` function via the ManlyMix package is used (default). Other options are: "kmeans" for k-means (via RcppArmadillo); "random" for random assignment; "kmeansf" for fuzzy c-means (via the function fcm of the package ppclust).
#' @param initm Character string specifying the method for the initialization of the variable-group membership. If "ucms", the multivariate model to be used for obtaining the variable-group membership estimated is the same model.name used for estimating the Parsimonious Ultrametric Manly Mixture Model (default); if "random", a random assignment is performed.
#' @param seed Numeric value specifying the seed (default: `seed = 123`).
#' @param parallel A logical value, specifying whether the models should be run in parallel.
#'
#' @return An object of class `puMmm` containing the results of the optimal -according to the model selection criteria- Parsimonious Ultrametric Manly Mixture Model estimation. \cr
#' @return `call` Matched call.
#' @return `X` Input data matrix.
#' @return `G` Number of components of the best model.
#' @return `m` Number of variable groups of the best model.
#' @return `label` Integer vector of dimension \eqn{n}, taking values in \eqn{\{1, \ldots, G\}}. It identifies the unit classification according to the maximum a posteriori of the best model.
#' @return `pp` Numeric vector of dimension \eqn{G} containing the prior probabilities for the best model.
#' @return `lambda` (\eqn{G \times p}) numeric matrix containing the component transformation vectors (by row) for the best model.
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
#' @return `BIC` BIC values for all the fitted models. If BIC is \eqn{NA}, the model has not been computed since its structure is equal to another model, while if BIC is \eqn{-Inf} the solution has a number of clusters \eqn{< G}.
#' @return `bic` BIC value of the best model.
#' @return `loglik` Log-likelihood of the best model.
#' @return `loop` Random start corresponding to the selected solution of the best model.
#' @return `iter` Number of iterations needed to estimate the best model.
#' @return `model.name` Character string denoting the PUGMM model name of the best model among the ones fitted.
#' @return `messages` Messages.
#' @details
#' The grouped coordinate ascent algorithm used for the estimation of PUMMMs parameters was demonstrated to be equivalent to an Expectation-Maximization (EM) algorithm in the GMM framework (Hathaway, 1986).
#' @references Cavicchia, C., Vichi, M., Zaccaria, G. (2024) Parsimonious ultrametric Gaussian mixture models. \emph{Statistics and Computing}, 34, 108.
#' @references Cavicchia, C., Vichi, M., Zaccaria, G. (2022) Gaussian mixture model with an extended ultrametric covariance structure. \emph{Advances in Data Analysis and Classification}, 16(2), 399-427.
#' @references Hathaway, R. (1986) Another interpretation of the EM algorithm for mixture distributions. \emph{Statistics and Probability Letters}, 4(2), 53-56.
#' @seealso [pugmm()], [pugmm_available_models()], [plot.pugmm()]
#' @examples
#' data(Harbour_metals)
#' x <- scale(Harbour_metals[,4:10])
#' \dontrun{
#' results <- puMmm(x, G = 1:4, m = 1:7, model = NULL, modelselect = "two-step")
#' results$G
#' results$m
#' results$model.name
#' table(Harbour_metals$Species, results$label)
#' plot.pugmm(results, what = c("BIC", "Path Diagram"))}
#' @export
#' puMmm
puMmm <-
  function(X,
           G = NULL,
           m = NULL,
           lambda = NULL,
           normalization = NULL,
           model = NULL,
           modelselect = "BIC",
           maxiter = 500,
           tol = 1e-6,
           stop = "aitken",
           rndstart = 1,
           initG = "ManlyMix",
           initm = "ucms",
           seed = 123,
           parallel = FALSE
  ) {

    if(modelselect == "BIC"){
      set.seed(seed)
      best <- main_puMmm(X,
                         G = G,
                         m = m,
                         model = model,
                         maxiter = maxiter,
                         tol = tol,
                         stop = stop,
                         rndstart = rndstart,
                         initG = initG,
                         initm = initm,
                         parallel = parallel)
    } else {
      set.seed(seed)
      step1 <- main_puMmm(X,
                          G = G,
                          m = m,
                          model = "FFFF",
                          maxiter = maxiter,
                          tol = tol,
                          stop = stop,
                          rndstart = rndstart,
                          initG = initG,
                          initm = initm,
                          parallel = parallel)
      set.seed(seed)
      step2 <- puMmm(X,
                     G = step1$G,
                     m = step1$m,
                     model = model,
                     maxiter = maxiter,
                     tol = tol,
                     stop = stop,
                     rndstart = rndstart,
                     initG = initG,
                     initm = initm,
                     parallel = parallel)
      best <- step2
    }

    return(best)
  }
