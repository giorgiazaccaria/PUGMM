#' Plotting method for `pugmm` object
#'
#' @description
#' Plots for Parsimonious Ultrametric Gaussian Mixture Models results, such as BIC and path diagrams.
#'
#' @param x Output from `pugmm`.
#' @param what A string specifying the type of graph requested. Available choices are: \cr
#' \describe{
#'   \item{\code{"BIC"}}{Plot of BIC values for the fitted models. For each \eqn{G}, the best BIC among the ones corresponding to different \eqn{m} is displayed.}
#'   \item{\code{"Path Diagram"}}{Path diagram representation of the extended ultrametric covariance matrix per component for the best model.}}
#' @param nrow Number of rows in the graphical window. A new graphical window is opened every 6 plots, i.e., components of `pugmm`.
#' @param ncol Number of columns in the graphical window. A new graphical window is opened every 6 plots, i.e., components of `pugmm`.
#' @param cluster_names String of dimension \eqn{G} with the clusters/components' name.
#' @param ... Other graphics parameters.
#' @return No return value since this is a plot method.
#' @seealso [pugmm()]
#' @examples
#' data(penguins)
##' x <- scale(penguins[, 2:5])
#' pugmm.penguins <- pugmm(x, 3, 1)
#' plot.pugmm(pugmm.penguins, what = c("BIC", "Path Diagram"))
#' @export plot.pugmm
#' @export
plot.pugmm <-
  function(x,
           what = NULL,
           nrow = NULL,
           ncol = NULL,
           cluster_names = NULL,
           ...)
  {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(oldpar))
    p <- dim(x$X)[2]
    if (!is.null(cluster_names) & length(cluster_names) != p){
      stop("length of cluster_names must be equal to G, which is the number of components of the best model.")
    }
    # what <- match.arg(what, several.ok = TRUE)
    if (is.null(what)) {
      what = c("BIC", "Path Diagram")
      choice <-
        utils::menu(what, graphics = FALSE, title = "Choose a plot")
      while (choice != 0)
      {
        if (what[choice] == "BIC")
          plot_pugmmBIC(x)
        if (what[choice] == "Path Diagram")
          plot_pugmmPathDiagram(x, nrow, ncol, cluster_names)
        # Represent menu waiting user choice
        choice <-
          utils::menu(what, graphics = FALSE, title = "Choose a plot")
      }
    }
    else{
      if ("BIC" %in% what) {
        plot_pugmmBIC(x)
      }
      if ("Path Diagram" %in% what) {
        plot_pugmmPathDiagram(x, nrow, ncol, cluster_names)
      }
    }
  }

plot_pugmmPathDiagram <- function(x, nrow = NULL, ncol = NULL, cluster_names = NULL, ...)
{
  if (x$model.name %in% c("EUUU", "EUUE", "EUEE", "EEEU", "EEEE")) {
    gg <- 1
    graphics::par(mfrow = c(1, 1))
    nrow <- 1
    ncol <- 1
  }
  else{
    gg <- length(x$V)
    if (is.null(nrow) && is.null(ncol)) {
      if (gg == 1) {
        graphics::par(mfrow = c(1, 1))
        nrow <- 1
        ncol <- 1
      }
      if (gg == 2) {
        graphics::par(mfrow = c(1, 2))
        nrow <- 1
        ncol <- 2
      }
      if (gg == 3 || gg == 4) {
        graphics::par(mfrow = c(2, 2))
        nrow <- 2
        ncol <- 2
      }
      if (gg >= 5) {
        graphics::par(mfrow = c(2, 3))
        nrow <- 2
        ncol <- 3
      }
    }
    if (is.null(nrow) && !is.null(ncol)){
      nrow <- 1
    }
    if (!is.null(nrow) && is.null(ncol)){
      ncol <- 1
    }
  }

  nm <- nrow*ncol
  number_windows <- ceiling(gg / nm)
  g <- 0


  for (windows in 1:number_windows) {
    graphics::par(mfrow = c(nrow, ncol))



    g_min <- g+1
    g_max <- min((windows - 1) * nm + nm, gg)


    for (g in g_min:g_max) {
      p <- dim(x$X)[2]
      m <- x$m
      V <- x$V[[g]]
      Sb <- x$Sb[[g]]
      Sw <- x$Sw[[g]]
      Sv <- x$Sv[[g]]
      column_names <- colnames(x$X)

      if (m == 1) {
        m <- p
        Sv <- as.numeric(Sv)*diag(m)
        Sb <- as.numeric(Sw)*(matrix(1, m, m) - diag(m))
        V <- diag(m)
        ultrcov.sim <- ultrcov_sim(Sb, V)
        Sw <- Sv
      } else if (m > 1) {
        ultrcov.sim <- ultrcov_sim(Sb, V)
      }
      A <- ultrcov.sim$A
      B <- ultrcov.sim$B
      levfus <- ultrcov.sim$levfus


      list <- list()
      for (i in 1:m) {
        list[[i]] <- list(i)
      }
      for (i in 1:length(A)) {
        for (j in 1:length(list)) {
          if (A[i] %in% list[[j]]) {
            z <- j
          }
        }
        for (j in 1:length(list)) {
          if (B[i] %in% list[[j]]) {
            y <- j
          }
        }
        if (B[i] == A[1]) {
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


      PathDiagram(p, m, V, outV, A, B, levfus, Sb, Sw, Sv, g, column_names, cluster_names)
    }
  }
  graphics::par(mfrow = c(1, 1))
}





plot_pugmmBIC <- function(x,
                          G = NULL,
                          modelNames = NULL,
                          symbols = NULL,
                          colors = NULL,
                          xlab = NULL,
                          ylab = "BIC",
                          legendArgs = list(
                            x = "bottomright",
                            ncol = 2,
                            cex = 1,
                            inset = 0.01
                          ),
                          ...)
{
  if (is.infinite(x$bic)) {
    warning("The solution has a number of clusters < G. The user can run models with a reduced number of clusters.")
  }
  else{
    graphics::par(mfrow = c(1, 1))


    # BIC
    BIC <- x$BIC
      z <- dim(BIC)[2]
    modelNames <- dimnames(BIC)[[2]]
    last <- length(dimnames(BIC)[[1]])

    string_first <- dimnames(BIC)[[1]][1]
    string_first <- gsub("[(]", "", string_first)
    string_first <- gsub("[)]", "", string_first)
    string_first <- gsub(",", " ", string_first)
    string_first_split <- strsplit(string_first, " ")
    g_first <- string_first_split[[1]][1]
    m_first <- string_first_split[[1]][2]

    string_last <- dimnames(BIC)[[1]][last]
    string_last <- gsub("[(]", "", string_last)
    string_last <- gsub("[)]", "", string_last)
    string_last <- gsub(",", " ", string_last)
    string_last_split <- strsplit(string_last, " ")
    g_last <- string_last_split[[1]][1]
    m_last <- string_last_split[[1]][2]

    G <- g_first:g_last
    m <- m_first:m_last
    M <- length(m)

    BIC[is.na(BIC)] <- -Inf
    values <- matrix(NA, z, length(G))
    for (i in 1:z) {
      for (j in 1:length(G)) {
        values[i, j] <- max(BIC[((j - 1) * M + 1):(j * M), i])
      }
    }
    values <- t(values)



    x <-
      matrix(
        as.vector(values),
        nrow = nrow(values),
        ncol = ncol(values),
        dimnames = list(G, modelNames)
      )


    n <- ncol(x)
    dnx <- dimnames(x)
    x <- matrix(as.vector(x), ncol = n)
    dimnames(x) <- dnx


    ylim <- range(as.vector(x[is.finite(x)]))

    # args <- list(...)
    if (is.null(xlab))
      xlab <- "G"



    colors <- c(
      "dodgerblue2",
      "red3",
      "green3",
      "slateblue",
      "darkorange",
      "skyblue1",
      "violetred4",
      "forestgreen",
      "steelblue4",
      "slategrey",
      "brown",
      "black",
      "darkseagreen",
      "darkgoldenrod3",
      "olivedrab",
      "royalblue",
      "tomato4",
      "cyan2",
      "springgreen2"
    )

    symbols = c(16, 0, 17, 3, 15, 4, 1, 8, 2, 7,
                5, 9, 6, 10, 11, 18, 12, 13, 14)

    graphics::matplot(
      G,
      x,
      type = "b",
      xaxt = "n",
      xlim = range(G),
      ylim = ylim,
      lty = 1,
      xlab = xlab,
      ylab = ylab,
      main = "",
      col = colors,
      pch = symbols
    )

    graphics::axis(side = 1, at = as.numeric(dnx[[1]]))
    graphics::legend(
      legend = modelNames,
      col = colors,
      pch = symbols,
      x = "bottomright",
      ncol = 2,
      cex = 1,
      inset = 0.01
    )
    if (!is.null(legendArgs))
    {
      do.call("legend", c(
        list(
          legend = modelNames,
          col = colors,
          pch = symbols
        ),
        legendArgs
      ))
    }

    invisible(symbols)
  }
}

