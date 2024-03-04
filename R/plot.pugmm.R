#' Plotting method for `pugmm` object
#'
#' @description
#' Plots for Parsimonious Ultrametric Gaussian Mixture Models results, such as BIC and path diagrams.
#'
#' @param x Output from `pugmm`.
#' @param ... Other graphics parameters.
#'
#' @details
#' Types of graphs:
#' \enumerate{
#' \item BIC: plot of BIC values for the fitted models.
#' \item Path Diagram: path diagram representation of the component extended ultrametric covariance matrices for the best model.
#' }
#' @seealso [pugmm()]
#' @examples
#' data(wine, package = "HDclassif")
#' x <- scale(wine[, -1])
#' pugmm.wine <- pugmm(x, 3, 5)
#' plot.pugmm(pugmm.wine)
#' @export plot.pugmm
#' @export
plot.pugmm <-
  function(x,...)
  {
    #what <- match.arg(what, several.ok = TRUE)
    what = c("BIC", "Path Diagram")
    choice <-
      utils::menu(what, graphics = FALSE, title = "Choose a plot")
    while (choice != 0)
    {
      if (what[choice] == "BIC")
      plot.pugmmBIC(x)
    if (what[choice] == "Path Diagram")
      plot.pugmmPathDiagram(x)
    # re-present menu waiting user choice
    choice <-
      utils::menu(what, graphics = FALSE, title = "Choose a plot")
  }
}

plot.pugmmPathDiagram <- function(x,...)
{
  if (x$model.name %in% c("EUUU", "EUUE", "EUEE", "EEEU", "EEEE")) {
    gg = 1
    graphics::par(mfrow = c(1, 1))
  }
  else{
    gg = length(x$V)
    if (gg == 2) {
      graphics::par(mfrow = c(1, 2))
    }

    if (gg == 3 || gg == 4) {
      graphics::par(mfrow = c(2, 2))
    }
    if (gg >= 5) {
      graphics::par(mfrow = c(2, 3))
    }
  }

  number_windows <- ceiling(gg / 6)
  g <- 0


  for (windows in 1:number_windows) {
    if (number_windows != 1) {
      graphics::par(mfrow = c(2, 3))
    }



    g_min <- g+1
    g_max <- min((windows - 1) * 6 + 6, gg)


    for (g in g_min:g_max) {
      p = dim(x$X)[2]
      m = x$m
      V = x$V[[g]]
      Sb = x$Sb[[g]]
      Sw = x$Sw[[g]]
      Sv = x$Sv[[g]]


      ultrcov_sim = ultrcov_sim(Sb, V)
      A = ultrcov_sim$A
      B = ultrcov_sim$B
      levfus = ultrcov_sim$levfus


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

      PathDiagram(p, m, V, outV, A, B, levfus, Sb, Sw, Sv, g)
    }
  }
  graphics::par(mfrow = c(1, 1))
}





plot.pugmmBIC <- function(x,
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

  #args <- list(...)
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
