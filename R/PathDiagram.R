PathDiagram <- function (p, m, V, outV, A, B, levfus, Sb, Sw, Sv, g, column_names = NULL, cluster_names = NULL) {
  nnode <- (m * 2) - 1 + p
  fnode <- m + p

  class <- V %*% (1:m)
  s <- matrix(0, 1, p)
  t <- matrix(0, 1, p)
  I <- matrix(0, 1, p)
  sizesq <- 0
  c <- 1
  for (q in 1:m) {
    qq <- outV[q]
    sq <- which(class == qq)
    s[(sizesq + 1):(sizesq + length(sq))] <- qq
    t[(sizesq + 1):(sizesq + length(sq))] <- c
    I[(sizesq + 1):(sizesq + length(sq))] <- sq
    sizesq <- sizesq + length(sq)
    c <- c + 1
  }

  tabd <- cbind(flipudr(A), flipudr(B), (fnode + 1 + m - 2):(fnode + 1))
  tabnode <- matrix(0, m - 1, 2)
  mintab <- matrix(c(min(tabd[m - 1, 1], tabd[m - 1, 2]), tabd[m - 1, 3]), 1, 2)
  tabnode[m - 1, ] <- tabd[m - 1, 1:2]
  if (m > 2) {
    for (q in (m - 2):1) {
      tr <- 0
      ft1 <- which(tabd[q, 1] == mintab[, 1])
      ft2 <- which(tabd[q, 2] == mintab[, 1])
      if (!(length(ft2) == 0)) {
        tabnode[q, 2] <- mintab[ft2, 2]
        mintab[ft2, 2] <- mintab[ft2, 2]
        tr <- tr + 1
      }
      if (!(length(ft1) == 0)) {
        tabnode[q, 1] <- mintab[ft1, 2]
        mintab[ft1, 2] <- max(mintab[ft1, 2], tabd[q, 3])
        tr <- tr + 1
      }
      if (tr == 1) {
        if (length(ft2) == 0) {
          tabnode[q, 2] <- tabd[q, 2]
        } else if (length(ft1) == 0) {
          tabnode[q, 1] <- tabd[q, 1]
          mintab <- rbind(mintab, c(tabd[q, 1], tabd[q, 3]))
        }
      }
      if (tr == 0) {
        mintab <- rbind(mintab, c(min(tabd[q, 1], tabd[q, 2]), tabd[q, 3]))
        tabnode[q, ] <- tabd[q, (1:2)]
      }
    }
  }

  tabnode1 <- matrix(0, m - 1, 2)
  for (qq in 1:m) {
    k = which(tabnode == outV[qq])
    tabnode1[k] <- qq
  }
  for (qq in  1:(m - 1)) {
    for (qqq in 1:2) {
      if (tabnode1[qq, qqq] == 0) {
        tabnode1[qq, qqq] <- tabnode[qq, qqq]
      }
    }
  }

  tt <- ((m + 1):fnode)


  levSv <- t(diag(((diag(
    Sv
  ))) * diag(diag(m))))
  levSv <- diag(levSv)

  outV <- c(outV)
  L <- levSv[outV]
  repV <- colSums(V[, outV])

  levSvf <- L[1] * matrix(1, 1, repV[1])
  for (i in 2:length(repV)) {
    levSvf <- cbind(levSvf, L[i] * matrix(1, 1, repV[i]))
  }

  levSw <- diag(Sw)
  levSw <- levSw[outV]




  smnode <- matrix(3, 1, m)


  # Initialization of the graph G

  G <- igraph::make_empty_graph()
  G <- igraph::add_vertices(G, fnode)
  for (i in 1:length(t)) {
    G <- igraph::add_edges(G, c(t[i], tt[i]))
  }

  kk <- 0
  xinit <- matrix(0, 1, m)
  for (k in 1:m) {
    xinit[1, k] <- (1 + sum(V[, outV[k]] == 1)) / 2 + kk
    kk <- kk + sum(V[, outV[k]] == 1)
  }


  coox <- cbind(xinit, matrix((1:p), 1, p))
  cooy <- cbind(matrix(levSw, 1, length(outV)), levSvf)
  coo <- cbind(t(coox), t(cooy))



  x0data <- cbind((1:m), t(xinit))
  x0data <- rbind(x0data, matrix(0, m - 1, dim(x0data)[2]))

  #Add node from m to 1
  d <- 1
  xxdata <- matrix(0, m - 1, 1)
  ydata <- matrix(0, m - 1, 1)
  #snode <- matrix(0,m-1,1)
  G <- igraph::add_vertices(G, m - 1)
  for (q in (m - 1):1) {
    G <- igraph::add_edges(G, c(m + p + d, tabnode1[q, 1], m + p + d, tabnode1[q, 2]))
    p1 <- which(x0data[, 1] == tabnode1[q, 1])
    p2 <- which(x0data[, 1] == tabnode1[q, 2])
    xxdata[d] <- (x0data[p1, 2] + x0data[p2, 2]) / 2
    ydata[d] <- levfus[d]
    x0data[m + d, 1] <- m + p + d
    x0data[m + d, 2] <- xxdata[d]
    d <- d + 1
  }


  #SizeNodes <- cbind(smnode,matrix(10,1,p),matrix(15,1,d-1))
  coox <- cbind(coox, t(xxdata))
  cooy <- cbind(cooy, t(ydata))
  coo <- cbind(t(coox), t(cooy))
  vlabels <- (1:length(coox))

  coo_min <- min(coo[, 2])
  coo_max <- max(coo[, 2])
  diff <- coo_max - coo_min
  digits <- 6

  diff_round <- round(diff, digits = digits)

  while (diff_round == 0) {
    digits <- digits + 1
    diff_round <- round(diff, digits = digits)
  }

  coo <- round(coo, digits = digits)


  vlabels[1:m] <- rep(NA, m)
  vlabels[(m + 1):(m + p)] <- I
  vlabels[(m + p + 1):(2 * m + p - 1)] <- rep(NA, m - 1)


  DELETE <- 0
  i <- 1
  while (i <= length(igraph::V(G)) - 1 && 0) {
    j <- i + 1
    while (j <= length(igraph::V(G))) {
      if (setequal(coo[i, ], coo[j, ])) {
        DELETE <- DELETE + 1

        mapping <- 1:length(igraph::V(G))
        mapping[i] <- j
        G <- igraph::contract(G, mapping)


        if (!is.na(vlabels[i])) {
          vlabels[j] <- vlabels[i]
        }

        if (!is.na(vlabels[i])) {
          vlabels[i] <- vlabels[j]
        }

        G <- igraph::delete.vertices(G, i)

        coo <- coo[-i, ]
        vlabels <- vlabels[-i]
      }
      j <- j + 1
    }
    i <- i + 1
  }


  G <- igraph::simplify(G)

  delete_one <- 1


  ALIGN <- 0

  vertices_to_delete <- c()


  while (delete_one == 1 ) {
    delete_one <- 0

    G <- igraph::simplify(G)
    yes_children <- c()
    for (i in igraph::V(G)) {
      x <- igraph::neighbors(G, i, mode = "out")
      if (length(x) > 0) {
        yes_children <- c(yes_children, i)
      }
    }


    for (v in yes_children) {
      vertices_to_delete <- c()


      if (v  %in% igraph::V(G)) {
        x <- igraph::neighbors(G, v, mode = "out")
        z <- which((coo[x, 2] == coo[v, 2]))

        if (length(z) > 0){
          for (i in z){
            if (!is.na(vlabels[x[i]])){
              z <- z[z != i]
            }
          }
        }

        if (length(z) > 0) {

          ALIGN <- ALIGN + 1

          delete_one <- 1
          vertices_to_delete <- c(vertices_to_delete, x[z])

          mapping <- 1:length(igraph::V(G))
          mapping[x[z]] <- v
          G <- igraph::contract(G, mapping)

          vertices_to_delete_special <- vertices_to_delete
          if ((length(igraph::V(G)) + 1) %in% vertices_to_delete_special) {
            vertices_to_delete_special <-
              vertices_to_delete[vertices_to_delete != (length(igraph::V(G)) + 1)]
          }

          G <-
            igraph::delete.vertices(G, as.numeric(vertices_to_delete_special))


          if (length(vertices_to_delete) > 0) {
            coo <- coo[-vertices_to_delete,]

            vlabels <- vlabels[-vertices_to_delete]
          }
        }
      }
    }
  }


  coo2 <-  sorting_no_children(G, g, coo, vlabels)


  if (!all(coo2 == coo)) {
    coo <- coo2
  }


  coo2 <- parent_in_the_middle(G, coo)

  if (!all(coo2 == coo)) {
    coo <- coo2
  }


#############################################################################
  selected <- c()

  for (i in igraph::V(G)){
    if (is.na(vlabels[i])){
      selected <- c(selected,i)
    }
  }

  selected_with_label <- c()

  for (i in igraph::V(G)) {
    if (!is.na(vlabels[i])) {
      parent_i <- igraph::neighbors(G, i, mode = c("in"))
      vertices_with_same_parent <- c()
      for (j in selected_with_label) {
        parent_j <- igraph::neighbors(G, j, mode = c("in"))
        if (parent_i == parent_j) {
          vertices_with_same_parent <- c(vertices_with_same_parent,j)
        }

        for (v in vertices_with_same_parent){
          if (coo[v,2] != coo[i,2])
            vertices_with_same_parent <- vertices_with_same_parent[vertices_with_same_parent != v]
        }

      }
      if (length(vertices_with_same_parent) <= 1) {
        selected_with_label <- c(selected_with_label, i)
        selected <- c(selected, i)
      }
      else{
        v_with_min_x <- vertices_with_same_parent[1]
        v_with_max_x <- vertices_with_same_parent[2]
        if (coo[v_with_min_x,1] > coo[v_with_min_x,1]){
          z <- v_with_min_x
          v_with_min_x <- v_with_max_x
          v_with_max_x <- z
        }
        if (coo[v_with_min_x,1] > coo[i,1]){
          selected_with_label <- selected_with_label[selected_with_label != v_with_min_x]
          selected <- selected[selected != v_with_min_x]
          selected_with_label <- c(selected_with_label, i)
          selected <- c(selected, i)
        }
        if (coo[v_with_max_x,1] < coo[i,1]){
          selected_with_label <- selected_with_label[selected_with_label != v_with_max_x]
          selected <- selected[selected != v_with_max_x]
          selected_with_label <- c(selected_with_label, i)
          selected <- c(selected, i)
        }


      }

    }
  }


######################################

    old_root <- length(igraph::V(G))
    no_children <- c()
    for (i in 1:length(igraph::V(G))) {
      x <- igraph::neighbors(G, i, mode = "out")
      if (length(x) == 0) {
        no_children <- c(no_children, i)
      }
    }


    max_y <- max(coo[, 2]) + 0.001

    for (i in no_children) {
      if (i %in% selected_with_label) {
        G <- igraph::add_vertices(G, 1)
        G <- igraph::add_edges(G, c(i, length(igraph::V(G))))
        coo <- rbind(coo, c(coo[i], max_y))
        vlabels <- c(vlabels, vlabels[i])
        selected <- c(selected, length(igraph::V(G)))
      }
    }





####################################

  maybe_overlapping <- get_edge_ids_subgraph(G, selected)

  swap <- 1
  while (swap > 0) {
    new_maybe_overlapping <-  c()
    swap <- 0
    for (i in maybe_overlapping) {
      for (j in maybe_overlapping) {
        if (i < j) {
          ei <- igraph::E(G)[i]
          ej <- igraph::E(G)[j]
          ui <- as.numeric(igraph::tail_of(G, ei))
          vi <- as.numeric(igraph::head_of(G, ei))
          uj <- as.numeric(igraph::tail_of(G, ej))
          vj <- as.numeric(igraph::head_of(G, ej))
          xui <- coo[ui, 1]
          yui <- coo[ui, 2]
          xvi <- coo[vi, 1]
          yvi <- coo[vi, 2]
          xuj <- coo[uj, 1]
          yuj <- coo[uj, 2]
          xvj <- coo[vj, 1]
          yvj <- coo[vj, 2]

          p1 <- c(xui, yui)
          p2 <- c(xvi, yvi)
          q1 <- c(xuj, yuj)
          q2 <- c(xvj, yvj)

          if (!identical(p1,q1) && !identical(p1,q2) && !identical(p2,q1) && !identical(p2,q2)) {
            if (IntersectEdges(p1, p2, q1, q2)) {
              lca <- lca(G, ui, uj)
              path_i <- igraph::all_simple_paths(G,
                                                 from = old_root,
                                                 to = vi,
                                                 #mode = c("in"),
                                                 cutoff = -1)

              x_i <- path_i[[1]][2]

              path_j <- igraph::all_simple_paths(G,
                                                 from = old_root,
                                                 to = vj,
                                                 #mode = c("in"),
                                                 cutoff = -1)
              x_j <- path_j[[1]][2]

              descendant_i <-
                as.numeric(igraph::subcomponent(G, x_i, mode = "out"))
              descendant_j <-
                as.numeric(igraph::subcomponent(G, x_j, mode = "out"))

              x_k <- x_i
              k <- i
              descendant <- descendant_i


              H <- G
              coo_H <- coo

              if (coo[x_j, 2] < coo[x_i, 2]) {
                x_k <- x_j
                k <- j
                descendant <- descendant_j
              }

              delta <- sign(coo[x_k, 1] - coo[lca, 1])

              delta <- delta * length(igraph::V(G))


              for (z in descendant) {
                coo[z, 1] <- coo[z, 1] + delta
              }

              coo <- sorting_no_children(G, g, coo, vlabels)
              swap <- swap + 1

              coo <- parent_vertical(G, coo)


              # Check if they still overlap

              xui <- coo[ui, 1]
              yui <- coo[ui, 2]
              xvi <- coo[vi, 1]
              yvi <- coo[vi, 2]
              xuj <- coo[uj, 1]
              yuj <- coo[uj, 2]
              xvj <- coo[vj, 1]
              yvj <- coo[vj, 2]

              p1 <- c(xui, yui)
              p2 <- c(xvi, yvi)
              q1 <- c(xuj, yuj)
              q2 <- c(xvj, yvj)
              if (IntersectEdges(p1, p2, q1, q2)) {
                G <- H
                coo <- coo_H

                descendant_ui <-
                  as.numeric(igraph::subcomponent(G, ui, mode = "out"))
                descendant_uj <-
                  as.numeric(igraph::subcomponent(G, uj, mode = "out"))

                for (z in descendant_ui) {
                  coo[z, 1] <- coo[z, 1] - (coo[ui, 1] - coo[uj, 1])
                }

                for (z in descendant_uj) {
                  coo[z, 1] <- coo[z, 1] - (coo[uj, 1] - coo[ui, 1])
                }

                sub_i <- igraph::subcomponent(G, ui, mode = "out")
                new_overlapping_i <- get_edge_ids_subgraph(G, sub_i)
                sub_j <- igraph::subcomponent(G, uj, mode = "out")
                new_overlapping_j <- get_edge_ids_subgraph(G, sub_j)
                new_maybe_overlapping <- c(new_maybe_overlapping, new_overlapping_i,
                                           new_overlapping_j)
              }
              else{
                sub_k <- igraph::subcomponent(G, x_k, mode = "out")
                new_overlapping <- get_edge_ids_subgraph(G, sub_k)

                new_maybe_overlapping <- c(new_maybe_overlapping, new_overlapping)
              }
            }
          }
        }
      }
    }

    maybe_overlapping <- new_maybe_overlapping
    maybe_overlapping <- intersect(maybe_overlapping, selected)

  }



#########################################################################

  number_vertices <- length(igraph::V(G))
  for (i in (old_root+1):number_vertices){
    G <- igraph::delete.vertices(G, length(igraph::V(G)))
  }


  element_to_delete <- (old_root+1):number_vertices
  coo <- coo[-element_to_delete, ]
  vlabels <- vlabels[-element_to_delete]

##################################################################



  y_min <- min(coo[, 2])
  y_max <- max(coo[, 2])

  x_min <- min(coo[, 1])
  x_max <- max(coo[, 1])

  y_min_round <- round(y_min, 2)
  y_max_round <- round(y_max, 2)

  round <- 2

  while (y_min_round == y_max_round) {
    round <- round + 1
    y_min_round <- round(y_min, round)
    y_max_round <- round(y_max, round)

  }


  c <- (y_max - y_min) / 10
  G <- igraph::as.undirected(G)

  color <- c()
  size <- c()

  for (i in 1:length(igraph::V(G))) {
    if (is.na(vlabels[i])) {
      color <- c(color, 'orange')
      size <- c(size, 3)
    } else {
      color <- c(color, "#1B9E77")
      size <- c(size, 1)
    }
  }



  a <- min(coo[, 2])
  b <- max(coo[, 2])

  a_round <- round(a, 2)
  b_round <- round(b, 2)

  round <- 2

  while (a_round == b_round) {
    round <- round + 1
    a_round <- round(a, round)
    b_round <- round(b, round)

  }

  a_x <- min(coo[, 1])
  b_x <- max(coo[, 1])



  igraph::plot.igraph(
    G,
    layout = coo,
    vertex.size = size,
    edge.color = "#1B9E77",
    edge.width = 3,
    vertex.frame.width = 0,
    vertex.color = color,
    vertex.label = NA,
    #frame.width = 1,
    xaxt = "n",
    yaxt = "n",
    #ylim = c(-1, 1),
    #xlim = c(-0.8, 0.8),
    #axes = TRUE,
    asp = 1
  )

  no_name <- 0
  for (i in 1:length(coo[,1])) {
    if (is.na(vlabels[i])){
      no_name <- no_name + 1
    }

    if (!is.na(vlabels[i])){

      if(is.null(column_names)){
        label_i <- vlabels[i]
      }
      else{
        label_i <- column_names[i-no_name]

      }




      len_i_word <- nchar(trimws(column_names[i-no_name], "both"))
      y_i <- coo[i,2]
      perc_i <- (y_i-a)/(b-a)

      if (perc_i > 0.9){
        label_i <-  substr(label_i, start = 1, stop = 11)
      }
      if (perc_i > 0.7 && perc_i <= 0.9){
        label_i <-  substr(label_i, start = 1, stop = 15)
      }
      graphics::text(x=-1-2*a_x*(1/(b_x-a_x))+(2/(b_x-a_x))*coo[i,1],
                     y=-1-2*a*(1/(b-a))+(2/(b-a))*coo[i,2]+0.025,
                     labels= label_i, adj = 0,
                     srt=90,xpd = NA,
                     cex = 1-0.07*log2(length(coo[,1])))
    }
  }


  graphics::mtext("Aggregation levels",
                  side = 2,
                  line = 2.7,
                  #cex = 1.2
                  )

  if (!is.null(cluster_names)){
  graphics::mtext(cluster_names[g],
                  side = 1,
                  line = 1,
                  #cex = 3
                 )
  }
  else{
    graphics::mtext(paste("Cluster", g),
                    side = 1,
                    line = 1,
                    #cex = 3
    )
  }

  digits <- 2
  if (abs(a) >= 100 || abs(b) >= 100){
    digits <- 1
  }
  if (abs(a) >= 1000 || abs(b) >= 1000){
    digits <- 0
  }

  graphics::axis(
    2,
    at = seq(-1, 1, by = 2 / 10),
    labels = c(round(seq(a, b, by = -(a - b) / 10), digits)),
    cex.axis = 1.1,
    las = 2,
    line = -1
  )
}


get_edge_ids_subgraph <- function(G, sub) {
  indices_ids <- c()
  for (i in 1:length(igraph::E(G))){
    e_i <- igraph::E(G)[i]
    extremal <- igraph::ends(G, e_i)
    if (extremal[1,1] %in% sub && extremal[1,2] %in% sub){
      indices_ids <- c(indices_ids, i)
    }
  }
  return(indices_ids)
}



IntersectEdges <- function(p1, p2, q1, q2) {
  b <- q1 - p1
  A <-
    cbind(p2 - p1, q1 - q2)

  if (abs(det(A)) < 0.001) {
    return(FALSE)
  }

  x <- solve(A, b)  # x = (s,t)

  if (all(0 < x & x < 1)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}



lca <- function(graph, ...) {
  dots = c(...)
  path = igraph::ego(graph,
             order = length(igraph::V(graph)),
             nodes = dots,
             mode = "in")
  max(Reduce(intersect, path))
}



sorting_no_children <- function(G, g, coo, vlabels) {
  min_x <- 1

  no_children <- c()
  for (i in 1:length(igraph::V(G))) {
    x <- igraph::neighbors(G, i, mode = "out")
    if (length(x) == 0) {
      no_children <- c(no_children, i)
    }
  }


  while (length(no_children) > 0) {
    z <- no_children[1]
    z_coo <- coo[z, 1]

    for (i in no_children) {
      if (coo[i, 1] < coo[z, 1]) {
        z <- i
        z_coo <- coo[z, 1]
      }
    }




    neigh_z <- igraph::neighbors(G, z, mode = "total")
    neigh_neigh_z <- igraph::neighbors(G, neigh_z, mode = "total")


    siblings_z <- c()
    for (i in neigh_neigh_z) {
      if (i %in% no_children && coo[i, 2] == coo[z, 2]) {
        siblings_z <- c(siblings_z, i)
      }
    }


    for (i in siblings_z) {
      coo[i, 1] <- min_x
      min_x <- min_x + 1
    }

    no_children <- no_children[!no_children %in% siblings_z]


  }

  coo2 <- parent_in_the_middle(G, coo)


  if (!all(coo2 == coo)) {
    coo <- coo2
  }

  return(coo)
}



parent_in_the_middle <- function(G, coo) {
  yes_children <- c()
  for (i in 1:length(igraph::V(G))) {
    x <- igraph::neighbors(G, i, mode = "out")
    if (length(x) > 0) {
      yes_children <- c(yes_children, i)
    }
  }

  align <- 1
  count_align <- 0
  while (align  == 1 && count_align < 10) {
    count_align <- count_align + 1
    align <- 0
    for (v in yes_children) {
      children <- igraph::neighbors(G, v, mode = "out")
      mean_x <- mean(coo[children, 1])
      if (mean_x != coo[v, 1]) {
        align <- 1
        coo[v, 1] <- mean_x
      }
    }
  }
  return(coo)
}



parent_vertical <- function(G, coo) {
  for (i in 1:(length(igraph::E(G)) - 1)) {
    for (j in i:length((igraph::E(G)))) {
      ei <- igraph::E(G)[i]
      ej <- igraph::E(G)[j]
      ui <- as.numeric(igraph::tail_of(G, ei))
      vi <- as.numeric(igraph::head_of(G, ei))
      uj <- as.numeric(igraph::tail_of(G, ej))
      vj <- as.numeric(igraph::head_of(G, ej))
      xui <- coo[ui, 1]
      yui <- coo[ui, 2]
      xvi <- coo[vi, 1]
      yvi <- coo[vi, 2]
      xuj <- coo[uj, 1]
      yuj <- coo[uj, 2]
      xvj <- coo[vj, 1]
      yvj <- coo[vj, 2]

      p1 <- c(xui, yui)
      p2 <- c(xvi, yvi)
      q1 <- c(xuj, yuj)
      q2 <- c(xvj, yvj)

      x_values <- c(xui, xvi, xuj, xvj)
      y_values <- c(yui, yvi, yuj, yvj)

      if (!identical(p1,q1) && !identical(p1,q2) && !identical(p2,q1) && !identical(p2,q2)) {
        if (IntersectEdges(p1, p2, q1, q2)) {
          if (yui < yuj) {
            coo[ui, 1] <- coo[vi, 1]
          }
          if (yui > yuj) {
            coo[uj, 1] <- coo[vj, 1]
          }
        }
      }
    }
  }
  return(coo)
}
