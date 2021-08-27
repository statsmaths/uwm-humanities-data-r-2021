assert <- function(statement, msg="")
{
  if (!statement)
  {
    stop(msg, call.=(msg==""))
  }
}

ifnull <- function(value, default)
{
  if (is.null(value)) { return(default) }
  
  return(value)
}

cmsg <- function(verbose, fmt, ...)
{
  if (verbose)
  {
    cat(sprintf(fmt, ...))
  }
}

#' Add network metadata and layout information
#'
#' @import        igraph
#'
#' @param edges       data frame describing the edges
#' @param node_name   name of the node id column; defaults to "id"
#' @param directed    logical; is the graph directed?
#'
#' @export
sm_graph_layout <- function(edges, node_name = "id", directed = FALSE)
{
  H <- igraph::graph.edgelist(as.matrix(edges[,1:2]), directed = directed)
  L <- igraph::layout_nicely(H)
  
  vs <- igraph::V(H)
  es <- igraph::get.edgelist(H)
  noms <- names(vs)
  cmp <- igraph::components(H)
  ids <- cbind(match(es[,1], noms),  match(es[,2], noms))
  
  # Properties from the whole graph
  if (directed)
  {
    node_out <- tibble::tibble(
      id = as.character(noms),
      x = L[,1],
      y = L[,2],
      degree_out = igraph::degree(H, mode = "out"),
      degree_in = igraph::degree(H, mode = "in"),
      degree_total = igraph::degree(H, mode = "total"),
      eigen = NA_real_,
      between = NA_real_,
      cluster = as.character(as.integer(
        igraph::membership(igraph::cluster_walktrap(H))
      )),
      component = as.integer(cmp$membership),
      component_size = cmp$csize[as.integer(cmp$membership)]
    )
  } else {
    node_out <- tibble::tibble(
      id = as.character(noms),
      x = L[,1],
      y = L[,2],
      degree = igraph::degree(H, mode = "all"),
      eigen = NA_real_,
      close = NA_real_,
      between = NA_real_,
      cluster = as.character(as.integer(
        igraph::membership(igraph::cluster_walktrap(H))
      )),
      component = as.integer(cmp$membership),
      component_size = cmp$csize[as.integer(cmp$membership)]
    )
  }
  
  names(node_out)[1] <- node_name
  
  # Properties by component
  membership <- as.integer(cmp$membership)
  for (i in unique(membership))
  {
    Hs <- induced_subgraph(H, membership == i)
    index <- which(node_out$component == i)
    node_out$eigen[index] <- igraph::eigen_centrality(Hs, directed = FALSE)$vector
    if (!directed) { node_out$close[index] <- igraph::closeness(Hs) }
    node_out$between[index] <- igraph::betweenness(Hs, directed = directed)
  }
  
  # Reorder components by size
  tab <- as.integer(names(sort(table(membership), decreasing = TRUE)))
  node_out$component <- match(node_out$component, tab)
  
  edge_out <- tibble::tibble(
    x = L[ids[,1],1],
    xend = L[ids[,2],1],
    y = L[ids[,1],2],
    yend = L[ids[,2],2]
  )
  
  list(node = node_out, edge = edge_out)
}

#' Add network metadata to edge and node lists
#'
#' @import        igraph
#'
#' @param edges       data frame describing the edges
#' @param node_name   name of the node id column; defaults to "id"
#'
#' @export
sm_graph_metadata <- function(edges, node_name = "id")
{
  H <- igraph::graph.edgelist(as.matrix(edges[,1:2]), directed = FALSE)
  
  vs <- igraph::V(H)
  es <- igraph::get.edgelist(H)
  noms <- names(vs)
  cmp <- igraph::components(H)
  ids <- cbind(match(es[,1], noms),  match(es[,2], noms))
  
  node_out <- tibble::tibble(
    id = as.character(noms),
    degree = igraph::degree(H, mode = "all"),
    eigen = igraph::eigen_centrality(H, directed = FALSE)$vector,
    close = igraph::closeness(H),
    between = igraph::betweenness(H),
    cluster = as.character(as.integer(
      igraph::membership(igraph::cluster_walktrap(H))
    )),
    component = as.integer(cmp$membership),
    component_size = cmp$csize[as.integer(cmp$membership)]
  )
  names(node_out)[1] <- node_name
  
  node_out
}

sm_text_tfidf <- function(
  object,
  min_df=0.1,
  max_df=0.9,
  max_features=1e4,
  doc_var="doc_id",
  token_var="lemma",
  vocabulary=NULL
) {
  
  assert(inherits(object, "data.frame"), "'input' must be a data frame.")
  assert(doc_var %in% names(object), "no valid 'doc_var' found")
  assert(token_var %in% names(object), "no valid 'token_var' found")
  
  doc_id <- token <- tf <- NULL
  x <- data.frame(
    doc_id = object[[doc_var]],
    token = object[[token_var]],
    stringsAsFactors=FALSE
  )
  
  N <- length(unique(x$doc_id))
  
  if (is.null(vocabulary)) {
    possible_vocab <- table(x[!duplicated(x),]$token) / N
    possible_vocab <- possible_vocab[
      possible_vocab >= min_df & possible_vocab <= max_df
      ]
    possible_vocab <- sort(possible_vocab, decreasing=TRUE)
    vocabulary <- names(possible_vocab[
      seq(1, min(max_features, length(possible_vocab)))
      ])
  }
  
  assert(length(vocabulary) >= 1, "vocabulary length is too small to continue")
  
  # create counts
  x <- x[x$token %in% vocabulary, ]
  tf_tibble <- dplyr::group_by(x, doc_id, token)
  tf_tibble <- dplyr::summarize(tf_tibble, tf = dplyr::n())
  tf_tibble <- dplyr::group_by(tf_tibble, token)
  tf_tibble <- dplyr::mutate(tf_tibble,
                             tfidf = (1 + log2(tf)) * log2(N / dplyr::n())
  )
  tf_tibble <- dplyr::ungroup(tf_tibble)
  
  return(tf_tibble)
}

sm_tidy_matrix <- function(
  x, rows_to = "document", cols_to = "term", values_to = "count"
) {
  if (is.null(rownames(x))) stop("input must have row names")
  if (is.null(colnames(x))) stop("input must have column names")
  
  x <- as.matrix(x)
  out <- tibble::tibble(
    var1 = rownames(x)[row(x)],
    var2 = colnames(x)[col(x)],
    var3 = as.numeric(x)
  )
  
  names(out) <- c(rows_to, cols_to, values_to)
  out
}

sm_tidy_distance <- function(x, item_name = "document")
{
  d <- as.matrix(stats::dist(as.matrix(x)))
  rownames(d) <- rownames(x)
  colnames(d) <- rownames(x)
  
  sm_tidy_matrix(
    d,
    sprintf("%s1", item_name),
    sprintf("%s2", item_name),
    "distance"
  )
}

sm_tidy_angle_distance <- function(x, item_name = "document")
{
  x <- as.matrix(x)
  sim <- x / sqrt(rowSums(x * x))
  sim <- sim %*% t(sim)
  
  out <- sm_tidy_matrix(
    sim, sprintf("%s1", item_name), sprintf("%s2", item_name), "distance"
  )
  out$distance[out$distance > 1] <- 1
  out$distance[out$distance < -1] <- -1
  out$distance <- acos(out$distance) / pi
  out
}

sm_paste <- function(x, name = NULL, sep = "; ") {
  res <- as.data.frame(paste(x, collapse = sep))
  cname <- ifelse(is.null(name), deparse(substitute(x)), name)
  names(res) <- sprintf("%s_paste", cname)
  res
}