



g <- und_intnet_chicago$graph
weight <- 'intensity'
node_id1 <- 'V115'
node_id2 <- 'V134'

if(!is.na(weight) && !(weight %in% igraph::edge_attr_names(g))){
  warning("The given weight doens't exist in the edge attributes, using default instead (NA)")
  weight = NA
}

if (is.na(weight)){
  path <- unlist(igraph::get.shortest.paths(g, node_id1, node_id2)$vpath)
}else{
  weight_vector <- igraph::edge_attr(g)[[weight]]
  path_data <- igraph::shortest_paths(graph = g, from = node_id1, to = node_id2, mode = 'all', weights = weight_vector, output = 'both')
  total_weight <- sum(igraph::edge_attr(graph = g, weight, index = igraph::E(g, path = unlist(path_data$vpath))))
}

#igraph::distances(graph = und_intnet_chicago$graph, v = node_id1, to = node_id2, weights = 'weight')
ret <- list(total_weight = total_weight, path = path_data$vpath)