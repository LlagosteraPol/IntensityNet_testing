calculateEventIntensities.intensitynetUnd = function(obj){
  g <- obj$graph
  intensities <- obj$intensities
  edge_counts <- c()
  counts <- c()
  
  pb = txtProgressBar(min = 0, max = gorder(g), initial = 0) 
  
  #TODO: Consider progressive calculation with the other methods
  
  # check if the intensities was previously calculated, if not, calculate them
  if(is.null(intensities)){
    cat("Calculating intensities...\n")
    for(node_id in V(g)){
      if(degree(g, node_id) > 0){
        
        ev_mat <- nodeIntensity(obj, node_id)$intensity
        
        #Adds result of Edgewise intenisty function to 'edge_counts'
        edge_counts[[node_id]] <- ev_mat   
        #Adds result of Nodewise mean intenisty function to 'counts'
        counts[[node_id]] <- Reduce('+', ev_mat)/degree(g, node_id)
        
        
        if(length(counts[[node_id]]) == 0){
          counts[[node_id]] <- 0
        }
      }
      
      # Counts for isolated nodes
      else{
        counts[[node_id]] <- 0
      }
      setTxtProgressBar(pb,node_id)
    }
    close(pb)
    
    intensities <- list(node_mean_intensity = counts, edge_intensity = edge_counts)
    
    g <- g %>% set_vertex_attr(name = "intensity", value = as.matrix(counts)) %>%
      set_edge_attr(name = "intensity", value = as.matrix(edge_counts))
    
    intnet <- list(graph = g, events = obj$events_mtx, graph_type = obj$graph_type, distances = obj$dist_mtx, intensities = intensities)
    attr(intnet, 'class') <- c("intensitynet", "intensitynetUnd")
    return(intnet)
  }
  else{
    print("The intensities of this net have already been calculated.")
    return(obj)
  }
}

edgeIntensity.intensitynet= function(obj,  node_id1, node_id2){
  g <- obj$graph
  distances_mtx <- obj$distances
  events_mtx <- obj$events
  
  # Note that the igraph library already handle the error when one of the node id's 
  # are not part of the graph. Also gives the proper information about it.
  edge_id <- get.edge.ids(g, c(node_id1, node_id2))
  
  # If the intensity of this edge was previously calculated, then return it
  if(edge_id != 0 & !is.null(edge_attr(g, "intensity", index=edge_id))){
    if(length(is.na(vertex_attr(g, "intensity", edge_id)))==0){
      return(edge_attr(g, 'intensity', index=edge_id))
    }
  }
  
  # Distance between the node and its neighbor
  res <- tryCatch(
    {
      abs(distances_mtx[node_id1, node_id2]) # Distance between the node and its neighbor 
    },
    # If the nodes are not part of the graph, give the proper information of the error.
    error=function(cond) {
      if(is.null(match(node_id1, V(g))) & is.null(match(node_id2, V(g)))){
        message("First and second vertices (node_id1, node_id2) doesn't exist in the graph.")
      }else if( is.null(match(node_id1, V(g))) ){
        message("First vertice ID (node_id1) doesn't exist in the graph.")
      }else if( is.null(match(node_id2, V(g))) ){
        message("Second vertice ID (v2) doesn't exist in the graph.")
      }else{
        neighbors_list <- neighbors(g, node_id1)
        if(! V(g)[node_id2] %in% neighbors_list){
          message("Second vertice (node_id2) it's not a neighbor of first vertice (node_id1)")
        }else{
          message(cond)
        }
      }
    }
  )    
  
  node1_coords <- list(xcoord = vertex_attr(g, "xcoord", node_id1),
                       ycoord = vertex_attr(g, "ycoord", node_id1))
  
  node2_coords <- list(xcoord = vertex_attr(g, "xcoord", node_id2),
                       ycoord = vertex_attr(g, "ycoord", node_id2))
  
  # Count events inside a window formed by the node and its neighbor
  
  # Defining event window
  x_min <- min(c(node1_coords$xcoord, node2_coords$xcoord))
  x_max <- max(c(node1_coords$xcoord, node2_coords$xcoord))
  
  y_min <- min(c(node1_coords$ycoord, node2_coords$ycoord))
  y_max <- max(c(node1_coords$ycoord, node2_coords$ycoord))
  
  indicator <- 0
  # Counting events
  for(row in 1:nrow(events_mtx)) {
    event_x <- events_mtx[row, 1]
    event_y <- events_mtx[row, 2]
    
    if(x_min <= event_x & x_max >= event_x & y_min <= event_y & y_max >= event_y){
      indicator <- indicator + 1
    } 
  }
  edge_intensity <- indicator/res
  
  edge_intensity
}

edgeIntensity.intensitynet= function(obj,  node_id1, node_id2, z=100){
  
  if(node_id1 == node_id2){
    stop("Both vertices cannot be the same.")
  }
  
  if(z <= 0){
    print("Warning: 'z' cannot be equal or less than 0, using default.")
    z <- 0
  }
  
  g <- obj$graph
  distances_mtx <- obj$distances
  events_mtx <- obj$events
  
  # Note that the igraph library already handle the error when one of the node id's 
  # are not part of the graph. Also gives the proper information about it.
  edge_id <- get.edge.ids(g, c(node_id1, node_id2))
  
  # If the intensity of this edge was previously calculated, then return it
  if(edge_id != 0 & !is.null(edge_attr(g, "intensity", index=edge_id))){
    if(length(is.na(vertex_attr(g, "intensity", edge_id)))==0){
      return(edge_attr(g, 'intensity', index=edge_id))
    }
  }
  
  # Distance between the node and its neighbor
  res <- tryCatch(
    {
      abs(distances_mtx[node_id1, node_id2]) # Distance between the node and its neighbor 
    },
    # If the nodes are not part of the graph, give the proper information of the error.
    error=function(cond) {
      neighbors_list <- neighbors(g, node_id1)
      if(! V(g)[node_id2] %in% neighbors_list){
        message("Second vertice (node_id2) it's not a neighbor of first vertice (node_id1)")
      }else{
        message(cond)
      }
    }
  )    
  
  node1_coords <- list(xcoord = vertex_attr(g, "xcoord", node_id1),
                       ycoord = vertex_attr(g, "ycoord", node_id1))
  
  node2_coords <- list(xcoord = vertex_attr(g, "xcoord", node_id2),
                       ycoord = vertex_attr(g, "ycoord", node_id2))
  
  # Count events inside a window formed by the node and its neighbor
  
  # Defining event window
  dx <- node1_coords$xcoord - node2_coords$xcoord
  dy <- node1_coords$ycoord - node2_coords$ycoord
  
  if(dx<dy){
    zx <- z/(sqrt(1+(dx/dy)^2))
    zy <- -(dx/dy)*zx
  }else{
    zy <- z/(sqrt(1+(dy/dx)^2))
    zx <- -(dy/dx)*zy
  }
  
  p1 <- c(node1_coords$xcoord - zx, node1_coords$ycoord - zy)
  p2 <- c(node1_coords$xcoord + zx, node1_coords$ycoord + zy)
  p3 <- c(node2_coords$xcoord - zx, node2_coords$ycoord - zy)
  p4 <- c(node2_coords$xcoord + zx, node2_coords$ycoord + zy)
  
  win_points <- list(p1, p2, p3, p4)
  anticlockwise <- orderPoints(x=c(p1[1], p2[1], p3[1], p4[1]), 
                               y=c(p1[2], p2[2], p3[2], p4[2]), clockwise = FALSE)
  
  win_points <- win_points[order(anticlockwise)]
  
  df <- data.frame(x=c(win_points[[1]][1], win_points[[2]][1], win_points[[3]][1], win_points[[4]][1]),
                   y=c(win_points[[1]][2], win_points[[2]][2], win_points[[3]][2], win_points[[4]][2]))
  obj_df <- list(df=df)
  class(obj_df) <- "netTools"
  
  if(clockwise(obj_df)){
    win <- owin(poly=list(x=rev(c(win_points[[1]][1], win_points[[2]][1], win_points[[3]][1], win_points[[4]][1])),
                          y=rev(c(win_points[[1]][2], win_points[[2]][2], win_points[[3]][2], win_points[[4]][2]))))
  }else{
    win <- owin(poly=list(x=c(win_points[[1]][1], win_points[[2]][1], win_points[[3]][1], win_points[[4]][1]),
                          y=c(win_points[[1]][2], win_points[[2]][2], win_points[[3]][2], win_points[[4]][2])))
  }
  
  
  indicator <- 0
  # Counting events
  for(row in 1:nrow(events_mtx)) {
    event_x <- events_mtx[row, 1]
    event_y <- events_mtx[row, 2]
    
    if(inside.owin(event_x, event_y, win)){
      indicator <- indicator + 1
    } 
  }
  edge_intensity <- indicator/res
  
  edge_intensity
}


GeoreferencedPlot.netTools = function(obj){
  g <- obj$graph
  distances_mtx <- obj$distances_mtx
  
  if(!is.null(distances_mtx)){
    norm_coords = layout.norm(matrix(cbind(vertex_attr(g)$xcoord, vertex_attr(g)$ycoord), ncol=2))
    plot(g,
         layout = norm_coords,
         vertex.label=NA,
         vertex.size=2,
         window=TRUE,
         axes=TRUE,
         edge.label = edge_attr(g)$intensity,
         edge.label.cex = 0.5)
  }
  else{
    plot(g, 
         vertex.label=NA, 
         vertex.size=2,
         vertex.size2=2)
  }
}


GeoreferencedGgplot2.netTools = function(obj, ...){
  arguments <- list(...)
  
  g <- obj$graph
  data_df <- obj$data_df
  mode <- obj$mode
  
  node_coords <- data.frame(xcoord = vertex_attr(g)$xcoord, ycoord = vertex_attr(g)$ycoord)
  #rownames(node_coords) <- sprintf("V%s",seq(1:nrow(node_coords)))
  rownames(node_coords) <- vertex_attr(g)$name
  #get edges, which are pairs of node IDs
  edgelist <- get.edgelist(g)
  #convert to a four column edge data frame with source and destination coordinates
  edges <- data.frame(node_coords[edgelist[,1],], node_coords[edgelist[,2],])
  colnames(edges) <- c("xcoord1","ycoord1","xcoord2","ycoord2")
  
  if(is.null(data_df$intensity) || is.na(data_df$heattype)){
    ggplot(data_df, aes(xcoord, ycoord), ...) + 
      geom_point(shape = 19, size = 1.5) +
      geom_segment(aes(x = xcoord1, y = ycoord1, xend = xcoord2, yend = ycoord2), 
                   data = edges, 
                   size = 0.5, 
                   colour = "grey") +
      scale_y_continuous(name = "y-coordinate") + 
      scale_x_continuous(name = "x-coordinate") + theme_bw() +
      theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
    
  }else{
    if(mode == 'moran_i') {
      ggplot(data_df, aes(xcoord, ycoord), ...) + 
        geom_point(aes(colour=as.factor(heattype)), shape=19, size=1.5) +
        geom_tile(aes(fill=as.factor(heattype)), show.legend = FALSE) + 
        labs(title = 'Moran-i Heatmap\n') +
        scale_color_manual(values = c("gray","skyblue", "yellow", "darkorange", "red4"), 
                           name = "", breaks=c(1,2,3,4,5), 
                           labels = c("insignificant","low-low","low-high","high-low","high-high")) +
        geom_segment(aes(x = xcoord1, y = ycoord1, xend = xcoord2, yend = ycoord2), 
                     data = edges, 
                     size = 0.5, 
                     colour = "grey") +
        scale_y_continuous(name = "y-coordinate") + 
        scale_x_continuous(name = "x-coordinate") + theme_bw()
    }else if(mode == 'geary'){
      ggplot(data_df, aes(xcoord, ycoord), ...) + 
        geom_point(aes(colour=as.factor(heattype)), shape=19, size=1.5) +
        geom_tile(aes(fill=as.factor(heattype)), show.legend = FALSE) + 
        labs(title = 'Geary-c Heatmap\n') +
        scale_color_manual(values=c("green", "gray", "red"), 
                           name="", breaks=c(1,2,3), labels=c("positive auto.","no auto.","negative auto.")) +
        geom_segment(aes(x=xcoord1, y=ycoord1, xend = xcoord2, yend = ycoord2), 
                     data=edges, 
                     size = 0.5, 
                     colour="grey") +
        scale_y_continuous(name="y-coordinate") + 
        scale_x_continuous(name="x-coordinate") + theme_bw() +
        theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
      # ggplot(data_df, aes(xcoord,ycoord), ...) +  
      #   geom_point(alpha = 0) + 
      #   geom_tile()+ 
      #   geom_text(aes(label=heattype),hjust=0, vjust=0, size=3, check_overlap = T) +
      #   scale_colour_grey(guide='none') + 
      #   geom_segment(aes(x = xcoord1, y = ycoord1, xend = xcoord2, yend = ycoord2), 
      #                data = edges, 
      #                size = 0.5, 
      #                colour = "grey") +
      #   scale_y_continuous(name="y-coordinate") + 
      #   scale_x_continuous(name="x-coordinate") + theme_bw() 
    }else if(mode=='intensity'){
      ggplot(data_df, aes(xcoord, ycoord, colour = heattype), ...) + 
        geom_point(shape=19, size=1.5,) +
        scale_fill_viridis() +
        labs(title = 'Intensity Heatmap\n', color = 'Norm. intensity') +
        geom_segment(aes(x=xcoord1, y=ycoord1, xend = xcoord2, yend = ycoord2), 
                     data=edges, 
                     size = 0.5, 
                     colour="grey") +
        scale_y_continuous(name="y-coordinate") + 
        scale_x_continuous(name="x-coordinate") + theme_bw() +
        theme(legend.title = element_text(face = "bold"),
              plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
    }
  }
}


CalculateEventIntensities.intensitynetUnd = function(obj){
  g <- obj$graph
  intensities <- obj$intensities
  #edge_counts <- c()
  counts <- c()
  
  # start_time <- Sys.time() # debug only
  # pb = txtProgressBar(min = 0, max = gsize(g), initial = 0) 
  # cat("Calculating edge intensities...\n")
  # for(edge_id in E(g)){
  #   setTxtProgressBar(pb, edge_id)
  #   if(is.null(edge_attr(g, 'intensity', edge_id))){
  #     #Adds result of Edgewise intenisty function to 'edge_counts'
  #     edge_counts[[edge_id]] <- EdgeIntensity(obj, ends(g, edge_id)[1], ends(g, edge_id)[2])
  #   }else if(is.na(edge_attr(g, 'intensity', edge_id))[1]){
  #     edge_counts[[edge_id]] <- 0
  #   }else{
  #     edge_counts[[edge_id]] <- edge_attr(g, 'intensity', edge_id)
  #   }
  # }
  # close(pb)
  # cat(paste0("Time: ", Sys.time() - start_time, "\n")) # debug only
  # 
  # # Encapsulate Edge intensities to pass them to 'MeanNodeIntensity' function to prevent its re-calculation
  # tmp_obj <- SetNetworkAttribute(obj = obj, where = 'edge', name = 'intensity', value = as.matrix(edge_counts))
  
  tmp_obj <- AllEdgeIntensities.intensitynet(obj)
  g <- tmp_obj$graph
  
  pb = txtProgressBar(min = 0, max = gorder(g), initial = 0) 
  cat("Calculating node intensities...\n")
  # check if the intensities was previously calculated, if not, calculate them
  for(node_id in V(g)){
    
    setTxtProgressBar(pb,node_id)
    
    if(is.null(vertex_attr(g, 'intensity', node_id))){
      if(igraph::degree(g, node_id) > 0){
        #Adds result of Nodewise mean intensity function to 'counts'
        counts[[node_id]]  <- MeanNodeIntensity(tmp_obj, node_id)
      }else{
        # Counts for isolated nodes or NA values
        counts[[node_id]] <- 0
      }
    }else if(is.na(vertex_attr(g, 'intensity', node_id))[1]){
      counts[[node_id]] <- 0
    }else{
      counts[[node_id]] <- vertex_attr(g, 'intensity', node_id)
    }
  }
  close(pb)
  
  g <- g %>% set_vertex_attr(name = "intensity", value = as.matrix(counts))
  
  intnet <- list(graph = g, events = obj$events, graph_type = obj$graph_type, distances_mtx = obj$distances_mtx)
  attr(intnet, 'class') <- c("intensitynet", "intensitynetUnd")
  return(intnet)
}


AllEdgeIntensities.intensitynet <- function(obj, z = 5){
  if(z <= 0){
    print("Warning: 'z' cannot be equal or less than 0, using default.")
    z <- 5
  }
  g <- obj$graph
  distances_mtx <- obj$distances_mtx
  event_coords <- obj$events
  edge_list <- get.edgelist(g)
  
  if(length(event_coords) == 0){
    return(NA)
  }
  
  edge_events <- as.numeric(gsub("V", "", edge_list[,1]))
  edge_events <- cbind(edge_events, as.numeric(gsub("V", "", edge_list[,2])))
  edge_events <- cbind(edge_events, 0)
  edge_events <- cbind(edge_events, 0)
  colnames(edge_events) <- c('from', 'to', 'n_events', 'intensity')
  
  node_coords <- as.numeric((gsub("V", "", as_ids(V(g)))))
  node_coords <- cbind(node_coords, vertex_attr(g, "xcoord"))
  node_coords <- cbind(node_coords, vertex_attr(g, "ycoord"))
  colnames(node_coords) <- c('node', 'xcoord', 'ycoord')
  
  start_time <- Sys.time() # debug only
  pb = txtProgressBar(min = 0, max = nrow(event_coords), initial = 0) 
  cat("Calculating edge intensities...\n")
  
  for(row in 1:nrow(event_coords)){
    setTxtProgressBar(pb, row)
    tmp_edge <- NULL
    shortest_d <- NULL
    for(edge_row in 1:nrow(edge_events)){
      # node1 <- node_coords[node_coords[,'node'] == edge_events[edge_row, 'from'],][2:3]
      # node2 <- node_coords[node_coords[,'node'] == edge_events[edge_row, 'to'],][2:3]
      
      # Faster but only works if the node ID is the same as its index
      node1 <- node_coords[edge_events[edge_row, 'from'],][2:3]
      node2 <- node_coords[edge_events[edge_row, 'to'],][2:3]
      
      ep <- event_coords[row, ]
      dist_obj <- list(p1 = node1, p2 = node2, ep = ep)
      class(dist_obj) <- 'netTools'
      d <- PointToSegment(dist_obj)
      
      # If the event is at a distance less or equal 'z' from the edge (segment)
      # connecting both given points (the road), then is counted as an event of that road
      if(d <= z){
        if(d == 0){
          tmp_edge <- edge_row
          break
        }
        if (is.null(shortest_d) || d < shortest_d){
          tmp_edge <- edge_row
        }
      }
    }
    if (!is.null(tmp_edge)){
      edge_events[tmp_edge, 'n_events'] <- edge_events[tmp_edge, 'n_events'] + 1
    }
  }
  close(pb)
  cat(paste0("Time: ", Sys.time() - start_time, "\n")) # debug only
  
  for (edge_row in 1:nrow(edge_events)) {
    # Distance between the node and its neighbor
    edge_dist <- abs(distances_mtx[edge_events[edge_row, 'from'], edge_events[edge_row, 'to']])
    edge_events[edge_row, 'intensity'] <-  edge_events[edge_row, 'n_events'] / edge_dist
  }
  SetNetworkAttribute(obj = obj, 
                      where = 'edge', 
                      name = 'intensity', 
                      value = as.matrix(edge_events[, 'intensity']))
}



# 07/03/22

#' Calculate all the edge intensities of the graph. It's more fast than using iteratively the 
#' function EdgeIntensity for all edges.
#' 
#' @name EdgeIntensitiesAndProportions.intensitynet
#' 
#' @param obj intensitynet object
#' 
#' @return intensitynet class object where the graph contains all the edge intensities as an attribute
#' 
EdgeIntensitiesAndProportions.intensitynet <- function(obj){
  if(obj$event_correction < 0){
    message("Warning: event correction value cannot be less than 0, using default.")
    z <- 5
  }
  else{
    z <- obj$event_correction
  }
  
  g <- obj$graph
  distances_mtx <- obj$distances_mtx
  early_discard_distance <- max(distances_mtx)/2 + z # Set up a limit for an early discard max distance
  event_data <- obj$events
  edge_list <- igraph::ends(g, igraph::E(g), names=FALSE)
  
  if(length(event_data) == 0){
    return(NA)
  }
  
  
  # Create a matrix of distances (medge_event_dist) from each event to each middle point
  # This matrix will be used later to an early event discartion to reduce the computation time
  from_coords <- igraph::vertex_attr(graph = g, name = "xcoord", index = edge_list[,1])
  from_coords <- cbind(from_coords, igraph::vertex_attr(graph = g, name = "ycoord", index = edge_list[,1]))
  
  to_coords <- igraph::vertex_attr(graph = g, name = "xcoord", index = edge_list[,2])
  to_coords <- cbind(to_coords, igraph::vertex_attr(graph = g, name = "ycoord", index = edge_list[,2]))
  
  # Coordinates of the middle point of all edges (vectors)
  edge_mid <- matrix(c( (from_coords[,1] + to_coords[,1]) / 2, (from_coords[,2] + to_coords[,2]) / 2), ncol = 2)
  
  # Distances from each event to each middle point
  medge_event_dist <- proxy::dist(x = as.matrix(event_data[,1:2]), y = edge_mid, method = 'euclidean' )
  medge_event_dist <- `dim<-`(c(medge_event_dist), dim(medge_event_dist))
  
  
  # Prepare structure to set information in the edges
  edge_events <- data.frame(from = edge_list[,1], 
                            to = edge_list[,2],
                            n_events = 0,
                            intensity = 0)
  
  # Set up the names for the covariates of the events
  if(ncol(event_data) > 2){
    for(i in 3:ncol(event_data)){
      if(is.numeric(event_data[,i])){
        edge_events[colnames(event_data[i])] <- 0 # Set event column name
      }else{
        edge_events[as.character(unique(event_data[[i]]))] <- 0 # Set event unique variables as names
      }
    }
  }
  
  node_coords <- as.numeric(igraph::V(g))
  node_coords <- cbind(node_coords, igraph::vertex_attr(g, "xcoord"))
  node_coords <- cbind(node_coords, igraph::vertex_attr(g, "ycoord"))
  colnames(node_coords) <- c('node', 'xcoord', 'ycoord')
  
  #start_time <- Sys.time() # debug only
  pb = utils::txtProgressBar(min = 0, max = nrow(event_data), initial = 0) 
  message("Calculating edge intensities...")
  
  e_count <- 0
  for(row in 1:nrow(event_data)){
    utils::setTxtProgressBar(pb, row)
    tmp_edge <- NULL
    shortest_d <- NULL
    
    ep <- event_data[row, ]
    dist_obj <- list(p1 = from_coords, p2 = to_coords, ep = ep[,1:2])
    class(dist_obj) <- 'netTools'
    event_seg_dist <- PointToSegment(dist_obj)
    
    for(edge_row in 1:nrow(edge_events)){
      
      # Check if the intensities are already calculated
      if(row == 1){
        if(!is.null(igraph::edge_attr(g, 'intensity', igraph::E(g)[edge_row]))){
          e_count <- e_count + 1
          next
        }
      }
      
      # Early event discartion
      if(medge_event_dist[row, edge_row] > early_discard_distance){
        next
      }
      
      # # Faster but only works if the node ID is the same as its index
      # node1 <- node_coords[edge_events[edge_row, 'from'],][2:3]
      # node2 <- node_coords[edge_events[edge_row, 'to'],][2:3]
      # 
      # ep <- event_data[row, ]
      # dist_obj <- list(p1 = node1, p2 = node2, ep = ep[,1:2])
      # class(dist_obj) <- 'netTools'
      # d <- PointToSegment_deprecated(dist_obj)
      
      
      d <- event_seg_dist[edge_row]
      
      
      # If the event is at a distance less or equal 'z' from the edge (segment) which
      # connects both given points (the road), then is counted as an event of that road
      if(d <= z){
        if(d == 0){
          tmp_edge <- edge_row
          break
        }
        if (is.null(shortest_d) || d < shortest_d){
          tmp_edge <- edge_row
        }
      }
    }
    # Set up the information (except intensity) to the edge_event DataFrame
    if ( !is.null(tmp_edge) ){
      edge_events[tmp_edge, 'n_events'] <- edge_events[tmp_edge, 'n_events'] + 1
      
      if(ncol(ep) > 2){
        for(i_col in 3:ncol(ep)){
          if( is.numeric(ep[,i_col]) ){
            tmp_str <- colnames(ep[i_col])
            edge_events[tmp_edge, tmp_str] <- edge_events[tmp_edge, tmp_str] + ep[,i_col]
          }else{
            tmp_str <-  as.character(ep[,i_col])
            edge_events[tmp_edge, tmp_str] <- edge_events[tmp_edge, tmp_str] + 1
          }
        }
      }
    }
    # If the intensity of all edges is already calculated return the object
    if(row == 1){
      if(e_count == length(igraph::E(g))){
        message("Intensities were already calculated.")
        return(obj)
      } 
    }
  }
  close(pb)
  #message(paste0("Time: ", Sys.time() - start_time)) # debug only
  
  #Calculate intensity and proportions
  for (edge_row in 1:nrow(edge_events)) {
    if(edge_events[edge_row, 'n_events'] > 0 ){
      # Distance between the node and its neighbor
      edge_dist <- abs(distances_mtx[edge_events[edge_row, 'from'], edge_events[edge_row, 'to']])
      edge_events[edge_row, 'intensity'] <-  edge_events[edge_row, 'n_events'] / edge_dist
      
      if(ncol(edge_events) > 4){
        for(i_col in 5:ncol(edge_events)){
          edge_events[edge_row, i_col] <- edge_events[edge_row, i_col] / edge_events[edge_row, 'n_events'] 
        }
      }
    }
  }
  
  # Save information from 'edge_events' to the edge attributes of the network
  for(i_col in 3:ncol(edge_events)){
    obj <- SetNetworkAttribute(obj = obj, 
                               where = 'edge', 
                               name = colnames(edge_events[i_col]), 
                               value = as.matrix(edge_events[, i_col]))
  }
  return(obj)
}