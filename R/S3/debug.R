plot_neighborhood <- function(obj, node_id){
  g <- obj$graph
  events <- obj$events
  w_margin <- 50
  
  nei <- neighbors(g, node_id)
  
  v_coords <- cbind(vertex_attr(g, 'xcoord', nei), vertex_attr(g, 'ycoord', nei))
  v_coords <- rbind(v_coords, c(vertex_attr(g, 'xcoord', node_id), vertex_attr(g, 'ycoord', node_id)))
  colnames(v_coords) <- c('xcoord', 'ycoord')
  rownames(v_coords) <- c(names(nei), node_id)
  
  window_coords <- list(min_x = min(v_coords[, 'xcoord']), min_y = min(v_coords[, 'ycoord']),
                        max_x = max(v_coords[, 'xcoord']), max_y = max(v_coords[, 'ycoord']))
  
  event_coords <- NULL
  for(row in 1:nrow(events)){
    if(events[row,'X'] >= window_coords$min_x - w_margin & 
       events[row,'X'] <= window_coords$max_x + w_margin & 
       events[row,'Y'] >= window_coords$min_y - w_margin & 
       events[row,'Y'] <= window_coords$max_y + w_margin){
      event_coords <- rbind(event_coords, events[row,])
    }
  }
  colnames(event_coords) <- c('xcoord', 'ycoord')
  
  # Plot vertices
  plot(v_coords, xlim = c(window_coords$min_x - w_margin , window_coords$max_x + w_margin), 
       ylim = c(window_coords$min_y - w_margin , window_coords$max_y + w_margin))
  text(x = v_coords[, 'xcoord'], y = v_coords[, 'ycoord'], c(names(nei), node_id), cex=1, col='blue')
  
  # Draw edges
  for(row in 1:nrow(v_coords)-1){
    lines(x = c(v_coords[nrow(v_coords),]['xcoord'], v_coords[row,]['xcoord']),
          y = c(v_coords[nrow(v_coords),]['ycoord'], v_coords[row,]['ycoord']),
          type = "l", lty = 1)
  }
  
  # Plot events
  points(x = event_coords[,'xcoord'], y = event_coords[,'ycoord'], col = 'red')
}

plot_neighborhood(intnet_und, 'V318')