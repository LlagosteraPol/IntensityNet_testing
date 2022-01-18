library(intensitynet)
library(spatstat)
data(chicago)
chicago_df <- as.data.frame(chicago[["data"]]) # Get as dataframe the data from Chicago

# Get the adjacency matrix. One way is to create an igraph object from the edge coordinates.
edges <- cbind(chicago[["domain"]][["from"]], chicago[["domain"]][["to"]])
chicago_net <- igraph::graph_from_edgelist(edges)

# And then use the igraph function 'as_adjacency_matrix'
chicago_adj_mtx <- as.matrix(igraph::as_adjacency_matrix(chicago_net))
chicago_node_coords <- data.frame(xcoord = chicago[["domain"]][["vertices"]][["x"]], 
                                  ycoord = chicago[["domain"]][["vertices"]][["y"]])

# Create a dataframe with the coordinates of the events 'assault'
chicago_assault <- chicago_df[chicago_df$marks == 'assault',]
assault_coordinates <- data.frame(xcoord = chicago_assault[,1],
                                  ycoord = chicago_assault[,2])

# Create the intensitynet object, in this case will be undirected 
intnet_chicago <- intensitynet(chicago_adj_mtx, 
                               node_coords = chicago_node_coords, 
                               event_coords = assault_coordinates)

intnet_chicago <- CalculateEventIntensities(intnet_chicago)

data_moran <- NodeLocalCorrelation(intnet_chicago, dep_type = 'moran', intensity = igraph::vertex_attr(intnet_chicago$graph)$intensity)
moran_i <- data_moran$correlation
intnet_chicago <- data_moran$intnet

PlotHeatmap(intnet_chicago, heattype = 'moran')
plot(intnet_chicago, enable_grid = TRUE)