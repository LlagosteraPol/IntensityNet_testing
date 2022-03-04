rm(list = ls())

#Set working directory
setwd("R/S3/")

source("./main.R")

# ---------------------------------------------DATA LOADING----------------------------------------------------
# load("../../Data/Castellon.RData")
# load("../../Data/nodes.RData")
# load("../../Data/crimes.RData")
# write.table(Castellon, file="../../Data/castellon.txt") # keeps the rownames
# write.table(nodes, file="../../Data/nodes.txt") # keeps the rownames
# write.table(crimes, file="../../Data/crimes.txt") # keeps the rownames


# Adjacency matrix (undirected): Segmenting locations of the traffic network treated as the vertex set of the network.
castellon <- read.table("../../Data/castellon.txt", header=TRUE, row.names=1) # says first column are rownames

# Node coordinates: Georeferenced coordinates from 'castellon' nodes
nodes <- read.table("../../Data/nodes.txt", header=TRUE, row.names=1) # says first column are rownames

# Event (crime coordinates)
crimes <- read.table("../../Data/crimes_corrected.txt", header=TRUE, row.names=1) # says first column are rownames
crimes_old <- read.table("../../Data/crimes.txt", header=TRUE, row.names=1) # says first column are rownames

#subset of events
crim <- crimes[11:111,] # From crimes, take 11 to 111 (both included)
#crim <- crimes
# --------------------------------------- INIT NETINTENSITY CLASS----------------------------------------------
intnet_und <- intensitynet(castellon, nodes, crim)

castellon_obj <-list(mtx = castellon)
class(castellon_obj) <- "netTools"
dir_castellon <-  Undirected2RandomDirectedAdjMtx(castellon_obj)
intnet_dir <- intensitynet(dir_castellon, nodes, crim, graph_type='directed')

intnet_mix <- intensitynet(dir_castellon, nodes, crim, graph_type='mixed')

# ---------------------------------- NET TOOLS CLASS: FUNCTION TESTING ----------------------------------------

# Choose the type of graph for testing
#intnet <- intnet_und
intnet <- intnet_dir
#intnet <- intnet_mix

class(intnet)

# DISTANCES MATRIX
dist_mtx_test <- intnet$distances

# SHORTEST DISTANCE
short_dist <- ShortestNodeDistance(intnet, node_id1 = 'V101', node_id2 = 'V701')
short_path_int <- ShortestPathIntensity(intnet, node_id1 = 'V101', node_id2 = 'V701')

# GEORREFERENCED PLOT
plot(intnet)

# NODE INTENSITY
node_int <- MeanNodeIntensity(intnet, 'V601') 

# EDGE INTENSITY
edge_int <- EdgeIntensity(intnet, 'V1', 'V9')

# PATH INTENSITY
int_path <- PathIntensity(intnet, short_dist$path)

# All intensities
intnet_all <- RelateEventsToNetwork(intnet)
g <- intnet_all$graph
igraph::edge_attr_names(g)
igraph::vertex_attr_names(g)

if(intnet_all$graph_type == 'undirected'){

  pdf("S3/Plots/area_with_grid_plot.pdf")
  plot(intnet_all, enable_grid = TRUE, axis=TRUE)
  dev.off()
  
  plot(intnet_all, node_label = 'intensity', edge_label='none', vertex.color='red')
  
  igraph::vertex_attr(g, 'intensity', igraph::V(g)['V1']) 
  
  for(node_id in igraph::V(g)){
    if(igraph::V(g)[node_id]$intensity>0) cat(node_id,": ", igraph::V(g)[node_id]$intensity, "\n")
  }
  
  gen_corr <- NodeGeneralCorrelation(intnet_all, dep_type = 'correlation', lag_max = 2, 
                                         intensity = igraph::vertex_attr(g)$intensity)
  
  gen_cov <- NodeGeneralCorrelation(intnet_all, dep_type = 'covariance', lag_max = 2, 
                                         intensity = igraph::vertex_attr(g)$intensity)
  
  data_moran <- NodeLocalCorrelation(intnet_all, dep_type = 'moran_i', intensity = igraph::vertex_attr(g)$intensity)
  moran_i <- data_moran$correlation
  intnet_all <- data_moran$intnet
  
  data_geary <- NodeLocalCorrelation(intnet_all, dep_type = 'geary', intensity = igraph::vertex_attr(g)$intensity)
  geary <- data_geary$correlation
  intnet_all <- data_geary$intnet
  
  # data_getis <- NodeLocalCorrelation(intnet_all, dep_type = 'getis', intensity = igraph::vertex_attr(g)$intensity)
  # getis <- data_getis$correlation
  # intnet_all <- data_getis$intnet
  
} else{
  igraph::vertex_attr(g, 'intensity_in', igraph::V(g)['V1']) 
  igraph::vertex_attr(g, 'intensity_out', igraph::V(g)['V1']) 
  
  pdf("Plots/area_with_grid_Dir_plot.pdf")
  plot(intnet_all, enable_grid = TRUE, vertex_intensity = 'intensity_out')
  dev.off()
  
  pdf("Plots/area_with_grid_Dir_gplot.pdf")
  gplot(intnet_all)
  dev.off()
  
  pdf("Plots/area_with_grid_Dir_moran_gplot.pdf")
  gplot(intnet_all, intensity = igraph::vertex_attr(intnet_all$graph)$intensity_in, heattype = 'moran_i')
  dev.off()
  
  pdf("Plots/area_with_grid_Dir_g_gplot.pdf")
  gplot(intnet_all, intensity = igraph::vertex_attr(intnet_all$graph)$intensity_in, heattype = 'geary')
  dev.off()
  
  for(node_id in igraph::V(g)){
    if(igraph::V(g)[node_id]$intensity_in>0) cat(node_id,": ", igraph::V(g)[node_id]$intensity_in, "\n")
  }
  
  for(node_id in igraph::V(g)){
    if(igraph::V(g)[node_id]$intensity_out>0) cat(node_id,": ", igraph::V(g)[node_id]$intensity_out, "\n")
  }
  
  if(intnet_all$graph_type == 'mixed'){
    igraph::vertex_attr(g, 'intensity_und', igraph::V(g)['V1']) 
    igraph::vertex_attr(g, 'intensity_all', igraph::V(g)['V1']) 
    
    for(node_id in igraph::V(g)){
      if(igraph::V(g)[node_id]$intensity_und>0) cat(node_id,": ", igraph::V(g)[node_id]$intensity_und, "\n")
    }
    
    for(node_id in igraph::V(g)){
      if(igraph::V(g)[node_id]$intensity_all>0) cat(node_id,": ", igraph::V(g)[node_id]$intensity_all, "\n")
    }
  }
  
  correlations <- NodeGeneralCorrelation(intnet_all, dep_type = 'correlation', lag_max = 2, 
                                         intensity = igraph::vertex_attr(g)$intensity_in)
  
  data_moran <- NodeLocalCorrelation(intnet_all, dep_type = 'moran_i', intensity = igraph::vertex_attr(g)$intensity_in)
  moran_i <- data_moran$correlation
  intnet_all <- data_moran$intnet
  
  data_geary <- NodeLocalCorrelation(intnet_all, dep_type = 'geary', intensity = igraph::vertex_attr(g)$intensity_in)
  geary <- data_geary$correlation
  intnet_all <- data_geary$intnet
}

for(edge_id in igraph::E(g)){
  if(igraph::E(g)[edge_id]$intensity>0) print(igraph::E(g)[edge_id]$intensity)
}

#-----------------------------INTENSITYNET CLASS: FUNCTION TESTING---------------------------------




#--------------------------------------------PLOTS------------------------------------------------
pdf("S3/Plots/area_with_grid_Dir.pdf")
plot(intnet_all, enable_grid = TRUE, vertex_intensity = 'intensity_in')
dev.off()

pdf("Plots/area_with_grid.pdf")
plot(intnet_all, enable_grid = FALSE, axis=TRUE)
dev.off()

plot(intnet_all, node_label = 'intensity', edge_label='none', vertex.color='red')

PlotHeatmap(intnet_all)
PlotHeatmap(intnet_all, heattype = 'moran_i')
PlotHeatmap(intnet_all, heattype = 'geary')
PlotHeatmap(intnet_all, heattype = 'v_intensity')
PlotHeatmap(intnet_all, heattype = 'e_intensity')
PlotHeatmap(intnet_all, heattype = 'intensity', net_vertices = igraph::V(g)[1:100])

short_dist <- ShortestNodeDistance(intnet, node_id1 = 'V1', node_id2 = 'V1501')
PlotHeatmap(intnet_all, heattype = 'intensity', net_vertices = short_dist$path)

#-------------------------------------------SUB-NETWORK--------------------------------------------
intnet_und <- intensitynet(castellon, nodes, crimes)
intnet_und_reduced <- ApplyWindow(intnet_und, x_coords = c(752500, 754500), y_coords = c(4429500, 4431500))
intnet_und_reduced <- RelateEventsToNetwork(intnet_und_reduced)

PlotHeatmap(intnet_und_reduced)
PlotHeatmap(intnet_und_reduced, heattype = 'v_intensity')

