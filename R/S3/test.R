
# to throw error: stop(sprintf("The object must be of class 'intensitynet', current class is %s", class(tt)[1]))


rm(list = ls())
setwd("R/S3/")
source("./main.R")


# Adjacency matrix (undirected): Segmenting locations of the traffic network treated as the vertex set of the network.
load("../../Data/Castellon.RData")

# Node coordinates: Georeferenced coordinates from 'castellon' nodes
load("../../Data/nodes.RData")

# Event (crime coordinates)
load("../../Data/crimes.RData")


#subset of events
crim <- crimes[11:111,] # From crimes, take 11 to 111 (both included)


init_graph <- function(obj){
  weighted_mtx = obj$adjacency_mtx * obj$distances_mtx
  g <- graph_from_adjacency_matrix(weighted_mtx, mode = obj$graph_type, weighted = TRUE)
  
  net_coords <- list(graph = g, node_coords = obj$node_coords)
  class(net_coords) <- "netTools"
  g <- setNetCoords(net_coords)
  
  g # return
}

setNetCoords <- function(obj){
  UseMethod("setNetCoords")
}

#TODO: Declare non-visible at namespace
setNetCoords.netTools = function(obj){
  x_coord_node <- obj$node_coords[1]
  y_coord_node <- obj$node_coords[2]
  
  # TODO: change x and y to long (longitude) and lat (latitude), beware of the coordinate system type (e.g WGS84)
  g <- obj$g %>% set_vertex_attr(name = "xcoord", value = x_coord_node) %>% 
    set_vertex_attr(name = "ycoord", value = y_coord_node)
  g
}

# ----------------------TESTING 'shortestDistance'-------------------------

shortestDistance = function(g, node_id1, node_id2, distances_mtx){
  
  weighted_path <- unlist(get.shortest.paths(g, node_id1, node_id2)$vpath)
  if(!is.null(distances_mtx)){
    weight_sum <- sum(E(g, path = unlist(weighted_path))$weight)
  }
  else{
    weight_sum <- length(weighted_path)
  }
  list(path = weighted_path, weight = weight_sum)  
}

nn <- as.matrix(nodes)
node_coords_obj <- list(node_coords=nn)
class(node_coords_obj) <- "netTools"
dist_mtx <- calculateDistancesMtx(node_coords_obj)

net_setup <- list(adjacency_mtx = as.matrix(Castellon), 
                  node_coords = nn, 
                  events = as.matrix(crimes), 
                  distances_mtx = dist_mtx, 
                  graph_type = 'undirected')

g <- init_graph(net_setup)
intnet <- intensitynet(Castellon, nodes, crim)
shortest <- shortestDistance(intnet$graph, 'V1', 'V20', intnet$distances_mtx)

# ---------------Setting vertice attributes-----------------------
vertex_attr(g, 'intensity', index='V1') 
g <- g %>% set_vertex_attr("intensity", index='V1', value = 333)
vertex_attr(g, 'intensity','V2') 

if(!is.null(vertex_attr(g, "intensity", index='V2'))){
  if(!is.na(vertex_attr(g, "intensity", index='V2')))  vertex_attr(g, "intensity", 'V1')
}

g <- g %>% set_edge_attr(name = "test", index = E(g)[1,2,3], value = c(10,20,30))
class(edge_attr(g)[1])
edge_attr(g)[2]
edge_attr_names(g)
edge_attr(g, "test", E(gt)[1,2,3])


library(sna)
g <- intnet_all$graph
nacf(g, vertex_attr(g, 'intensity'))


# ----------------------------------Window between two points--------------------
node1 <- c(2, 5)
node2 <- c(7, 1)
z <- 0.5

dx <- node1[1] - node2[1]
dy <- node1[2] - node2[2]

if(dx<dy){
  zx <- z/(sqrt(1+(dx/dy)^2))
  zy <- -(dx/dy)*zx
}else{
  zy <- z/(sqrt(1+(dy/dx)^2))
  zx <- -(dy/dx)*zy
}

p1 <- c(node1[1] - zx, node1[2] - zy)
p2 <- c(node1[1] + zx, node1[2] + zy)
p3 <- c(node2[1] - zx, node2[2] - zy)
p4 <- c(node2[1] + zx, node2[2] + zy)

plot(1, type = "n",                         # Remove all elements of plot
     xlab = "", ylab = "",
     xlim = c(749000, 750000), ylim = c(4430500, 4431000))
lines(c(node1[1], node2[1]), c(node1[2], node2[2]), type = "l", lty = 1)
points(c(p1[1], p2[1], p3[1], p4[1]), c(p1[2], p2[2], p3[2], p4[2]))

lines(c(p1[1], p3[1]), c(p1[2], p3[2]), type = "l", lty = 1)
lines(c(p2[1], p4[1]), c(p2[2], p4[2]), type = "l", lty = 1)

lines(c(p1[1], p2[1]), c(p1[2], p2[2]), type = "l", lty = 1)
lines(c(p3[1], p4[1]), c(p3[2], p4[2]), type = "l", lty = 1)

library(contoureR)
win_points <- list(p1, p2, p3, p4)
anticlockwise <- orderPoints(x=c(p1[1], p2[1], p3[1], p4[1]), 
                             y=c(p1[2], p2[2], p3[2], p4[2]), clockwise = FALSE)

win_points <- win_points[order(anticlockwise)]

win <- tryCatch(
  {
    owin(poly=list(x=c(win_points[[1]][1], win_points[[2]][1], win_points[[3]][1], win_points[[4]][1]),
                   y=c(win_points[[1]][2], win_points[[2]][2], win_points[[3]][2], win_points[[4]][2])))
  },
  error=function(cond) {
    owin(poly=list(x=rev(c(win_points[[1]][1], win_points[[2]][1], win_points[[3]][1], win_points[[4]][1])),
                   y=rev(c(win_points[[1]][2], win_points[[2]][2], win_points[[3]][2], win_points[[4]][2]))))
  }
)


plot(win)


df <- data.frame(x=c(win_points[[1]][1], win_points[[2]][1], win_points[[3]][1], win_points[[4]][1]),
                 y=c(win_points[[1]][2], win_points[[2]][2], win_points[[3]][2], win_points[[4]][2]))
clockwise(df)
clockwise <- function(x) {
  
  x.coords <- c(x[[1]], x[[1]][1])
  y.coords <- c(x[[2]], x[[2]][1])
  
  double.area <- sum(sapply(2:length(x.coords), function(i) {
    (x.coords[i] - x.coords[i-1])*(y.coords[i] + y.coords[i-1])
  }))
  
  double.area > 0
} 

dist_obj <- list(p1= node1, p2= node2, ep=c(3,3))
class(dist_obj) <- 'netTools'
perpenidularDistance(dist_obj)

#-----------------------------SNA library testings------------------------------------


intnet <- Intensitynet(Castellon, nodes, crim)
intnet_all <- RelateEventsToNetwork(intnet)
g <- intnet_all$graph


event_correlation <- function(g, dep_type, lag_max){
  g_sna <- intergraph::asNetwork(g)
  nacf(g_sna, vertex_attr(g, "intensity"), type = dep_type, mode = "graph", lag.max = lag_max)
}

# Manual graph set-up
#test_g <- igraph::graph(c(1,2, 1,3, 2,3, 3,4, 3,5), directed = FALSE)
test_g <- igraph::graph(c(1,3, 2,3, 3,4, 3,5), directed = FALSE)
plot(test_g)
test_g <- test_g %>% set_vertex_attr(name = "intensity", value = cbind(c(30, 25, 27, 20, 22)))
test_g <- test_g %>% set_edge_attr(name = "intensity", index=E(test_g), value = cbind(c(10, 15, 20, 25)))


adj_mtx <- as_adj(graph = test_g, attr = 'intensity')
m_adj_listw <- mat2listw(adj_mtx, style="M")
w_adj_listw <- mat2listw(adj_mtx, style="W")
b_adj_listw <- mat2listw(adj_mtx, style="B")

nb <- m_adj_listw$neighbours
w_listw <- nb2listw(nb, style="W", zero.policy=T) 
b_listw <- nb2listw(nb, style="B", zero.policy=TRUE) 

# Moran I
auto <- event_correlation(test_g, 'moran', 2)
gen_moran <- moran(x = vertex_attr(test_g)$intensity, listw = w_listw, n=length(nb), S0=Szero(w_listw))

# Local Moran I
node_locmoran <- localmoran(x = vertex_attr(test_g)$intensity, listw = w_listw, zero.policy=FALSE, na.action = na.omit)
#edge_locmoran <- localmoran(x = edge_attr(test_g)$intensity, listw = w_listw, zero.policy=FALSE, na.action = na.omit)

# Local Moran I using a sub-graph
sub_test_g <- induced_subgraph(test_g, 1:3)
sub_adj_mtx <- as_adj(graph = sub_test_g, attr = 'intensity')
sub_m_adj_listw <- mat2listw(sub_adj_mtx, style="M")
sub_nb <- sub_m_adj_listw$neighbours
sub_w_listw <- nb2listw(sub_nb, style="W", zero.policy=T) 

sub_auto <- event_correlation(sub_test_g, 'moran', 2)

sub_node_locmoran <- localmoran(x = vertex_attr(sub_test_g)$intensity, listw = sub_w_listw, zero.policy=FALSE, na.action = na.omit)
sub_edge_locmoran <- localmoran(x = edge_attr(sub_test_g)$intensity, listw = sub_w_listw, zero.policy=FALSE, na.action = na.omit)


#-------------------------------------------------Matrices testing-------------------------------------------------
test_mtx <- matrix(data=c(2,4,6,8, 10,12,14,16, 18,20,22,24, 0,28,30,32), nrow = 4)
diag_test_mtx <- diag(test_mtx)
upp_tri_test_mtx <- upper.tri(test_mtx)
#upp_tri_test_mtx[lower.tri(upp_tri_test_mtx)] <- 0

low_tri_test_mtx <- test_mtx
low_tri_test_mtx[upper.tri(low_tri_test_mtx)] <- 0

tmp_mtx <- test_mtx #1 * upp_tri_test_mtx
for(row in 1:nrow(tmp_mtx)) {
  for(col in row:ncol(tmp_mtx)) {
    if(tmp_mtx[row, col] != 0){
      if(runif(1, 0.0, 1.0) <= 0.15) tmp_mtx[row, col] <- 0
    }
  }
}
random_directed_matrix <- tmp_mtx

#------------------------------------------------------Ploting------------------------------------------------------
node_coords <- matrix(cbind(vertex_attr(g)$xcoord, vertex_attr(g)$ycoord), ncol=2)

min_x <- min(node_coords[,1])
max_x <- max(node_coords[,1])
min_y <- min(node_coords[,2])
max_y <- max(node_coords[,2])


x_dist <- max_x - min_x
y_dist <- max_y - min_y

plot(g, layout=node_coords, vertex.label=NA, vertex.size=2, window=TRUE, axes=FALSE, 
     edge.label = '', edge.label.cex = 0.2)

rect(-1.05,-1.05,-1.05,1.05)
rect(-1.05,-1.05, 1.05,-1.05)
rect(-1.05,1.05, 1.05,1.05)
rect(1.05,-1.05, 1.05,1.05)

mtext(expression(bold("x-coordinate")), at=0, line = -28, cex=0.70)
mtext(floor(min_x + (1 * (x_dist/6))), at=-0.67, line = -27, cex=0.70)
mtext(floor(min_x + (2 * (x_dist/6))), at=-0.34, line = -27, cex=0.70)
mtext(floor(min_x + (3 * (x_dist/6))), at=0, line = -27, cex=0.70)
mtext(floor(min_x + (4 * (x_dist/6))), at=0.34, line = -27, cex=0.70)
mtext(floor(min_x + (5 * (x_dist/6))), at=0.67, line = -27, cex=0.70)

mtext(expression(bold("y-coordinate")), at=0, line=-15, cex=0.70, side = 2)
mtext(floor(min_y + (1 * (y_dist/6))), at=-0.67, line = -16, cex=0.70, side = 2)
mtext(floor(min_y + (2 * (y_dist/6))), at=-0.34, line = -16, cex=0.70, side = 2)
mtext(floor(min_y + (3 * (y_dist/6))), at=0, line = -16, cex=0.70, side = 2)
mtext(floor(min_y + (4 * (y_dist/6))), at=0.34, line = -16, cex=0.70, side = 2)
mtext(floor(min_y + (5 * (y_dist/6))), at=0.67, line = -16, cex=0.70, side = 2)


rect(min(node_coords[,1]),min(node_coords[,2]), max(node_coords[,1]),min(node_coords[,2]),border="black",lwd=1,col="black")
rect(min(node_coords[,1]),min(node_coords[,2]), min(node_coords[,1]),min(node_coords[,2]),border="black",lwd=1,col="black")
rect(max(node_coords[,1]),min(node_coords[,2]), max(node_coords[,1]),max(node_coords[,2]),border="black",lwd=1,col="black")
rect(min(node_coords[,1]),max(node_coords[,2]), max(node_coords[,1]),max(node_coords[,2]),border="black",lwd=1,col="black")


#--------------------------------------------------ggplot2-----------------------------------------------------

x <- c(1,3,4,6,9)
y <- x^2
coords = paste(x,y,sep=",")

df = data.frame(x,y)

ggplot(df,aes(x,y))+geom_point(col="blue")+
  geom_label(aes(x+.5,y+0.5,label=coords))


#-------------------------------------------------visNetwork----------------------------------------------

intnet <- intensitynet(Castellon, nodes, crim)
intnet <- RelateEventsToNetwork(intnet)
intnet <- NodeLocalCorrelation(intnet, 'moran')
intnet <- NodeLocalCorrelation(intnet, 'g')

g <- intnet$graph

nodes <- data.frame(id = paste(vertex_attr(g)$name),
                    xcoord = vertex_attr(g)$xcoord,
                    ycoord = vertex_attr(g)$ycoord,
                    intensity = vertex_attr(g)$intensity)

edges <- data.frame(from = get.edgelist(g)[,1],
                    to = get.edgelist(g)[,2],
                    distance = edge_attr(g)$weight,
                    intensity = edge_attr(g)$intensity)

nodes <- data.frame(id = paste(vertex_attr(g)$name),
                    x = vertex_attr(g)$xcoord,
                    y = vertex_attr(g)$ycoord,
                    label = paste(round(vertex_attr(g)$intensity, 4)))

edges <- data.frame(from = get.edgelist(g)[,1],
                    to = get.edgelist(g)[,2],
                    label = paste(round(edge_attr(g)$intensity, 4)))


visNetwork(nodes, edges)  %>%
  visIgraphLayout() %>%
  visEvents(selectNode = "function(properties) {
      alert('Node Properties: ' + this.body.data.nodes.get(properties.nodes[0]).id);}") %>% 
  visOptions(nodesIdSelection = TRUE,  highlightNearest = TRUE)


tt <- paste(round(vertex_attr(g)$getis_g, 4))
nodes <- data.frame(id = paste(vertex_attr(g)$name),
                    label = paste(round(vertex_attr(g)$getis_g, 4)))

#------------------------------------------Coloring plot-------------------------------------------------------
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("graph")

intnet <- intensitynet(Castellon, nodes, crim)
intnet <- RelateEventsToNetwork(intnet)
intnet <- NodeLocalCorrelation(intnet, 'moran')
intnet <- NodeLocalCorrelation(intnet, 'g')

g <- intnet$graph


library("graph")
g1  <- igraph.to.graphNEL(intnet$graph) 
adj_mtx1 <- as(g1, "sparseMatrix") 
nb1 <- mat2listw(adj_mtx1)$neighbours
w_listw1 <- nb2listw(nb1, style="W", zero.policy=T) 

g <- intnet$graph
adj_mtx <- as_adj(graph = g)
adj_listw <- mat2listw(adj_mtx)
nb <- adj_listw$neighbours
w_listw <- nb2listw(nb, style="W", zero.policy=TRUE) 
locmoran <- localmoran(x = vertex_attr(g)$intensity, listw = w_listw, zero.policy=TRUE)
summary(locmoran)

node_coords <- data.frame(xcoord = vertex_attr(g)$xcoord, ycoord = vertex_attr(g)$ycoord)
rownames(node_coords) <- sprintf("V%s",seq(1:nrow(node_coords)))
#get edges, which are pairs of node IDs
edgelist <- get.edgelist(g)
#convert to a four column edge data frame with source and destination coordinates
edges <- data.frame(node_coords[edgelist[,1],], node_coords[edgelist[,2],])
colnames(edges) <- c("xcoord1","ycoord1","xcoord2","ycoord2")

# Calculate deviations
node_int_deviation <- vertex_attr(g)$intensity - mean(vertex_attr(g)$intensity)  
locmoran_deviation <- locmoran[,1] - mean(locmoran[,1])  # get the 'li' component 

# create a new variable identifying the moran plot quadrant for each observation, dismissing the non-significant ones
quad_sig <- NA
significance <- 0.5

# non-significant observations
quad_sig[(locmoran[, 5] > significance)] <- 0 # "insignificant"  
# low-low quadrant
quad_sig[(node_int_deviation < 0 & locmoran_deviation < 0) & (locmoran[, 5] <= significance)] <- 1 # "low-low"
# low-high quadrant
quad_sig[(node_int_deviation < 0 & locmoran_deviation > 0) & (locmoran[, 5] <= significance)] <- 2 # "low-high"
# high-low quadrant
quad_sig[(node_int_deviation > 0 & locmoran_deviation < 0) & (locmoran[, 5] <= significance)] <- 3 # "high-low"
# high-high quadrant
quad_sig[(node_int_deviation > 0 & locmoran_deviation > 0) & (locmoran[, 5] <= significance)] <- 4 # "high-high"


length(which(quad_sig==0)) # To count values

data_df <- data.frame(intensity = vertex_attr(g)$intensity , 
                      xcoord = node_coords$xcoord, 
                      ycoord = node_coords$ycoord, 
                      heattype = NA)

data_df_i <- data.frame(intensity = vertex_attr(g)$intensity , 
                        xcoord = node_coords$xcoord, 
                        ycoord = node_coords$ycoord, 
                        heattype = quad_sig)

b_listw <- nb2listw(nb, style="B", zero.policy=TRUE) 
# local net G
locg_all <- localG(x = vertex_attr(g)$intensity, listw = b_listw)
locg <- unlist(as.list(round(locg_all, 1)))

data_df_g <- data.frame(intensity = vertex_attr(g)$intensity, 
                      xcoord = node_coords$xcoord, 
                      ycoord = node_coords$ycoord, 
                      heattype = locg)


ggplot(data_df, aes(xcoord,ycoord)) + 
  geom_point(shape=19, size=1.5) + #show.legend = FALSE
  geom_segment(aes(x=xcoord1, y=ycoord1, xend = xcoord2, yend = ycoord2), 
               data=edges, 
               size = 0.5, 
               colour="grey") +
  scale_y_continuous(name="y-coordinate") + 
  scale_x_continuous(name="x-coordinate") + theme_bw()


ggplot(data_df_i, aes(xcoord,ycoord)) + 
  geom_point(aes(colour=as.factor(heattype)), shape=19, size=1.5) + 
  geom_tile(aes(fill=as.factor(heattype)), show.legend = FALSE)+ 
  scale_color_manual(values=c("gray","skyblue", "yellow", "darkorange", "red4"), 
                     name="", breaks=c(0,1,2,3,4), labels=c("insignificant","low-low","low-high","high-low","high-high")) +
  geom_segment(aes(x = xcoord1, y = ycoord1, xend = xcoord2, yend = ycoord2), 
               data = edges, 
               size = 0.5, 
               colour = "grey") +
  scale_y_continuous(name="y-coordinate") + 
  scale_x_continuous(name="x-coordinate") + theme_bw()  #+ guides(colour = FALSE)


ggplot(data_df_g, aes(xcoord,ycoord)) +  
  geom_point(alpha = 0) + 
  geom_tile()+ 
  geom_text(aes(label=heattype),hjust=0, vjust=0, size=3, check_overlap = T) +
  scale_colour_grey(guide='none') + 
  geom_segment(aes(x = xcoord1, y = ycoord1, xend = xcoord2, yend = ycoord2), 
               data = edges, 
               size = 0.5, 
               colour = "grey") +
  scale_y_continuous(name="y-coordinate") + 
  scale_x_continuous(name="x-coordinate") + theme_bw() 
  
  
ggplot_net(intnet, heattype='locmoran')


#----------------------------------COMPARE DATA WITH A REFERENCE-------------------------------------
load("../../Data/meanintense.rdata") # Reference data (gives 'net' list object)
intnet_und <- intensitynet(castellon, nodes, crim)
intnet_und <- RelateEventsToNetwork(intnet_und)

ref_g <- net$graph
g <- intnet_und$graph

ref_nodemean <- net$mean.intensity
ref_edgeint <- net$edge.intensity

nodemean <- vertex_attr(g, 'intensity')
edgeint <- edge_attr(g, 'intensity')


ref_neiV318 <- neighbors(ref_g, 'V318') # V318 -> V295 V300 V344
neiV318 <- neighbors(g, 'V318') # V318 -> V295 V300 V344

ref_nodemeanV318 <- ref_nodemean[[318]]
ref_edgeintV318 <- ref_edgeint[[318]]

nodemeanV318 <- nodemean[[318]]
edgeintV318 <- edge_attr(g, 'intensity', E(g)[get.edge.ids(g, c(318,295, 318,300, 318,344))])

ref_edgeint_318_344 <- ref_edgeint[[318]][3]
edgeint_318_300 <- edge_attr(g, 'intensity', E(g)[get.edge.ids(g, c(318,300))])
tt <- edge_attr(g, 'intensity', E(g, c(318,300))) # <-- correct, incorrect --> edge_attr(g, 'intensity', E(g)[318,300])


formatted_edgeint <- list()
for (node_id in V(g)) {
  formatted_edgeint[[node_id]] <- edge_attr(g, 'intensity', E(g)[get.edge.ids(g, c(318,295, 318,300, 318,344))])
}

# Plotting subgraph
PlotNeighborhood(intnet_und, 'V300')
PlotNeighborhood(intnet_und, 'V39')

v39v61_data <- list(p1=c(vertex_attr(g, 'xcoord', 'V39'), vertex_attr(g, 'ycoord', 'V39')), 
                    p2=c(vertex_attr(g, 'xcoord', 'V61'), vertex_attr(g, 'ycoord', 'V61')), 
                    ep=c(750422, 4430854))
class(v39v61_data) <- c(class(v39v61_data), 'netTools')
PointToSegment(v39v61_data)

ref_nodemean[[39]]
nodemean[[39]]

ref_edgeint[[39]]
edge_attr(g, 'intensity', E(g)[get.edge.ids(g, c(39,61))])

plot(nodes)
points(crimes, col='red')

coordsTest <- read.table("../../Data/coordsTEST.txt", header=TRUE) # says first column are rownames
plot(nodes)
points(test_df, col="red")
points(coordsTest, col="blue")

node_coords <- as.matrix(nodes)
rownames(node_coords) <- sprintf("V%s",seq(1:nrow(node_coords)))
data_df <- data.frame(xcoord = round(nodes$cx), 
                      ycoord = nodes$cy)
edgelist <- get.edgelist(g)
#convert to a four column edge data frame with source and destination coordinates
edges <- data.frame(node_coords[edgelist[,1],], node_coords[edgelist[,2],])
colnames(edges) <- c("xcoord1","ycoord1","xcoord2","ycoord2")




# New data
load("../../Data/network-CS.RData")

test_df <- data.frame(xcoord = round(c1213CS$data[,1]$x) - 14, ycoord = round(c1213CS$data[,2]$y) - 77)
test_df2 <- data.frame(xcoord = round(c1213CS$data[,1]$x), ycoord = round(c1213CS$data[,2]$y))

write.table(test_df, file="../../Data/crimes_corrected.txt") # keeps the rownames

t1 <- test_df[test_df$ycoord <= 4428800,] # Node coordinate reference point
t2 <- data_df[data_df$ycoord <= 4428800,] # Event coordinate reference point

test_mtx <- cbind(c1213CS[["domain"]][["lines"]][["marks"]][["aa_utmx0"]])
test_mtx <- cbind(test_mtx, c1213CS[["domain"]][["lines"]][["marks"]][["aa_utmy0"]])
# Test_mtx == test_df
test_mtx2 <- cbind(c1213CS[["domain"]][["lines"]][["marks"]][["aa_utmx1"]])
test_mtx2 <- cbind(test_mtx2, c1213CS[["domain"]][["lines"]][["marks"]][["aa_utmy1"]])

head(c1213CS$data)

data_df$cat <- "df1"
test_df$cat <- "df2"
test_df2$cat <- "df2"
all_df <- rbind(data_df, test_df)
ggplot(all_df, aes(x=xcoord,y=ycoord, color=cat)) + 
  geom_point(shape=19, size=1.5) +
  geom_segment(aes(x=xcoord1, y=ycoord1, xend = xcoord2, yend = ycoord2), 
               data=edges, 
               size = 0.5, 
               colour="grey") +
  scale_y_continuous(name="y-coordinate") + 
  scale_x_continuous(name="x-coordinate") + theme_bw()


#-----------------------------------SPATSTAT CHICAGO DATA---------------------------------------------
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


node_coords_obj <- list(node_coords = chicago_node_coords)
class(node_coords_obj) <- "netTools"
dist_mtx <- CalculateDistancesMtx(node_coords_obj)

#-------UNDIRECTED:

und_intnet_chicago <- intensitynet(chicago_adj_mtx, 
                               node_coords = chicago_node_coords, 
                               event_data = chicago_df)

start_time <- Sys.time()
und_intnet_chicago <- RelateEventsToNetwork(und_intnet_chicago)
end_time <- Sys.time()
print(end_time - start_time)

g <- und_intnet_chicago$graph

gen_cov <- NodeGeneralCorrelation(und_intnet_chicago, 
                                  dep_type = 'covariance', 
                                  lag_max = 2, 
                                  intensity = igraph::vertex_attr(g)$intensity)


gen_corr <- NodeGeneralCorrelation(und_intnet_chicago, 
                                   dep_type = 'correlation', 
                                   lag_max = 2, 
                                   intensity = igraph::vertex_attr(g)$intensity)

data_moran <- NodeLocalCorrelation(und_intnet_chicago, 
                                   dep_type = 'moran', 
                                   intensity = igraph::vertex_attr(und_intnet_chicago$graph)$intensity)
moran_i <- data_moran$correlation
intnet_mi <- data_moran$intnet

data_geary <- NodeLocalCorrelation(und_intnet_chicago, 
                                   dep_type = 'geary', 
                                   intensity = igraph::vertex_attr(und_intnet_chicago$graph)$intensity)
geary <- data_geary$correlation
intnet_gy <- data_geary$intnet


data_getis <- NodeLocalCorrelation(und_intnet_chicago, 
                                   dep_type = 'getis', 
                                   intensity = igraph::vertex_attr(und_intnet_chicago$graph)$intensity)
getis <- data_getis$correlation
intnet_gs <- data_getis$intnet

chicago_g <- und_intnet_chicago$graph

short_path_und <- ShortestPath(und_intnet_chicago, node_id1 = 'V1', node_id2 = 'V300', weight = 'intensity')





#-------DIRECTED:
dir_intnet_chicago <- intensitynet(chicago_adj_mtx, 
                                   node_coords = chicago_node_coords, 
                                   event_data = chicago_df,
                                   graph_type='directed')
dir_intnet_chicago <- RelateEventsToNetwork(dir_intnet_chicago)

gen_cov <- NodeGeneralCorrelation(dir_intnet_chicago, 
                                  dep_type = 'covariance', 
                                  lag_max = 2, 
                                  intensity = igraph::vertex_attr(dir_intnet_chicago$graph)$intensity_in)


data_moran <- NodeLocalCorrelation(dir_intnet_chicago, 
                                   dep_type = 'moran_i', 
                                   intensity = igraph::vertex_attr(dir_intnet_chicago$graph)$intensity_in)
moran_i <- data_moran$correlation
intnet <- data_moran$intnet

data_geary <- NodeLocalCorrelation(dir_intnet_chicago, 
                                   dep_type = 'geary', 
                                   intensity = igraph::vertex_attr(dir_intnet_chicago$graph)$intensity_in)
geary <- data_geary$correlation
intnet <- data_geary$intnet

short_path_dir <- ShortestPath(dir_intnet_chicago, node_id1 = 'V1', node_id2 = 'V300', weight = 'intensity')

#-------MIXED:
mix_intnet_chicago <- intensitynet(chicago_adj_mtx, 
                                   node_coords = chicago_node_coords, 
                                   event_data = chicago_df,
                                   graph_type='mixed')
mix_intnet_chicago <- RelateEventsToNetwork(mix_intnet_chicago)

#--------------------------------------Plot Chicago events-------------------------------------------
PlotNeighborhood(und_intnet_chicago, node_id = 'V300')

PlotHeatmap(und_intnet_chicago, heattype='moran', show_events = TRUE)
PlotHeatmap(und_intnet_chicago, heattype='geary')
PlotHeatmap(und_intnet_chicago, heattype='intensity')

chicago_edge_coords <- data.frame(xcoord1 = chicago[["domain"]][["lines"]][["ends"]][["x0"]],
                                  ycoord1 = chicago[["domain"]][["lines"]][["ends"]][["y0"]],
                                  xcoord2 = chicago[["domain"]][["lines"]][["ends"]][["x1"]],
                                  ycoord2 = chicago[["domain"]][["lines"]][["ends"]][["y1"]]  )
chicago_node_coords$cat <- "nodes"
assault_coordinates$cat <- "events"

chicago_plot <- rbind(chicago_node_coords, assault_coordinates)
ggplot(chicago_plot, aes(x=xcoord,y=ycoord, color=cat)) + 
  geom_point(shape=19, size=1.5) +
  geom_segment(aes(x=xcoord1, y=ycoord1, xend = xcoord2, yend = ycoord2), 
               data=chicago_edge_coords, 
               size = 0.5, 
               colour="grey") +
  scale_y_continuous(name="y-coordinate") + 
  scale_x_continuous(name="x-coordinate") + theme_bw()

sub_intnet_chicago <- ApplyWindow(dir_intnet_chicago, 
                                  x_coords = c(300, 900), 
                                  y_coords = c(500, 1000))
PlotHeatmap(sub_intnet_chicago, heattype='geary')


#-----------------------------------------TESTING WINDOW SUB_GRAPH--------------------------------------

sub_intnet <- ApplyWindow(intnet_und, x_coords = c(752500, 754500), y_coords = c(4429500, 4431500))
