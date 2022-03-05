


ep <- matrix(c(639.1747, 1190.523), nrow=1)
p1 <- rbind(c(0.3894739,    1253.8027), c(109.6830137,	1251.7715), c(109.6830137,	1251.7715))
p2 <- rbind(c(109.6830137,	1251.7715), c(111.1897363,	1276.5601), c(197.9987626,	1251.1532))

res1 <- matrix()
for(i in 1:nrow(p1)){
  dist_obj1 <- list(p1 = p1[i,], p2 = p2[i,], ep = ep)
  class(dist_obj1) <- 'netTools'
  res1<- rbind(res1, PointToSegments(dist_obj1))
}

dist_obj2 <- list(p1 = p1, p2 = p2, ep = ep)
class(dist_obj2) <- 'netTools'
res2 <- PointToSegment(dist_obj2)