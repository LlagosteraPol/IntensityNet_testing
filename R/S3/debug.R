# PointToSegment <- function(obj){
#   p1 <- obj$p1
#   p2 <- obj$p2
#   ep <- obj$ep
#   
#   # Return minimum distance between line segment p1p2 and point ep
#   l2 <- sqrt(sum((p1 - p2) ^ 2))
#   
#   if (l2 == 0.0) return( sqrt(sum((ep - p2) ^ 2)) )  # p1 == p2 case
#   # Consider the line extending the segment, parameterized as p1 + t (p2 - p1).
#   # We find projection of point p onto the line. 
#   # It falls where t = [(ep-p1) . (p2-p1)] / |p2-p1|^2
#   # We clamp t from [0,1] to handle points outside the segment p1p2.
#   t <- ((ep[1] - p1[1])*(p2[1] - p1[1]) + (ep[2] - p1[2])*(p2[2] - p1[2])) / l2
#   t <- max(0, min(1, t))
#   
#   projection2 <- c( p1[1] + t * (p2[1] - p1[1]), p1[2] + t * (p2[2] - p1[2]))
#   
#   projection = p1 + t * (p2 - p1)  # Projection falls on the segment
#   return( sqrt(sum((ep - projection) ^ 2)) )
# }

PointToSegment <- function(obj) {
  p1 <- obj$p1
  p2 <- obj$p2
  ep <- obj$ep
  A <- ep[1] - p1[1]
  B <- ep[2] - p1[2]
  C <- p2[1] - p1[1]
  D <- p2[2] - p1[2]
  
  dot <- A * C + B * D
  len_sq <- C * C + D * D
  param <- -1
  if (len_sq != 0){
    param <- dot / len_sq # in case of 0 length line
  } 
  
  
  if (param < 0) {
    xx <- p1[1]
    yy <- p1[2]
  }
  else if (param > 1) {
    xx <- p2[1]
    yy <- p2[2]
  }
  else {
    xx <- p1[1] + param * C
    yy <- p1[2] + param * D
  }
  
  dx <- ep[1] - xx
  dy <- ep[2] - yy
  return(sqrt(dx * dx + dy * dy))
}

p1 <- c(3,4)
p2 <- c(7,4)
ep <- c(2,4)
PointToSegment(list(p1=p1,p2=p2,ep=ep))