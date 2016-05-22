#' Fits a spline to the set of points provided. 
#' 
#' @param points The points which make up the digitised curve
#' @return 
fitspline <- NULL

#' Finds the center of a circle given three points
#' Equations taken from http://www.ambrsoft.com/TrigoCalc/Circle3D.htm
#' 
#' @param points that make up the three points on the circumference
#' @return a vector [x,y] that is the center of the circle 
findCenter <- function(x1, y1, x2, y2, x3, y3){
  A = x1 * (y2 - y3) - 
      y1 * (x2 - x3) + 
      x2 * y3 - 
      x3 * y2

  B = (x1 * x1 + y1 * y1)*(y3 - y2) + 
      (x2 * x2 + y2 * y2)*(y1 - y3) +
      (x3 * x3 + y3 * y3)*(y2 - y1) 

  C = (x1 * x1 + y1 * y1)*(x2 - x3) + 
      (x2 * x2 + y2 * y2)*(x3 - x1) +
      (x3 * x3 + y3 * y3)*(x1 - x2)

  x = (-B) / (2 * A)
  y = (-C) / (2 * A)

  return(c(x,y))
}

#' Compute Euclidean distance between two points
#'
#' @param two points to find the distance between
#' @return Euclidean distance between the two points
distance <- function(x1, y1, x2, y2){
  dx = x1 - x2
  dy = y1 - y2
  sqDist = dx * dx + dy * dy
  return(sqrt(sqDist))
}

#' Find that tangent angle to a point on a circle given it's center and 
#' a point on its circumference
#' 
#' @param xc, yc coordinates of the circle's center
#' @param x,y coordinates of a point on the cicle's circumference
#' @return tangent angle to the circle
findTangentAngle <- function(xc, yc, x, y){
  x = x - xc 
  y = y - yc
  theta = atan2(y, x) - (pi / 2);
  theta = ifelse(theta < 0, 2*pi + theta, theta) 
  # Tangent is 90 degrees of the line itself
  return(theta) 
}
