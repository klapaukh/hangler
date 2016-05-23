#' Fits a spline to the set of points provided. 
#' 
#' @param points The points which make up the digitised curve
#' @return Coeffiecients for the spline 
#' @export
fitSpline <- NULL

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


#' Compute the value of the tangent angle spline at s.
#'
#' @param s The arc distance along the spline to compute the tangent at
#' @param si The arc distance at point measured point i. 
#'            si < s and si > sk for all other sk < s.
#' @param  sj The arc distance at the point after i.
#'            sj > s and sj < sh for all other sh > s.
#' @param thetai The tangent to the curve at distance si.
#' @param thetaj The tangent to the curve at distance sj.
#' @param deli The uncircleness of the spline between si and sj
#' @return The tangent to the spline at distance s along the spline
#' @export
computeSpline <- function(s, si, sj, thetai, thetaj, deli){
  t = (s - si) / (sj - si)
  previousContrib = thetai * (1 - t) 
  nextContrib     = thetaj * t
  errorContrib    = 0.5*deli*t*(1-t)

  return(previousContrib + nextContrib + errorContrib) 
}

#' Numerically estimate the value of deli for a curve segment of 
#' the spline.
#'
estimateDeli <- NULL


#' Secant method root finding
#' 
#' @param f The function to optimise to zero. Must take exactly 1 parameter
#' @param x1 First guess
#' @param x2 Second guess
#' @param maxIter Maximum number of iterations to run for
#' @param targetError Maximum error before the solution will be accepted as correct
#' @return X value which returns zero. If root finding fails, returns NA
#' @export
secantMethod <- function(f, x1, x2, maxIter, targetError){
  x = Reduce(function(guesses, iter){
      fx1 = f(guesses[1])
      if(abs(fx1) < targetError) return(guesses)
      fx2 = f(guesses[2])
      xNext = (guesses[2]*fx1 - guesses[1]*fx2) / (fx1 - fx2)
      return(c(xNext,guesses[1]))
    }, 1:maxIter, c(x1,x2));
  x = x[1]
  return(ifelse(abs(f(x)) <= targetError, x, NA))  
}
