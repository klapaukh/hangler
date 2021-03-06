#' Turn a shape into a set of tangents as offsets from the unit
#' circle.
#' 
#' @param x The x coordinates for the points which make up the digitised curve
#' @param y The y coordinates for the points which make up the digitised curve
#' @return Tangents 
#' @export
computeTangents <- function(x,y){
  if(length(x) != length(y)){
    stop(paste0("x (", length(x), " elements) and y (",
                length(y), " elements) must be the same length but are not"))
  }

  nPoints = length(y)

  # Start by finding the center of every triple of points as
  # a circle
  centers = lapply(1:nPoints, function(i){
    ip = i - 1
    inext = i + 1

    ip = ifelse(ip < 1, nPoints - ip, ip) 
    inext = ifelse(inext > nPoints, inext - nPoints, inext) 
    findCenter(x[ip],y[ip],x[i],y[i],x[inext],y[inext])
  })

  #extract the x and y coordinates
  xc = sapply(centers, function(c) c[1])
  yc = sapply(centers, function(c) c[2])

  #Compute the tangents
  tangents = sapply(1:nPoints, function(i) findTangentAngle(xc[i], yc[i], x[i], y[i]))

  # Smooth the values to avoid discontinuities from trig functions 
  tangents = sequenceAngles(tangents)

  # Reset the curve to start from zero, by turning it into a sequence of 
  # step changes and recombining
  diffs = diff(tangents)

  #diffs = sapply(diffs, shortestAngle)

  samples = cumsum(c(0,diffs))

  return(samples)
}

#' Resample the tangents curve to a given number of samples.
#' Uses the Hangle spline function. This function has a lot of 
#' numerical approximation in it. You should always compare its
#' output to the original to make sure nothing when wrong.
#'
#' @param x The x coordinates for the points which make up the digitised curve
#' @param y The y coordinates for the points which make up the digitised curve
#' @param tangents A vector of tangents
#' @param length.out The length of the output array
#' @return A resampled tangent vector with length length.out
resampleTangents <- function(x,y,tangents,length.out=1024, 
                             integral.iter=5, solver.error= 1e-10,
                             solver.max.iter=1000){
  nPoints = length(tangents)

  # Try approximate the deli values for the spline
  deli = sapply(1:nPoints, function(i) {
    inext = i + 1
    inext = ifelse(inext > nPoints, inext - nPoints, inext) 
  
    tryCatch(solveDeli(x[inext] - x[i], y[inext] - y[i], tangents[i], tangents[inext],
                       maxIter=solver.max.iter, targetError = solver.error, 
                       integralIter=integral.iter), 
      error = function(e){ 
                stop(paste(e,"\nFailed to compute deli for spline between",i, 
                           "and",inext))
              }
    )
  })

  # Using the deli values try approximate the ds values
  ds = sapply(1:nPoints, function(i) {
    inext = i + 1
    inext = ifelse(inext > nPoints, inext - nPoints, inext) 
 
    tryCatch(solveDs(x[inext] - x[i], deli[i], tangents[i], tangents[inext],
                     integralIter = integral.iter), 
      error = function(e){ 
                stop(paste(e, "\nFailed to compute ds for spline between",i,
                          "and", inext))
              }
    )
  })

  # Need them to be positive
  ds = abs(ds)

  # si are the cumulative sum of ds value
  si = cumsum(c(0,ds))

  # Recompute the right number of tangents from the spline
  samples = sapply(seq(0,si[length(si)],length = length.out), function(s){
    i = which(si < s)
    if(length(i)==0){
      i = 1
    }else{
      i = max(i)
    }

    inext = i + 1
    sinext = ifelse(inext > length(si), inext-length(si), inext)
    tinext = ifelse(inext > length(tangents), inext-length(tangents), inext)

    computeSpline(s, si[i], si[sinext], tangents[i], tangents[tinext], deli[i])
  })

  return(samples)
}

#' Make a sequence of angles in radians all follow each other in with a value
#' as close to the previous as possible
#'
#' @param angles An ordered vector of angles (in radians)
#' @return a vector of angles where the adjacent angles are as close to each
#'  other in magnitude as possible
#' @export
sequenceAngles <- function(angles){
  Reduce(function(soFar, angle){
    if(length(soFar) == 0) { 
      return(angle)
    }
    lastAngle = soFar[length(soFar)]
    while(angle < lastAngle){
      angle = angle + pi
    }
    while(angle - pi > lastAngle){
      angle = angle - pi
    }
    nangle = angle - pi
    nangle = ifelse(abs(nangle-lastAngle) < abs(angle-lastAngle), nangle,angle)
    return(c(soFar, nangle))
  }, angles, c())
}
  
#' Convert an angle in radians into a standard form.
#' Takes any angle in radians and returns an angle beween pi and -pi.
#'
#' @param Theta an angle in radians
#' @return an angle in radians between pi and -pi
boundAngle <- function(theta){
  if(theta >= 2*pi){
    theta = theta %% (2*pi)
  } else if(theta <= -2*pi){
    theta = theta  %% (-2*pi)
  }
  if(theta > pi){
    return(theta - 2*pi)
  }else if(theta < -pi){
    return(2*pi + theta)
  }else{
    return(theta)
  }
}

#' Finds the center of a circle given three points
#' Equations taken from http://www.ambrsoft.com/TrigoCalc/Circle3D.htm
#' 
#' Sometimes this doesn't work because the circle has infinite size. 
#' In that case make a up center that will give the correct angle.
#'
#' @param points that make up the three points on the circumference
#' @return a vector [x,y] that is the center of the circle 
#' @export
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

  if(is.finite(x) & is.finite(y)){
    return(c(x,y))
  }

  #make a unit circle with the straight line tangent
  tangentX = x3 - x1
  tangentY = y3 - y1

  #Now rotate it by 90 degrees
  temp = tangentX
  tangentX = tangentY
  tangentY = -temp

  #Add it back to the middle one to get the circle center
  x = x2 + tangentX
  y = y2 + tangentY
  
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
#' @export
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
  computeSplineT(t, thetai, thetaj, deli)
}

#' Compute the value of the tangent angle spline at s.
#'
#' @param t The arc distance along the spline to compute the tangent at
#' @param thetai The tangent to the curve at distance si.
#' @param thetaj The tangent to the curve at distance sj.
#' @param deli The uncircleness of the spline between si and sj
#' @return The tangent to the spline at distance s along the spline
#' @export
computeSplineT <- function(t, thetai, thetaj, deli){
  if(abs(thetai - thetaj) > pi){
    #E.g. 6.2 & 0.2 so the linear interpolation will be going the long way round
    if(thetai < thetaj){
      thetai = thetai + 2*pi
    } else {
      thetaj = thetaj + 2*pi
    }
  }

  previousContrib = thetai * (1 - t) 
  nextContrib     = thetaj * t
  errorContrib    = 0.5*deli*t*(1-t)

  return(previousContrib + nextContrib + errorContrib) 
}

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
      if(length(guesses) == 1) return(guesses)
      fx1 = f(guesses[1])
      if(abs(fx1) < targetError) return(guesses[1])
      fx2 = f(guesses[2])
      xNext = (guesses[2]*fx1 - guesses[1]*fx2) / (fx1 - fx2)
      return(c(xNext,guesses[1]))
    }, 1:maxIter, c(x1,x2));
  x = x[1]
  return(ifelse(abs(f(x)) <= targetError, x, NA))  
}

#' False position method root finding
#' 
#' @param f The function to optimise to zero. Must take exactly 1 parameter
#' @param x1 First guess
#' @param x2 Second guess
#' @param maxIter Maximum number of iterations to run for
#' @param targetError Maximum error before the solution will be accepted as correct
#' @return X value which returns zero. If root finding fails, returns NA
#' @export
falsePositionMethod <- function(f, x1, x2, maxIter, targetError){
  x = Reduce(function(guesses, iter){
      if(length(guesses) == 1) return(guesses)
      fx1 = f(guesses[1])
      if(abs(fx1) < targetError) return(guesses[1])
      fx2 = f(guesses[2])
      xNext = (guesses[1]*fx2 - guesses[2]*fx1) / (fx2 - fx1)
      if(is.nan(xNext)){
           stop(paste("iter:",iter,"- xNext is NaN, f(x1): f(",guesses[1],") =",fx1,
                      "f(x2): f(",guesses[2]),") =",fx2)
      }
      if(xNext == guesses[2]) {
        print(paste("Adjusting on iter:",iter))
        xNext = xNext*0.99  
      }
      return(c(guesses[2],xNext))
    }, 1:maxIter, c(x1,x2));
  x = x[length(x)]
  return(ifelse(abs(f(x)) <= targetError, x, NA))  
}


#' Newton-Raphson root finding
#'
#' This version of Newton-Raphson has been modified to not move more than
#' 1 unit at a time (because the shallow bits at the top of the curve are 
#' liable to vastly estimate distance), and to not move to a place that is
#' a worse solution than the current (instead it bisects back until it 
#' finds such a place.
#'
newtonRaphsonMethod <- function(f, guess=0, maxIter=1000, targetError=1e-6) {
  derivativeDistance = 1e-6
  x = Reduce(function(guess, iter){
    x0  = guess[1]
    fx0 = guess[2] # Note that fx0 is f(x0)
    if(abs(fx0) < targetError) { 
      return(c(x0,fx0))
    }
    fx0e = f(x0 + derivativeDistance) # f(x0 + epsilon)
    dfx0 = (fx0e - fx0) / derivativeDistance # f'(x0)
    shift = fx0 / dfx0
    if(abs(shift) > 1){
      shift = sign(shift)
    }  else if(shift == 0) {
      xNew = x0 - sign(x0) * derivativeDistance
      fxNew = f(xNew)
      return (c(xNew, fxNew))
    }
    xNew = x0 - shift
    fNew = f(xNew)
    while(abs(fNew) > abs(fx0)) {
      shift = shift / 2
      xNew = x0 - shift
      fNew = f(xNew)
    }
    return (c(xNew,fNew))
  } , 1:maxIter, c(guess, f(guess)));
  return(ifelse(abs(x[2]) <= targetError, x[1], NA)) 
}

#' Numerical integration using Simpson's rule
#'
#' @export
simpsonsRule <- function(f, a, b, nBins){
  cells = seq(a ,b, length.out = nBins + 1)
  Reduce(function(soFar, nextCell){
    current = simpsonsRuleCell(f, soFar[2], nextCell)
    return(c(soFar[1] + current, nextCell))
   }, cells[-1], c(0, a))[1]
}

#' Numerical integration of a single cell using
#' Simpson's rule. 
#'
#' @param f The function to integrate. Must take exactly 1 parameter
#' @param a Lower region bound
#' @param b Upper region bound
#' @return Numerical approximation of the integral of f on [a,b]
simpsonsRuleCell <- function(f, a, b){
  (b - a) * (f(a) + 4* f((a+b)/2) + f(b)) / 6
}

#' Find the delta i values along the spiline
#'
#' @export
solveDeli <- function(dx, dy, thetai, thetaj, maxIter=1000, targetError = 1e-10, integralIter=5) {
  deli = newtonRaphsonMethod(function(deli){

     cosInt = simpsonsRule(function(x) {cos(computeSplineT(x, thetai, thetaj, deli))},0,1,integralIter)
     sinInt = simpsonsRule(function(x) {sin(computeSplineT(x, thetai, thetaj, deli))},0,1,integralIter)
     if(!is.finite(cosInt)){
       stop(paste("cosInt is not finite for deli:",deli))
     }
     if(!is.finite(sinInt)){
       stop(paste("sinInt is not finite for deli: ",deli))
     }
     return((dy * cosInt) - (dx * sinInt))
   }, 0, maxIter = maxIter, targetError = targetError )
  return (deli)
}

#' Find the delta s values along the spline
#'
#' @export
solveDs <- function(dx, deli, thetai, thetaj, integralIter=5){
  dx / simpsonsRule(function(x) { cos(computeSplineT(x, thetai, thetaj, deli)) },0,1,integralIter)
}

#' Filter a digitised curve to reduce noise
#'
#' @param points Points to filter (e.g. the x coordinates)
#' @param outputLength The number of resampled tangent points that will be fed into the FFT
#' @param iterations The number of smoothing iterations to do
#' @export
digitisationFilter <- function(points, 
                               outputLength, 
                               iterations = ceiling(2*((length(points) / outputLength)**2))){
  Reduce(function(input, iteration) {
    n          = length(points)
    rightShift = c(input[n],input[-n])
    leftShift  = c(input[-1],input[1])
    newPoints  = 0.25*rightShift + 0.5* input + 0.25 * leftShift
    return(newPoints)
  }, 1:iterations, points)
}

#' Remove the circle effect from the tangents
#'
#' Removes the effect of the circle on the tangents. This basically 
#' Turns the tangents into the deviations from the circle, which 
#' can be run through an FFT
#' @param tangents The tangents to flatting (These should run from 0 to 2*pi)
#' @return The tangents with the circle component removed
flattenTangents <- function(tangents) {
  circle = seq(0, 2*pi, length.out=length(tangents))
  return(tangents - circle)
}

# vim: expandtab sw=2 ts=2 tw=80  
