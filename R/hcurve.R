#' Convert a list of curve tangents to xy points for plotting.
#' Each step along the curve is assumed to be equally spaced. 
#'
#' @param tangents Equally spaced tangents along the curve
#' @return A list contaning a vector x points and a vector of y points.
#' @export
tangentsToXY <- function(tangents){

  points = Reduce(function(soFar, t){
    xc = cos(t)
    yc = sin(t)

    x = soFar[1] + xc
    y = soFar[2] + yc

    c(x,y)
    },tangents, c(0,0), accumulate=TRUE) 

  list(x= lapply(points, function(x) x[1]), 
       y= lapply(points, function(x) x[2]))

}


#' Convert fourier coefficients to Ak, sk values. 
#' 
#' @param Coeffs The output of fft
#' @return Ak and sk values for the reverse fourier
#' @export
fftoToCoeffs <- function(fftCoeffs){
  coeffs = lapply(1:length(fftCoeffs) , function(i){
    k = i-1
    ak = Re(fftCoeffs[i])
    bk = Im(fftCoeffs[i])

    Ak = 2 * sqrt( ak * ak + bk * bk)
    sk = -atan(bk/ak)/k

    return(c(Ak,sk))
  })

  return(list(Ak = sapply(coeffs, function(i) i[1]), 
              sk = sapply(coeffs, function(i) i[2])))
}


#' Get tangent at curve. Length is normalised to 2pi
#'
#' @export
getTangentAtS <- function(s, coeffs){
  maxK = length(coeffs$sk)
  s + sum(sapply(2:maxK, function(idx){
    k = idx - 1
    Ak = coeffs$Ak[idx]
    sk = coeffs$sk[idx]
    return(Ak*cos(k*(s-sk)))
  }))
}
