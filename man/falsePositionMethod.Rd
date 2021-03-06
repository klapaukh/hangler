% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spline.R
\name{falsePositionMethod}
\alias{falsePositionMethod}
\title{False position method root finding}
\usage{
falsePositionMethod(f, x1, x2, maxIter, targetError)
}
\arguments{
\item{f}{The function to optimise to zero. Must take exactly 1 parameter}

\item{x1}{First guess}

\item{x2}{Second guess}

\item{maxIter}{Maximum number of iterations to run for}

\item{targetError}{Maximum error before the solution will be accepted as correct}
}
\value{
X value which returns zero. If root finding fails, returns NA
}
\description{
False position method root finding
}
