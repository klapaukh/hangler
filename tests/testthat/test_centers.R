library(hangler)

context("Circle centre detection")

test_that("findCenter returns 2 values", {
  x = sin(runif(min = 0, max = 1, n = 9))
  y = cos(runif(min = 0, max = 1, n = 9))
  expect_equal(length(findCenter(x[1], y[1], x[2], y[2], x[3], y[3])), 2)
  expect_equal(length(findCenter(x[4], y[4], x[5], y[5], x[6], y[6])), 2)
  expect_equal(length(findCenter(x[7], y[7], x[8], y[8], x[9], y[9])), 2)
 })


test_that("Distances make sense", {
  expect_equal(distance(0,1,0,0), 1)
  expect_equal(distance(1,0,0,0), 1)

  theta = runif(min = 0, max = 2*pi, n = 9)
  x = sin(theta)
  y = cos(theta)

  expect_equal(distance(x[1],y[1],0,0), 1)
  expect_equal(distance(x[2],y[2],0,0), 1)
  expect_equal(distance(x[3],y[3],0,0), 1)
  expect_equal(distance(x[4],y[4],0,0), 1)
  expect_equal(distance(x[5],y[5],0,0), 1)
  expect_equal(distance(x[6],y[6],0,0), 1)
  expect_equal(distance(x[7],y[7],0,0), 1)
  expect_equal(distance(x[8],y[8],0,0), 1)
  expect_equal(distance(x[9],y[9],0,0), 1)
 })

test_that("findCenter radii are correct", {
  theta = runif(min = 0, max = 2*pi, n = 9)
  x = sin(theta)
  y = cos(theta)

  p1 = findCenter(x[1], y[1], x[2], y[2], x[3], y[3])
  p2 = findCenter(x[4], y[4], x[5], y[5], x[6], y[6])
  p3 = findCenter(x[7], y[7], x[8], y[8], x[9], y[9])

  expect_equal(p1[1], 0)
  expect_equal(p1[2], 0)

  expect_equal(p2[1], 0)
  expect_equal(p2[2], 0)

  expect_equal(p3[1], 0)
  expect_equal(p3[2], 0)
 })

test_that("findCenter radii are the same", {
  theta = runif(min = 0, max = 2*pi, n = 9)
  x = sin(theta)
  y = cos(theta)

  p1 = findCenter(x[1], y[1], x[2], y[2], x[3], y[3])
  p2 = findCenter(x[4], y[4], x[5], y[5], x[6], y[6])
  p3 = findCenter(x[7], y[7], x[8], y[8], x[9], y[9])

  expect_equal(distance(p1[1], p1[2], x[1], y[1]), distance(p1[1], p1[2], x[2], y[2]))
  expect_equal(distance(p1[1], p1[2], x[1], y[1]), distance(p1[1], p1[2], x[3], y[3]))
                                                                                
  expect_equal(distance(p2[1], p2[2], x[4], y[4]), distance(p2[1], p2[2], x[5], y[5]))
  expect_equal(distance(p2[1], p2[2], x[4], y[4]), distance(p2[1], p2[2], x[6], y[6]))
                                                                                
  expect_equal(distance(p3[1], p3[2], x[7], y[7]), distance(p3[1], p3[2], x[8], y[8]))
  expect_equal(distance(p3[1], p3[2], x[7], y[7]), distance(p3[1], p3[2], x[9], y[9]))
 })

test_that("findTangentAngle solutions make sense", {
  theta = seq(0,2*pi, 0.01)
  tangent = theta - (pi / 2)
  tangent = ifelse(tangent < 0, 2*pi + tangent, tangent)

  xc = 0
  yc = 0 
  x = cos(theta)
  y = sin(theta)

  sapply(1:length(theta), function(idx){
    expect_equal(findTangentAngle(xc,yc,x[idx],y[idx]), tangent[idx])
  })
 })


test_that('secant Method search sometimes works', {
  f1 <- sin
  f2 <- tan
  f3 <- function(x)  x**3 + x*x + 5*x + 7

  sapply(c(f1,f2,f3), function(f) { 
    expect_less_than(abs(f(secantMethod(f, 0.4, 0.2, 500, 1e-4))), 1e-4)
  });

 })

test_that('false position Method search sometimes works', {
  f1 <- sin
  f2 <- tan
  f3 <- function(x)  x**3 + x*x + 5*x + 7

  sapply(c(f1,f2,f3), function(f) { 
    expect_less_than(abs(f(falsePositionMethod(f, 0.4, 0.2, 500, 1e-4))), 1e-4)
  });

 })

test_that('simpsons rule integration works', {
  f1 = function(x) 9 
  f2 = sin
  f3 = function(x) x + 2

  expect_equal(simpsonsRule(f1, 0, 5, 1000), 45)
  expect_equal(simpsonsRule(f2, 0, pi, 1000), 2)
  expect_equal(simpsonsRule(f3, 0, 5, 1000), 22.5)
  
 })
