context("credible gradient")

test_that("Test Input Errors", {
  set.seed(1)
  n<-200 ## generate from clustered, skewed, data distribution
  d <- data.frame(y = rnorm(n), x = rnorm(n))
  d$y[1:30] <- 0
  
  set.seed(1)
  
  N = 200
  x<-seq(1,10*pi, length.out = N) # confounder
  y<-rnorm(n = length(x), sin(.5*x), .07*x )
  d <- data.frame(x=x, y=y)
  d$x <- as.numeric(scale(d$x))
  d$y <- as.numeric(scale(d$y))
  
  d_test = data.frame(x=seq(max(d$x), max(d$x+1 ), .01 ))
  
  
  res = fDPMix(d_train = d, d_test = d_test, formula = y ~ x,
               iter=100, burnin=50, tau_x = c(.01, .001) )
  
  plot(d$x,d$y, pch=20, xlim=c(min(d$x), max(d$x)+1 ), col='gray')
  
  expect_error( credible_gradient(d$x, 'character') )
  expect_error( credible_gradient(d$x, res$train, col_gradient = 1) )
  expect_error( credible_gradient(d$x, res$train, col_mean_line = 1) )
  expect_error( credible_gradient(d$x[1:5], res$train) )
})

