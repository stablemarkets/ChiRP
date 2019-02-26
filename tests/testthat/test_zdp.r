context("ZDP Tests")

test_that("Test Iteration Errors", {
  set.seed(1)
  n<-200 ## generate from clustered, skewed, data distribution
  d <- data.frame(y = rnorm(n), x = rnorm(n))
  d$y[1:30] <- 0

  expect_error( ZDPMix(d_train = d, formula = y ~ x, burnin=200, iter=100) )
  expect_error( ZDPMix(d_train = d, formula = y ~ x, burnin=200, iter=200) )
  expect_error( ZDPMix(d_train = d, formula = y ~ x, burnin=-200, iter=-200) )
  expect_error( ZDPMix(d_train = d, formula = y ~ x, burnin=1.5, iter=4.5) )
  
  d <- data.frame(y = rnorm(n), x = rnorm(n))
  expect_error( ZDPMix(d_train = d, formula = y ~ x, burnin=100, iter=200) )
})

test_that("Test Init K Errors", {
  set.seed(1)
  n<-200 ## generate from clustered, skewed, data distribution
  d <- data.frame(y = rnorm(n), x = rnorm(n))
  d$y[1:30] <- 0
  
  expect_error( ZDPMix(d_train = d, formula = y ~ x, init_k = -10) )
  expect_error( ZDPMix(d_train = d, formula = y ~ x, init_k = 1.5) )
  expect_error( ZDPMix(d_train = d, formula = y ~ x, init_k = 300) )
  expect_error( ZDPMix(d_train = d, formula = y ~ x, init_k = 'test') )
  expect_error( ZDPMix(d_train = d, formula = y ~ x, init_k = c(1,2)) )
})


test_that("Test d_train inputs", {
  set.seed(1)
  n<-200 ## generate from clustered, skewed, data distribution
  d <- data.frame(y = rnorm(n), x = rnorm(n), z=sample(c(1,2), size = n, replace = T))
  d$y[1:30] <- 0
  
  expect_error( ZDPMix(d_train = d, formula = y ~ x + z, init_k = 10) )
  expect_error( ZDPMix(d_train = d, formula = y ~ x + z + A, init_k = 10) )
  expect_error( ZDPMix(d_train = data.frame(), formula = y ~ x + z + A, init_k = 10) )
})

test_that("Test d_test inputs", {
  set.seed(1)
  n<-50 ## generate from clustered, skewed, data distribution
  d <- data.frame(y = rnorm(n), x = rnorm(n), z=sample(c(1,0), size = n, replace = T))
  d$y[1:10] <- 0
  
  n<-50 ## generate from clustered, skewed, data distribution
  dt <- data.frame(y = rnorm(n), x = rnorm(n), z=sample(c(1,0), size = n, replace = T))
  dt$y[1:10] <- 0
  
  expect_error( ZDPMix(d_train = d, d_test = dt[,-2], formula = y ~ x + z, init_k = 10) )
  expect_error( ZDPMix(d_train = d, d_test = data.frame(), formula = y ~ x + z + A, init_k = 10) )
})

test_that("Test ZDP outputs", {
  set.seed(1)
  n<-50 ## generate from clustered, skewed, data distribution
  d <- data.frame(y = rnorm(n), x = rnorm(n), z=sample(c(1,0), size = n, replace = T))
  d$y[1:10] <- 0
  
  n<-30 ## generate from clustered, skewed, data distribution
  dt <- data.frame(y = rnorm(n), x = rnorm(n), z=sample(c(1,0), size = n, replace = T))
  dt$y[1:10] <- 0
  
  res <- ZDPMix(d_train = d, d_test = dt, formula = y ~ x + z, init_k = 10)
  
  expect_length( res, 2 )
  expect_length( res$cluster_inds, 2 )
  expect_length( res$predictions, 2 )
  
  expect_equal( dim(res$predictions$train), c(50, 900) )
  expect_equal( dim(res$predictions$test),  c(30, 900) )
  
  expect_equal( dim(res$cluster_inds$train), c(50, 900) )
  expect_equal( dim(res$cluster_inds$test),  c(30, 900) )
})

