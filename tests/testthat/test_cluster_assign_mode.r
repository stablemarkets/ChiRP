context("Cluster Assign Mode Tests")

test_that("Test Cluster Assignments", {
  set.seed(1)
  n<-200 ## generate from clustered, skewed, data distribution
  d <- data.frame(y = rnorm(n), x = rnorm(n))
  d$y[1:30] <- 0

  res <- ZDPMix(d_train = d, formula = y ~ x, burnin=10, iter=50)
  
  train_clus <- cluster_assign_mode(res$cluster_inds$train)
  
  expect_equal( dim(train_clus$adjmat), c(n, n)  )
  expect_length( train_clus$class_mem,  n)
})

