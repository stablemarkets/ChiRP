class_update_train_ndp <- function(n, K , alpha, name_new, uniq_clabs, clabs,
                         y, x, x_cat_shell, x_num_shell,
                         cat_idx, num_idx,
                         beta_shell,  psi_shell,
                         beta_new, psi_new, cat_new, num_new){

  # shells
  c_shell <- matrix(NA, nrow = n, ncol = K + 1)
  colnames(c_shell) <- c(uniq_clabs, name_new)

  ## prior for existing cluster
  pr_exist <- matrix(NA, nrow = n, ncol = K)
  colnames(pr_exist) <- uniq_clabs

  pr_new <- numeric(length = n)

  clabs <- factor(clabs, levels = uniq_clabs)
  for(j in 1:n){
    n_min_j <- table(clabs[-j] )
    pr_exist[j,] <- log( (alpha*(n_min_j==0) + n_min_j)/(j + alpha - 1) )
    pr_new[j] <-  log( (alpha) / (j + alpha - 1) )
  }


  ## existing clusters
  for(k in uniq_clabs){

    lk_exist <- numeric(length = n)

    for(p in cat_idx){
      lk_exist <- lk_exist + dbinom(x[ ,p], 1, x_cat_shell[[p]][1,k], T)
    }

    for(p in num_idx){
      lk_exist <- lk_exist + dnorm(x[, p],
                                   x_num_shell[[p]][[1]][,k],
                                   sqrt(x_num_shell[[p]][[2]][,k]), T )
    }

    lk_exist <- lk_exist + dnorm(x = y,
                                 mean = x %*% beta_shell[,k, drop=F],
                                 sd = sqrt(psi_shell[,k]), T)

    c_shell[,k] <- lk_exist + pr_exist[,k]
  }


  ## New clusters
  lk_new <- numeric(length = n)

  for(p in cat_idx){
    lk_new <- lk_new + dbinom(x[ ,p], 1, cat_new[p], T)
  }

  for(p in num_idx){
    lk_new <- lk_new + dnorm(x[, p], num_new[1,p], sqrt(num_new[2,p]), T)
  }

  lk_new <- lk_new + dnorm(x = y,
                           mean = x %*% t(beta_new),
                           sd = sqrt(psi_new), T)

  c_shell[,name_new] <- lk_new + pr_new

  weights <- t(apply(c_shell, 1, function(x) exp(x)/sum(exp(x))  ))

  return(weights)
}

class_update_test_ndp <- function(n, K , alpha, name_new, uniq_clabs, clabs, x,
                             x_cat_shell, x_num_shell,
                             cat_idx, num_idx,
                             cat_new, num_new){

  # shells
  c_shell <- matrix(NA, nrow = n, ncol = K + 1)
  colnames(c_shell) <- c(uniq_clabs, name_new)

  ## prior for existing cluster
  clabs <- factor(clabs, levels = uniq_clabs)
  pr_exist <- log(table(clabs)/(length(clabs)+alpha))

  ## existing clusters
  for(k in uniq_clabs){

    lk_exist <- numeric(length = n)

    for(p in cat_idx){
      lk_exist <- lk_exist + dbinom(x[ ,p], 1, x_cat_shell[[p]][1,k], T)
    }

    for(p in num_idx){

      lk_exist <- lk_exist + dnorm(x[, p],
                                   x_num_shell[[p]][[1]][,k],
                                   sqrt(x_num_shell[[p]][[2]][,k]), T )
    }

    c_shell[,k] <- lk_exist + pr_exist[k]
  }


  ## New clusters
  lk_new <- numeric(length = n)

  for(p in cat_idx){
    lk_new <- lk_new + dbinom(x[ ,p], 1, cat_new[p], T)
  }

  for(p in num_idx){
    lk_new <- lk_new + dnorm(x[, p], num_new[1,p], sqrt(num_new[2,p]), T)
  }

  c_shell[,name_new] <- lk_new + log(alpha/(length(clabs)+alpha))

  weights <- t(apply(c_shell, 1, function(x) exp(x)/sum(exp(x))  ))

  return(weights)
}

post_pred_draw_train_ndp <- function(n, x, pc, beta_shell, psi_shell){

  y_pp <- numeric(length = n)

  clabs <- unique(pc)

  for( k in clabs ){
    y_pp[pc==k] <- rnorm(n = sum(pc==k),
                         mean = x[pc==k,, drop=F] %*% beta_shell[,k,drop=F],
                         sd = sqrt(psi_shell[1,k]) )
  }

  return(y_pp)
}


post_pred_draw_test_ndp <- function(n, x, pc, beta_shell, psi_shell,
                           name_new, beta_new, psi_new){

  y_pp <- numeric(length = n)

  clabs <- unique(pc)

  for( k in setdiff(clabs, name_new) ){
    y_pp[pc==k] <- rnorm(n = sum(pc==k),
                         mean = x[pc==k,, drop=F] %*% beta_shell[,k,drop=F],
                         sd = sqrt(psi_shell[1,k]) )
  }

  y_pp[pc==name_new] <- rnorm(n = sum(pc==name_new),
                       mean = x[pc==name_new,, drop=F] %*% t(beta_new),
                       sd = sqrt(psi_new) )

  return(y_pp)
}
