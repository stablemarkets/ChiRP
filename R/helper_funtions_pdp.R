class_update_train_pdp <- function(n, K , alpha, name_new, uniq_clabs, clabs,
                                   y, x, x_cat_shell, x_num_shell,
                                   cat_idx, num_idx,
                                   beta_shell, beta_new, cat_new, num_new){
  
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
    
    eta <- x %*% beta_shell[,k,drop=F]
    eta <- ifelse(eta < -10, -10, ifelse(eta>10, 10, eta))
    
    lk_exist <- lk_exist + dbinom(y, 1,prob = LaplacesDemon::invlogit(eta), T )
    
    
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
  
  
  eta <- x %*% t(beta_new)
  eta <- ifelse(eta < -10, -10, ifelse(eta>10, 10, eta))
  
  lk_exist <- lk_exist + dbinom(y, 1,prob = LaplacesDemon::invlogit( eta ), T )
  
  c_shell[,name_new] <- lk_new + pr_new
  
  weights <- t(apply(c_shell, 1, function(x) exp(x)/sum(exp(x))  ))
  
  return(weights)
}

class_update_test_pdp <- function(n, K , alpha, name_new, uniq_clabs, clabs, x,
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

post_pred_draw_train_pdp <- function(n, x, pc, beta_shell){
  
  y_pp <- numeric(length = n)
  
  clabs <- unique(pc)
  
  for( k in clabs ){
    y_pp[pc==k] <- LaplacesDemon::invlogit( x[pc==k,, drop=F] %*% beta_shell[,k,drop=F] )
  }
  
  return(y_pp)
}


post_pred_draw_test_pdp <- function(n, x, pc, beta_shell,
                                    name_new, beta_new){
  
  y_pp <- numeric(length = n)
  
  clabs <- unique(pc)
  
  for( k in setdiff(clabs, name_new) ){
    y_pp[pc==k] <- LaplacesDemon::invlogit( x[pc==k,, drop=F] %*% beta_shell[,k,drop=F] )
  }
  
  y_pp[pc==name_new] <- LaplacesDemon::invlogit( x[pc==name_new, , drop=F] %*% t(beta_new) )
  
  return(y_pp)
}