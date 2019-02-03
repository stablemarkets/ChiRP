#' Function for posterior sampling of DP mixture of zero-inflated regressions.
#'
#' This function takes input dataset and performs posterior sampling. It returns posterior predictive draws, posterior clustering, and posterior parameter draws.
#' @param d_train Input dataset.
#' @keywords Dirichlet
#' @export
PDPMix<-function(d_train, formula, d_test, burnin, iter,
                 phi_y=c(shape=5, rate=1000),
                 beta_prior_mean=NULL, beta_prior_var=NULL,
                 prop_sigma_b = diag(rep(.025, nparams)),
                 init_k=10, beta_var_scale=1000, mu_scale=1, tau_scale=1){
  
  ###------------------------------------------------------------------------###
  #### 0 - Parse User Inputs                                                ####
  ###------------------------------------------------------------------------###
  x <- all.vars(formula[[3]]) # covariate names
  y <- all.vars(formula[[2]]) # outcome name
  
  xt <- model.matrix(data=d_test,
                     object= as.formula(paste0('~ ',paste0(x, collapse = '+'))))
  nt <- nrow(xt)
  
  x_type <- vector(length = length(x))
  for(v in 1:length(x)){
    if( length(unique(d_train[, x[v] ]))==2 ){
      x_type[v] <- 'binary'
    }else if(length(unique(d_train[,x[v]]))>2){
      x_type[v] <- 'numeric'
    }
  }
  
  store_l <- iter - burnin
  y <- d_train[,y]
  x_names <- x
  x <- model.matrix(data=d_train, object = formula )
  
  nparams <- ncol(x)
  n<-nrow(x)
  
  xall_names <- x_names
  nparams_all <- length(xall_names)
  xall <- d_train[,xall_names]
  
  all_type_map <- c(x_type)
  names(all_type_map) <- c(x_names)
  xall_type <- all_type_map[as.character(xall_names)]
  
  xall_names_bin <- xall_names[xall_type == 'binary']
  xall_names_num <- xall_names[xall_type == 'numeric']
  
  # parse numeric versus binary variables for outcome model
  n_bin_p <- sum(x_type=='binary')
  n_num_p <- sum(x_type=='numeric')
  
  ## calibrate priors
  if(is.null(beta_prior_mean)){
    reg <- glm(data=d_train, formula = formula, family = binomial('probit'))
    beta_prior_mean <- reg$coefficients
  }
  
  if(is.null(beta_prior_var)){
    reg <- glm(data=d_train, formula = formula, family = binomial('probit'))
    beta_prior_var <- beta_var_scale*diag(vcov(reg))
  }
  
  g1=phi_y[1]
  b1=phi_y[2]
  
  prior_means <- apply(x[,xall_names_num, drop=F], 2, mean)
  names(prior_means) <- xall_names_num
  
  prior_var <- mu_scale*apply(x[,xall_names_num, drop=F], 2, var)
  names(prior_var) <- xall_names_num
  
  g2 <- rep(2, n_num_p)
  names(g2) <- xall_names_num
  
  b2 <- tau_scale*apply(x[,xall_names_num, drop=F], 2, var)
  names(b2) <- xall_names_num
  
  K <- init_k
  a <- 1
  
  ###------------------------------------------------------------------------###
  #### 1 - Create Shells for storing Gibbs Results                          ####
  ###------------------------------------------------------------------------###
  
  alpha <- numeric(length = store_l)
  
  curr_class_list <- 1:K
  class_names <- paste0('c',curr_class_list)
  
  c_shell<-matrix(NA, nrow=n, ncol=1) # shell for indicators
  
  class_shell <- matrix(NA, nrow=n, ncol=store_l)
  
  beta_shell <- matrix(NA, nrow = nparams, ncol = K)
  colnames(beta_shell) <- class_names
  rownames(beta_shell) <- colnames(x)
  
  c_b_shell <- vector(mode = 'list', length = n_bin_p) # -1 for trt
  names(c_b_shell) <- xall_names_bin
  
  c_n_shell <- vector(mode = 'list', length = n_num_p)
  names(c_n_shell) <- xall_names_num
  
  for( p in xall_names_bin ){
    c_b_shell[[p]] <- matrix(NA, nrow = 1, ncol = K)
    
    c_b_shell[[p]][1,] <- rep(.5, K)
    
    colnames(c_b_shell[[p]]) <- class_names
  }
  
  for( p in xall_names_num ){
    c_n_shell[[p]] <- list(matrix(NA, nrow=1, ncol=K),
                           matrix(NA, nrow=1, ncol=K))
    
    c_n_shell[[p]][[2]][1, ] <- rep(0, K)
    c_n_shell[[p]][[2]][1, ] <- rep(20^2, K)
    
    colnames(c_n_shell[[p]][[1]]) <- class_names
    colnames(c_n_shell[[p]][[2]]) <- class_names
  }
  
  c_b_new <- numeric(length = n_bin_p)
  names(c_b_new) <- xall_names_bin
  
  c_n_new <- matrix(NA, nrow = 2, ncol = n_num_p)
  colnames(c_n_new) <- xall_names_num
  
  c_b_new_test <- numeric(length = n_bin_p)
  names(c_b_new_test) <- xall_names_bin
  
  c_n_new_test <- matrix(NA, nrow = 2, ncol = n_num_p)
  colnames(c_n_new_test) <- xall_names_num
  
  predictions <- vector(mode = 'list', length = 2)
  predictions[[1]] <- matrix(nrow = n, ncol = store_l)
  predictions[[2]] <- matrix(nrow = nt, ncol = store_l)
  names(predictions) <- c('train','test')
  
  cluster_inds <- vector(mode = 'list', length = 2)
  cluster_inds[[1]] <- matrix(nrow = n, ncol = store_l)
  cluster_inds[[2]] <- matrix(nrow = nt, ncol = store_l)
  names(cluster_inds) <- c('train','test')
  
  ###------------------------------------------------------------------------###
  #### 2 - Set Initial Values                                               ####
  ###------------------------------------------------------------------------###
  
  # initially assign everyone to one of K clusters randomly with uniform prob
  c_shell[,1] <- sample(x = class_names, size = n,
                        prob = rep(1/K,K), replace = T)
  
  beta_start <- matrix(0, nrow=nparams, ncol=K)
  colnames(beta_start) <- class_names
  
  for(k in class_names){
    beta_shell[,k] <- beta_prior_mean
  }
  
  for(p in xall_names_bin ){
    c_b_shell[[p]][1,] <- rep(1, K)
  }
  
  for(p in xall_names_num ){
    c_n_shell[[p]][[2]][ 1 , ] <- rep(0, K)
    c_n_shell[[p]][[2]][ 1 , ] <- rep(20^2, K)
  }
  
  
  ###------------------------------------------------------------------------###
  #### 3 - Gibbs Sampler                                                    ####
  ###------------------------------------------------------------------------###
  message_seq <- seq(2, iter, round(iter/20) )
  cat(paste0('Iteration ',2,' of ', iter,' (burn-in) \n'))
  message_step <- (iter)/10
  for(i in 2:iter){
    
    ## print message to user
    if( i %% message_step == 1 ) cat(paste0('Iteration ',i,' of ',
                                            iter,
                                            ifelse(i>burnin,' (sampling) ',
                                                   ' (burn-in)') ,
                                            '\n'))
    
    # compute size of each existing cluster
    class_ind<-c_shell[,1]
    nvec<-table(factor(class_ind,levels = class_names) )
    
    ###----------------------------------------------------------------------###
    #### 3.1 - update Concentration Parameter                               ####
    ###----------------------------------------------------------------------###
    
    a_star <- rnorm(n = 1, mean =  a, sd = 1)
    if(a_star>0){
      r_num <- cond_post_alpha(a_star, n, K, nvec)
      r_denom <- cond_post_alpha(a, n, K, nvec)
      r <- exp(r_num - r_denom)
      accept <- rbinom(1,1, min(r,1) )==1
      if(accept){ a <- a_star }
    }
    
    ###----------------------------------------------------------------------###
    #### 3.2 - update parameters conditional on classification              ####
    ###----------------------------------------------------------------------###
    K <- length(unique(class_names))
    for(k in class_names){ # cycle through existing clusters - update parameters
      
      ck_ind <- class_ind == k
      
      y_ck <- y[ck_ind]
      
      x_ck <- x[ck_ind,, drop=FALSE]
      
      beta_shell[ ,k] <- metrop_hastings(x_0 = beta_shell[ ,k], iter = 1,
                                         log_post_density = cond_post_beta_pdp,
                                         y= y, xm = x_ck,
                                         prop_sigma = prop_sigma_b,
                                         beta_prior_mean = beta_prior_mean,
                                         beta_prior_var = beta_prior_var)[[1]]
      
      for( p in xall_names_bin ){
        c_b_shell[[p]][1, k] <- rcond_post_mu_trt(x_ck = x_ck[, p])
      }
      
      for( p in xall_names_num ){
        c_n_shell[[p]][[1]][1,k] <- rcond_post_mu_x(x_ck = x_ck[, p],
                                                    phi_x = c_n_shell[[p]][[2]][1,k],
                                                    lambda = prior_means[p],
                                                    tau = prior_var[p])
        
        c_n_shell[[p]][[2]][1,k] <- rcond_post_phi(mu_x = c_n_shell[[p]][[1]][1, k],
                                                   x = x_ck[ , p],
                                                   g2 = g2[p], b2 = b2[p])
      }
      
    }
    
    ###----------------------------------------------------------------------###
    #### 3.3 - update classification conditional on parameters              ####
    ###----------------------------------------------------------------------###
    
    ## draw parameters for potentially new cluster from priors
    ### draw beta, gamma, phi, covariate parameters
    
    beta_new <- mvtnorm::rmvnorm(n = 1, mean = beta_prior_mean,
                                 sigma = diag(beta_prior_var) )

    for(p in xall_names_bin){
      c_b_new[p] <- rbeta(n = 1, shape1 = 1, shape2 = 1)
    }
    
    for(p in xall_names_num){
      c_n_new[1, p] <- rnorm(n = 1, prior_means[p], sqrt(prior_var[p]) )
      c_n_new[2, p] <- invgamma::rinvgamma(n = 1, shape = g2[p], rate = b2[p])
    }
    
    name_new <- paste0('c', max(as.numeric(substr(class_names, 2,10))) + 1)
    
    ## compute post prob of membership to existing and new cluster,
    ## then classify
    
    weights <- class_update_train_pdp(n = n, K = length(class_names), alpha = a,
                                      name_new =  name_new,
                                      uniq_clabs = colnames(beta_shell),
                                      clabs = class_ind,
                                      
                                      y = y, x = x,
                                      x_cat_shell = c_b_shell , x_num_shell = c_n_shell,
                                      
                                      ## col number of each covar...index with base 0
                                      cat_idx=xall_names_bin, num_idx=xall_names_num,
                                      
                                      beta_shell = beta_shell,
                                      
                                      beta_new=beta_new,
                                      cat_new = c_b_new, num_new=c_n_new)
    
    c_shell[,1] <- apply(weights, 1, FUN = sample,
                         x=c(colnames(beta_shell), name_new), size=1, replace=T)
    
    ###----------------------------------------------------------------------###
    #### 4.0 - update shell dimensions as clusters die/ are born            ####
    ###----------------------------------------------------------------------###
    
    new_class_names <- unique(c_shell[,1])
    K_new <- length(new_class_names)
    
    # Grow shells by adding labels/parameters of new clusters
    new_class_name <- setdiff(new_class_names, class_names)
    if(length(new_class_name)>0){ # new cluster was born
      
      beta_shell <- cbind(beta_shell, t(beta_new) )
      colnames(beta_shell)[K+1] <- new_class_name
      
      for(p in xall_names_bin ){
        c_b_shell[[p]] <- cbind( c_b_shell[[p]], c_b_new[p] )
        colnames(c_b_shell[[p]])[K+1] <- new_class_name
      }
      
      for(p in xall_names_num ){
        c_n_shell[[p]][[1]] <- cbind(c_n_shell[[p]][[1]], c_n_new[1, p])
        colnames(c_n_shell[[p]][[1]])[K+1] <- new_class_name
        
        c_n_shell[[p]][[2]] <- cbind(c_n_shell[[p]][[2]], c_n_new[2, p])
        colnames(c_n_shell[[p]][[2]])[K+1] <- new_class_name
      }
    }
    
    # shrink shells by deleting labels of clusters that have died.
    dead_classes <- setdiff(class_names, new_class_names)
    if(length(dead_classes) > 0){ # clusters have died
      beta_shell <- beta_shell[,! colnames(beta_shell) %in% dead_classes, drop=F]
      
      for(p in xall_names_bin){
        c_b_shell[[p]] <- c_b_shell[[p]][, ! colnames(c_b_shell[[p]]) %in% dead_classes, drop=F]
      }
      
      for(p in xall_names_num){
        c_n_shell[[p]][[1]] <- c_n_shell[[p]][[1]][, ! colnames(c_n_shell[[p]][[1]]) %in% dead_classes, drop=F]
        c_n_shell[[p]][[2]] <- c_n_shell[[p]][[2]][, ! colnames(c_n_shell[[p]][[2]]) %in% dead_classes, drop=F]
      }
    }
    
    K <- K_new
    class_names <- new_class_names
    
    if( i > burnin ){
      
      ###--------------------------------------------------------------------###
      #### 5.0 - Predictions on a Training and Test set                     ####
      ###--------------------------------------------------------------------###
      
      name_new_test<-paste0('c', max(as.numeric(substr(class_names,2,10)))+1)
      
      for(p in xall_names_bin){
        c_b_new_test[p] <- rbeta(n = 1, shape1 = 1, shape2 = 1)
      }
      
      for(p in xall_names_num){
        c_n_new_test[1, p] <- rnorm(n = 1, prior_means[p], sqrt(prior_var[p]) )
        c_n_new_test[2, p] <- invgamma::rinvgamma(n = 1,shape = g2[p], rate = b2[p])
      }
      
      beta_new <- mvtnorm::rmvnorm(n = 1, mean = beta_prior_mean,
                                   sigma = diag(beta_prior_var) )

      weights_test <- class_update_test_pdp(n = nt, K = length(class_names) ,
                                            alpha =a,
                                            name_new = name_new_test,
                                            uniq_clabs = colnames(beta_shell),
                                            clabs = c_shell[,1],
                                            x = xt,
                                            x_cat_shell = c_b_shell,
                                            x_num_shell = c_n_shell,
                                            cat_idx = xall_names_bin,
                                            num_idx = xall_names_num,
                                            cat_new = c_b_new_test,
                                            num_new = c_n_new_test)
      
      c_shell_test <- apply(weights_test, 1, FUN = sample,
                            x=c(colnames(beta_shell), name_new_test),
                            size=1, replace=T)
      
      test_pred <- post_pred_draw_test_pdp(n=nt, x=xt, pc=c_shell_test,
                                           beta_shell = beta_shell,
                                           name_new = name_new_test,
                                           beta_new = beta_new)
      
      train_pred<-post_pred_draw_train_pdp(n=n, x=x, pc=c_shell[,1],
                                           beta_shell = beta_shell)
      
      ###--------------------------------------------------------------------###
      #### 6.0 - Store Results                                              ####
      ###--------------------------------------------------------------------###
      
      # cluster indicators
      cluster_inds[['train']][, i-burnin] <- c_shell[,1]
      cluster_inds[['test']][, i-burnin] <- c_shell_test
      
      # predictions
      predictions[['train']][, i-burnin ] <- train_pred
      predictions[['test']][, i-burnin  ] <- test_pred
    }
    
    
  }
  
  results <- list(cluster_inds=cluster_inds,
                  predictions = predictions)
  cat('Sampling Complete.')
  return(results)
}