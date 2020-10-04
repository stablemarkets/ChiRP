#' Function for posterior sampling of DP mixture of zero-inflated regressions.
#'
#' This function takes in a training data.frame and optional testing data.frame and performs posterior sampling. It returns posterior predictions and posterior clustering for training and test sets.
#' The function is built for zero-inflated, but otherwise continuous, outcomes.
#' 
#' Please see \url{https://stablemarkets.github.io/ChiRPsite/index.html}for examples and detailed model and parameter descriptions.
#' 
#' Please see \url{https://arxiv.org/abs/1810.09494} for a methodological reference.
#' 
#' @import stats
#' 
#' @param d_train A \code{data.frame} object with outcomes and model covariates/features. All features must be \code{as.numeric} - either continuous or binary with binary variables coded using \code{1} and \code{0}. Categorical features are not supported. We recommend standardizing all continuous features. NA values are not allowed and each row should represent a single subject, longitudinal data is not supported.
#' @param d_test Optional \code{data.frame} object containing a test set of subjects containing all variables specifed in \code{formula}. All the same rules for \code{d_train} apply to \code{d_test}.
#' @param formula Specified in the usual way, e.g. for \code{p=2} covariates, \code{y ~ x1 + x2}. All covariates - continuous and binary - must be \code{as.numeric} , with binary variables coded as \code{1} or \code{0}. We recommend standardizing all continuous features. NA values are not allowed and each row should represent a single subject, longitudinal data is not supported.
#' @param burnin integer specifying number of burn-in MCMC draws. 
#' @param iter integer greater than \code{burnin} specifying how many total MCMC draws to take.
#' @param init_k Optional. integer specifying the initial number of clusters to kick off the MCMC sampler.
#' @param prop_sigma_z Optional. If you specified \code{p} covariates in \code{formula}, \code{p+1} regression parameters are sampled for the probability of the outcome being zero using a Metropolis step.  \code{prop_sigma_z} is a \code{p+1} by \code{p+1} covariance matrix for the Metropolis proposal distribution.
#' @param phi_y Optional. Length two \code{as.numeric} vector specifying the shape and rate, respectively, of the Inverse Gamma hyper-prior placed on the outcome variance. 
#' @param beta_prior_mean Optional. If there are \code{p} covariates, it is a length \code{p+1} \code{as.numeric} vector specifying mean of the Gaussian prior on the outcome model's conditional mean parameter vector. Default is regression coefficients from running OLS on positive outcomes.
#' @param beta_prior_var Optional. If there are \code{p} covariates, a length \code{p+1} \code{as.numeric} vector specifying variance of the Gaussian prior on the outcome model's conditional mean parameter vector. The full covarince of the prior is set to be diagonal. This vector specifies the diagonal enteries of this prior covariance. Default is estimated variances from running OLS on positive outcomes.
#' @param beta_var_scale Optional. A multiplicative constant that scales \code{beta_prior_var}. If you leave \code{beta_prior_mean} and \code{beta_prior_var} at their defaults, This constant toggles how wide new cluster parameters are dispersed around the observed data parameters, larger values implies wider distribution.
#' @param gamma_prior_mean Optional. If there are \code{p} covariates, a length \code{p+1} \code{as.numeric} vector specifying mean of the Gaussian prior on the zero probability logistic model's conditional mean parameter vector. Default is a vector of 0s - i.e., null-centered prior mean.
#' @param gamma_prior_var Optional. If there are \code{p} covariates, a length \code{p+1} \code{as.numeric} vector specifying variance of the Gaussian prior on the zero probability logistic model's conditional mean parameter vector. Default is vector of 2s - moderately flat on the odds ratio scale.
#' @param mu_scale Optional. An numeric, scalar constant that controls how widely distributed new cluster continuous covariate means are distributed around the empirical covariate mean. Specifically, all continuous covariates are assumed to have Gaussian likelihood with Gaussian prior on their means. \code{mu_scale=2} specifies that the variance of the Gaussian prior is twice as large as the empirical variance.
#' @param tau_scale Optional. An numeric, scalar constant that controls how widely distributed new cluster continuous covariate variances are distributed around the empirical variance. Specifically, all continuous covariates are assumed to have Gaussian likelihood with Inverse Gamma prior on their variance. \code{tau_scale=2} specifies that the rate of the InvGamma prior is twice as large as the empirical variance.
#' @return Returns \code{predictions$train} and \code{cluster_inds$train}. \code{predictions$train} returns an \code{nrow(d_train)} by \code{iter - burnin} matrix of posterior predictions. \code{cluster_inds$train} returns an \code{nrow(d_train)} by \code{iter - burnin} matrix of cluster assignment indicators, which can be input into the function \code{cluster_assign_mode()} to compute posterior mode assignment. \code{predictions$test} and \code{cluster_inds$test} are returned if \code{d_test} is specified.
#' @keywords Dirichlet Process Zero Inflated
#' @examples
#' set.seed(1)
#' n<-200 ## generate from clustered, skewed, data distribution
#' X11 <- rnorm(n = n, mean = 10, sd = 3)
#' X12 <- rnorm(n = n, mean = 0, sd = 2)
#' X13 <- rnorm(n = n, mean = -10, sd = 4)
#' 
#' Y1 <- rnorm(n = n, mean = 100 + .5*X11, 20)*(1-rbinom(n, 1, prob = pnorm( -10 + 1*X11 ) ))
#' Y2 <- rnorm(n = n, mean = 200 + 1*X12, 30)*(1-rbinom(n, 1, prob = pnorm( 1 + .05*X12 ) ))
#' Y3 <- rnorm(n = n, mean = 300 + 2*X13, 40)*(1-rbinom(n, 1, prob = pnorm( -3 -.2*X13 ) ))
#' 
#' d <- data.frame(X1=c(X11, X12, X13), Y = c(Y1, Y2, Y3))
#' 
#' d$X1 <- scale(d$X1)
#' 
#' ids <- sample(1:600, size = 500, replace = FALSE )
#' d_train <- d[ids,]
#' d_test <- d[-ids, ]
#' 
#' res <- ChiRP::ZDPMix(d_train = d_train, d_test = d_test, formula = Y ~ X1,
#'                      burnin=100, iter=200, init_k = 5, phi_y = c(10, 10000))
#' @export
ZDPMix<-function(d_train, formula, d_test=NULL, burnin=100, iter=1000,
                 phi_y=c(shape=5, rate=1000),
                 beta_prior_mean=NULL, beta_prior_var=NULL,
                 gamma_prior_mean=NULL, gamma_prior_var=NULL,
                 init_k=10, beta_var_scale=1000, mu_scale=1, tau_scale=1,
                 prop_sigma_z = diag(rep(.025, nparams)) ){

  ###------------------------------------------------------------------------###
  #### 0 - Parse User Inputs                                                ####
  ###------------------------------------------------------------------------###
  
  # error checking user inputs
  if( missing(d_train) ){ stop("ERROR: must specify a training data.frame.") }
  
  x <- all.vars(formula[[3]]) # covariate names
  y <- all.vars(formula[[2]]) # outcome name
  
  nparams <- length(x) + 1
  
  func_args<-mget(names(formals()),sys.frame(sys.nframe()))
  error_check(func_args, 'ZDP')
  
  if(!is.null(d_test)){
    xt <- model.matrix(data=d_test,
                       object= as.formula(paste0('~ ',paste0(x, collapse = '+'))))
    nt <- nrow(xt)
  }
  
  x_type <- vector(length = length(x))
  for(v in 1:length(x)){
    if( length(unique(d_train[, x[v] ]))==2 ){
      x_type[v] <- 'binary'
    }else if(length(unique(d_train[,x[v]]))>2){
      x_type[v] <- 'numeric'
    }
  }

  store_l <- iter - burnin
  y <- d_train[[y]]
  z <- as.numeric(y==0)
  x_names <- x
  x <- model.matrix(data=d_train, object = formula )

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
    reg <- lm(data=d_train[y>0,], formula = formula)
    beta_prior_mean <- reg$coefficients
  }

  if(is.null(beta_prior_var)){
    reg <- lm(data=d_train[y>0,], formula = formula)
    beta_prior_var <- beta_var_scale*diag(vcov(reg))
  }

  if(is.null(gamma_prior_mean)){
    gamma_prior_mean <- rep(0, nparams)
  }

  if(is.null(gamma_prior_var)){
    gamma_prior_var <- rep(2, nparams)
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

  psi_shell <- matrix(NA, nrow= 1, ncol = K) # shell for Y's cond. var
  colnames(psi_shell) <- class_names

  beta_shell <- matrix(NA, nrow = nparams, ncol = K)
  colnames(beta_shell) <- class_names
  rownames(beta_shell) <- colnames(x)

  gamma_shell <- matrix(NA, nrow = nparams, ncol = K)
  colnames(gamma_shell) <- class_names
  rownames(gamma_shell) <- colnames(x)

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

  if(!is.null(d_test)){
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
  }else{
    predictions <- vector(mode = 'list', length = 1)
    predictions[[1]] <- matrix(nrow = n, ncol = store_l)
    names(predictions) <- c('train')
    
    cluster_inds <- vector(mode = 'list', length = 1)
    cluster_inds[[1]] <- matrix(nrow = n, ncol = store_l)
    names(cluster_inds) <- c('train')
  }
  
  ###------------------------------------------------------------------------###
  #### 2 - Set Initial Values                                               ####
  ###------------------------------------------------------------------------###

  # initially assign everyone to one of K clusters randomly with uniform prob
  c_shell[,1] <- sample(x = class_names, size = n,
                        prob = rep(1/K,K), replace = T)

  beta_start <- matrix(0, nrow=nparams, ncol=K)
  colnames(beta_start) <- class_names

  gamma_start <- matrix(0, nrow=nparams, ncol=K)
  colnames(gamma_start) <- class_names

  for(k in class_names){
    beta_shell[,k] <- beta_prior_mean
    gamma_shell[,k] <- gamma_start[,k]
  }

  for(p in xall_names_bin ){
    c_b_shell[[p]][1,] <- rep(1, K)
  }

  for(p in xall_names_num ){
    c_n_shell[[p]][[2]][ 1 , ] <- rep(0, K)
    c_n_shell[[p]][[2]][ 1 , ] <- rep(20^2, K)
  }

  psi_shell[1,1:K] <- rep(500000, K)

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

      psi_shell[1, k]<-rcond_post_psi(beta = beta_shell[,k, drop=F],
                                      y = y_ck[y_ck>0], xm = x_ck[y_ck>0, ,drop=F],
                                      g1 = g1, b1 = b1)

      beta_shell[,k] <- rcond_post_beta(y=y_ck[y_ck>0],
                                        xm=x_ck[y_ck>0, , drop=F],
                                        psi = psi_shell[1,k],
                                        beta_prior_mean = beta_prior_mean,
                                        beta_prior_var = beta_prior_var)

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

      gamma_shell[ ,k] <- metrop_hastings(x_0 = gamma_shell[ ,k], iter = 1,
                                          log_post_density = cond_post_pz,
                                          y= y_ck==0, xp = x_ck,
                                          prop_sigma = prop_sigma_z,
                                          gamma_prior_mean = gamma_prior_mean,
                                          gamma_prior_var = gamma_prior_var)[[1]]
    }

    ###----------------------------------------------------------------------###
    #### 3.3 - update classification conditional on parameters              ####
    ###----------------------------------------------------------------------###

    ## draw parameters for potentially new cluster from priors
    ### draw beta, gamma, phi, covariate parameters

    beta_new <- mvtnorm::rmvnorm(n = 1, mean = beta_prior_mean,
                                 sigma = diag(beta_prior_var) )
    psi_new <- invgamma::rinvgamma(n = 1, shape = g1, rate = b1)
    gamma_new <- mvtnorm::rmvnorm(n = 1, mean = gamma_prior_mean,
                                  sigma = diag(gamma_prior_var))

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

    weights <- class_update_train(n = n, K = length(class_names), alpha = a,
                                  name_new =  name_new,
                                  uniq_clabs = colnames(beta_shell),
                                  clabs = class_ind,

                                  y = y, x = x, z = z,
                                  x_cat_shell = c_b_shell , x_num_shell = c_n_shell,

                                  ## col number of each covar...index with base 0
                                  cat_idx=xall_names_bin, num_idx=xall_names_num,

                                  beta_shell = beta_shell,
                                  psi_shell = psi_shell,
                                  gamma_shell = gamma_shell,

                                  beta_new=beta_new, psi_new=psi_new,
                                  cat_new = c_b_new, num_new=c_n_new,
                                  gamma_new=gamma_new)

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

      psi_shell <- cbind(psi_shell, psi_new)
      colnames(psi_shell)[K+1] <- new_class_name

      gamma_shell <- cbind(gamma_shell, t(gamma_new) )
      colnames(gamma_shell)[K+1] <- new_class_name

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
      gamma_shell <- gamma_shell[,! colnames(gamma_shell) %in% dead_classes, drop=F]

      psi_shell <- psi_shell[, ! colnames(psi_shell) %in% dead_classes, drop=F]

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
      
      if(!is.null(d_test)){
        ###------------------------------------------------------------------###
        #### 5.0 - Predictions on a Training and Test set and Store Results ####
        ###------------------------------------------------------------------###
        
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
        psi_new <- invgamma::rinvgamma(n = 1, shape = g1, rate = b1)
        gamma_new <- mvtnorm::rmvnorm(n = 1, mean = gamma_prior_mean,
                                      sigma = diag(gamma_prior_var))
        
        weights_test <- class_update_test(n = nt, K = length(class_names) ,
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
        
        test_pred <- post_pred_draw_test(n=nt, x=xt, pc=c_shell_test,
                                         beta_shell = beta_shell,
                                         psi_shell=psi_shell,
                                         gamma_shell=gamma_shell,
                                         name_new = name_new_test,
                                         beta_new = beta_new, psi_new=psi_new,
                                         gamma_new=gamma_new)
        
        train_pred<-post_pred_draw_train(n=n, x=x, pc=c_shell[,1],
                                         beta_shell = beta_shell,
                                         psi_shell = psi_shell,
                                         gamma_shell = gamma_shell)

        cluster_inds[['train']][, i-burnin] <- c_shell[,1]
        cluster_inds[['test']][, i-burnin] <- c_shell_test
        
        predictions[['train']][, i-burnin ] <- train_pred
        predictions[['test']][, i-burnin  ] <- test_pred  
     
      }else{
        train_pred<-post_pred_draw_train(n=n, x=x, pc=c_shell[,1],
                                         beta_shell = beta_shell,
                                         psi_shell = psi_shell,
                                         gamma_shell = gamma_shell)

        cluster_inds[['train']][, i-burnin] <- c_shell[,1]
        predictions[['train']][, i-burnin ] <- train_pred
      }

    }


  }

  results <- list(cluster_inds=cluster_inds,
                  predictions = predictions)
  cat('Sampling Complete.')
  return(results)
}
