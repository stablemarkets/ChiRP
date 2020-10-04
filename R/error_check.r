error_check <- function(func_args, model_type){
  
  x <- all.vars(func_args$formula[[3]]) # covariate names
  y <- all.vars(func_args$formula[[2]]) # outcome name
  
  nparams <- length(x) + 1
  
  ## check init_k - initial number of clusters
  if(class(func_args$init_k)!='numeric' | length(func_args$init_k)!=1 ){
    stop("ERROR: init_k must be length 1 numeric integer such that 0 < init_k <= n .")
  }else if(func_args$init_k<0){
    stop("ERROR: init_k must be 0 < init_k <= n .")
  }else if( func_args$init_k %% 1 !=0 ){
    stop('ERROR: init_k must be an integer.')
  }else if(  func_args$init_k > nrow(func_args$d_train) ){
    stop('ERROR: init_k cannot be greater than sample size.')
  }
    
  ## Check burnin and iter settings
  if( ! (func_args$burnin %% 1 ==0 & func_args$iter %% 1 ==0) ){
    stop("ERROR: burnin and iter must be integers.")
  }else if( func_args$burnin < 0 | func_args$iter < 0){
    stop("ERROR: burnin and iter must be positive.")
  }else if(func_args$burnin==0){
    warning("Warning (non-fatal): burnin==0. Should specify 0 < burnin < iter.")
  }else if(func_args$iter==0){
    stop('ERROR: iter==0. Should specify 0 < burnin < iter.')
  }else if(! func_args$burnin < func_args$iter){
    stop("ERROR: iter must be greater than burnin.")      
  }
  
  # check training data
  if( !is.data.frame(func_args$d_train) ){
    stop("ERROR: d_train is not a data.frame object.")
  }else if( nrow(func_args$d_train)==0 ){
    stop("ERROR: d_train has no rows.")
  }else if( ncol(func_args$d_train)==0 ){
    stop("ERROR: d_train has no columns.")
  }else if( nrow(func_args$d_train)<=2  ){
    stop('ERROR: d_train should contain at least 3 observations.')
  }else if( nrow(func_args$d_train) < func_args$init_k ){
    stop("ERROR: nrow(d_train) < init_k.")
  }else if( max( ! x %in% colnames(func_args$d_train) )==1 ){
    stop("ERROR: Some covariate specified in formula are not in d_train")
  }else if(  max(apply(func_args$d_train[,x, drop=F], 2, class)!="numeric") ){
    stop("ERROR: All columns of d_train must be numeric.")
  }
  
  for( var in x ){
    unvals <- unique(func_args$d_train[,var])
    nu <- length(unvals)
    if( nu==2 & sum(unvals %in% c(0,1)) != 2 ){
      stop(paste0('ERROR: ',var," in d_train is binary, but not coded as 0-1"))
    }
  }
  
  # check testing data set if supplied.
  if(! is.null(func_args$d_test) ){
    
    if( !is.data.frame(func_args$d_test) ){
      stop("ERROR: d_test is not a data.frame object.")
    }else if( nrow(func_args$d_test)==0 ){
      stop("ERROR: d_test has no rows.")
    }else if( ncol(func_args$d_test)==0 ){
      stop("ERROR: d_test has no columns.")
    }else if( max( ! x %in% colnames(func_args$d_test) )==1 ){
      stop("ERROR: Some covariate specified in formula are not in d_test.")
    }else if(  max(apply(func_args$d_test[,x, drop=F], 2, class)!="numeric") ){
      stop("ERROR: All columns of d_test must be numeric.")
    }
    
  }
  
  if(!is.null(func_args$beta_prior_mean)){
    if( class(func_args$beta_prior_mean) != 'numeric' ){
      stop('ERROR: beta_prior_mean should be a numeric vector.')
    }else if(length(func_args$beta_prior_mean) != (length(x)+1) ){
      stop("ERROR: beta_prior_mean should be length p+1, where p is number of covariates.")
    }
  }

  if(!is.null(func_args$beta_prior_var)){
    if( class(func_args$beta_prior_var) != 'numeric' ){
      stop('ERROR: beta_prior_var should be a numeric vector.')
    }else if(length(func_args$beta_prior_var) != (length(x)+1) ){
      stop("ERROR: beta_prior_var should be length p+1, where p is number of covariates.")
    }else if( max(func_args$beta_prior_var<0)==1 ){
      stop("ERROR: beta_prior_var should contain non-zero values since it is a vector of variances.")
    }
  }
  
  if(length(func_args$beta_var_scale)!=1){
    stop("ERROR: beta_var_scale should be a scalar, length 1 numeric.")
  }else if(class(func_args$beta_var_scale)!='numeric'){
    stop("ERROR: beta_var_scale should be numeric.")
  }else if(func_args$beta_var_scale<0){
    stop("ERROR: beta_var_scale should be >0.")
  }
  
  if(!is.numeric(func_args$mu_scale) | length(func_args$mu_scale)!=1 | func_args$mu_scale<=0 ){
    stop("ERROR: mu_scale should be length 1 numeric values greater than 0.")
  }
  
  if(model_type!='fDP'){
    if(!is.numeric(func_args$tau_scale) | length(func_args$tau_scale)!=1 | func_args$tau_scale<=0 ){
      stop("ERROR: tau_scale should be length 1 numeric values greater than 0.")
    }
  }
  
  if(! class(func_args$d_train[[y]]) %in% c('numeric', 'integer') ){
    stop("ERROR: outcome is not numeric")
  }
  
  ## model-specific checks
  if(model_type=='ZDP'){
    
    if( length(func_args$phi_y)!=2 ){
      stop("ERROR: phi_y should be a length 2 vector of positive real numbers.")
    }else if( max(func_args$phi_y<0)==1 ){
      stop("ERROR: phi_y contains negative values. Should be positive real numbers")
    }
    
    if( sum(func_args$d_train[,y]==0)==0 ){
      stop("ERROR: outcome does not contain any zeros. ZDPMix() is inappropriate.")
    }
    
    if(!is.null(func_args$gamma_prior_mean)){
      if( class(func_args$gamma_prior_mean) != 'numeric' ){
        stop('ERROR: gamma_prior_mean should be a numeric vector.')
      }else if(length(func_args$gamma_prior_mean) != (length(x)+1) ){
        stop("ERROR: gamma_prior_mean should be length p+1, where p is number of covariates.")
      }
    }
    
    if(!is.null(func_args$gamma_prior_var)){
      
      if( class(func_args$gamma_prior_var) != 'numeric' ){
        stop('ERROR: gamma_prior_var should be a numeric vector.')
      }else if(length(func_args$gamma_prior_var) != (length(x)+1) ){
        stop("ERROR: gamma_prior_var should be length p+1, where p is number of covariates.")
      }else if( max(func_args$gamma_prior_var<0)==1 ){
        stop("ERROR: gamma_prior_var should contain non-zero values since it is a vector of variances.")
      }
      
    }
    
    if( is.null(func_args$prop_sigma_z) ){
      stop('ERROR: Proposal covariance is null.')
    }else if( class(func_args$prop_sigma_z)!='matrix' ){
      stop('ERROR: prop_sigma_z should be a (p+1) X (p+1) numeric covariance matrix')
    }else if(!is.numeric(func_args$prop_sigma_z)){
      stop("ERROR: prop_sigma_z is not numeric.")
    }else if( !isSymmetric(func_args$prop_sigma_z) ){
      stop("ERROR: prop_sigma_z is not symmetric.")
    }else if( det(func_args$prop_sigma_z)==0 ){
      stop("ERROR: prop_sigma_z is not positive definite.")
    }
    
    
  }else if(model_type=='NDP'){
    
    if( length(func_args$phi_y)!=2 ){
      stop("ERROR: phi_y should be a length 2 vector of positive real numbers.")
    }else if( max(func_args$phi_y<0)==1 ){
      stop("ERROR: phi_y contains negative values. Should be positive real numbers")
    }
    
  }else if(model_type=='fDP'){
    
    if( !is.null(func_args$phi_y) & length(func_args$phi_y)!=2 ){
      stop("ERROR: phi_y should be a length 2 vector of positive real numbers.")
    }else if( max(func_args$phi_y<0)==1 ){
      stop("ERROR: phi_y contains negative values. Should be positive real numbers")
    }
    
    if( !is.null(func_args$tau_x) & length(func_args$tau_x)!=2 ){
      stop("ERROR: tau_x should be a length 2 vector of positive real numbers.")
    }else if( max(func_args$phi_y<0)==1 ){
      stop("ERROR: tau_x contains negative values. Should be positive real numbers")
    }
    
  }else if(model_type=='PDP'){
    
    uy <- unique(func_args$d_train[,y])
    if(length(uy)!=2){
      stop("ERROR: outcome is not binary.")
    }else if( sum(uy %in% c(0,1))!=2  ){
      stop("ERROR: outcome must be coded as numeric vector with 1 and 0. 1 indicates events.")
    }
    
    if( is.null(func_args$prop_sigma_b) ){
      stop('ERROR: Proposal covariance is null.')
    }else if( class(func_args$prop_sigma_b)!='matrix' ){
      stop('ERROR: prop_sigma_b should be a (p+1) X (p+1) numeric covariance matrix')
    }else if(!is.numeric(func_args$prop_sigma_b)){
      stop("ERROR: prop_sigma_b is not numeric.")
    }else if( !isSymmetric(func_args$prop_sigma_b) ){
      stop("ERROR: prop_sigma_b is not symmetric.")
    }else if( det(func_args$prop_sigma_b)==0 ){
      stop("ERROR: prop_sigma_b is not positive definite.")
    }
    
  }
  
}


error_check_credible_gradient = function(func_args){
  
  if( class(func_args$x) != 'numeric' ){
    stop("ERROR: x is not numeric")
  }else if( class(func_args$post_draws) != "matrix" ){
    stop("ERROR: post_draws is not a matrix. Should be matrix with dimensions length(x) by #posterior draws ")
  }else if( !is.character(func_args$col_gradient) | !is.character(func_args$col_mean_line) ){
    stop("ERROR: either col_gradient or col_mean_line is not a character vector specifying a particular color")
  }
  
}

