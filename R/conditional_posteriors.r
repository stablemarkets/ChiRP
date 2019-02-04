rcond_post_beta <- function(y, xm, psi, beta_prior_mean, beta_prior_var){
  
  mu_beta <- beta_prior_mean
  v_beta <- diag(beta_prior_var)
  v_beta_inv <- diag(1/beta_prior_var)
  
  xtx <- t(xm)%*%xm
  
  if(length(y)>0 ){
    
    post_cov <- chol2inv(chol( v_beta_inv + (1/psi)*xtx))
    post_mean <- post_cov %*% (v_beta_inv %*% mu_beta + (1/psi)*t(xm)%*%y )
    
  }else{
    post_mean <- mu_beta
    post_cov <- v_beta
  }
  
  draw <- mvtnorm::rmvnorm(n = 1, mean = post_mean, sigma = post_cov)
  return(draw)  
}

cond_post_beta_pdp <- function(beta_vec, y, xm, beta_prior_mean, beta_prior_var){
  
  eta <- xm %*% beta_vec
  eta <- ifelse(eta < -30, -30, ifelse(eta>5, 5, eta))
  
  lk <- sum(dbinom(x = y, size = 1, prob = LaplacesDemon::invlogit(eta), log = T))
  pr <- mvtnorm::dmvnorm(x = beta_vec, beta_prior_mean, diag(beta_prior_var), log = T)
  
  return(lk + pr)
}

rcond_post_psi <- function(beta, y, xm, g1, b1){
  n_k <- length(y)
  if(n_k==0){ 
    shape_k <- g1
    rate_k <- b1
  }else{
    shape_k <- g1 + n_k/2
    rate_k <- .5*(  sum( (y - xm %*% beta)^2 ) )  + b1
  }
  draw <- invgamma::rinvgamma(n = 1, shape = shape_k, rate = rate_k)
  return(draw)
}

rcond_post_mu_x <- function(x_ck, phi_x, lambda, tau ){
  nvec_k <- length(x_ck)
  mu_mean <- (1/(1/tau + nvec_k/phi_x))*(lambda/tau + sum(x_ck)/phi_x)
  mu_sd <- sqrt( (1/tau + nvec_k/phi_x)^(-1) )
  draw <- rnorm(n = 1, mean = mu_mean, sd =  mu_sd)
  return(draw)
}

rcond_post_phi <- function(mu_x, x, g2, b2){
  n_k <- length(x)
  if(n_k==0){ 
    shape_k <- g2
    rate_k <- b2
  }else{
    shape_k <- g2 + n_k/2
    rate_k <- .5*(  sum( (x - mu_x )^2 ) )  + b2
  }
  phi_post <- invgamma::rinvgamma(n = 1, shape = shape_k, rate = rate_k)
  return(phi_post)
}

rcond_post_mu_trt <- function(x_ck){
  draw<-rbeta(n = 1,shape1 = 1 + sum(x_ck), 
              shape2 = 1 + length(x_ck) - sum(x_ck))
  return(draw)
}

cond_post_alpha <- function(alpha, n, K, nvec){
  lik <- lgamma(alpha) - lgamma(n + alpha) + sum(lgamma(nvec + alpha/K)) - K*lgamma(alpha/K)
  pr <- invgamma::dinvgamma(x = alpha, 1, 1,log = T)
  return(lik + pr)
}

cond_post_pz <- function(gamma, y, xp, 
                         gamma_prior_mean,gamma_prior_var){
  
  eta <- xp %*% gamma
  eta <- ifelse(eta < -10, -10, ifelse(eta>10, 10, eta))
  p <- LaplacesDemon::invlogit(eta)
  
  lk <- sum(dbinom(x = y, size = 1, prob = p, log = T))
  pr <- mvtnorm::dmvnorm(x = gamma, gamma_prior_mean, diag(gamma_prior_var), log = T)
  
  return(lk + pr)
}