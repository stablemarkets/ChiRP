#' Function for computing the posterior mode cluster assignment of subjects in test and training data sets.
#'
#' This function takes the \code{cluster_inds} - which is a posterior matrix of cluster assignments output by \code{NDPMix()}, \code{PDPMix()}, and \code{ZDPMix()} and computers the posterior mode cluster assignment while implementing a deterministic re-labeling of subjects.
#'  
#' Please see \url{https://stablemarkets.github.io/ChiRPsite/index.html} for examples and detailed model and parameter descriptions.
#' 
#' Please see \url{https://arxiv.org/abs/1810.09494}, section on hyperparameters and label switching for more information.
#' 
#' @param c_shell set this to \code{cluster_inds} (for either training or testing sets). \code{cluster_inds} is output by \code{NDPMix()}, \code{PDPMix()}, and \code{ZDPMix()}.
#' @return This function returns a list of two objects: \code{adjmat} and \code{class_mem}. For \code{n} subjects, the former is an \code{n} X \code{n} adjacency matrix showing posterior probability subject i and j being clustered together. This can be visualized using a network diagram. \code{class_mem} is a vector of length \code{n} giving posterior mode clsuter membership for each of the \code{n} subjects. This can be used for cluster-specific analysis, for example.
#' @examples
#' # simulate data 
#' 
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
#' # split into training and testing
#' ids <- sample(1:600, size = 500, replace = F )
#' 
#' d$X1 <- scale(d$X1)
#' 
#' d_train <- d[ids,]
#' d_test <- d[-ids, ]
#' 
#' # run zero-inflated model #' 
#' res <- ChiRP::ZDPMix(d_train = d_train, d_test = d_test, formula = Y ~ X1,
#'                      burnin=2000, iter=3000, init_k = 5, phi_y = c(10, 10000))
#' 
#' # compute the posterior model cluster assignment for training subjects
#' train_clus <- ChiRP::cluster_assign_mode(res$cluster_inds$train)
#' @export
cluster_assign_mode <- function(c_shell){
  c_shell <- c_shell
  iter<-ncol(c_shell)
  n <- nrow(c_shell)
  
  adjmat<-matrix(0, nrow=n, ncol=n)
  
  ## compute frequency matrix 
  for(i in 1:iter){
    adjmat_i <- matrix(0, nrow=n, ncol=n)
    class <- c_shell[,i]
    
    for(j in 1:n){
      adjmat_i[j,] <- class==class[j]
    }
    adjmat <- adjmat + adjmat_i
  }
  
  ## choose best cluster assignment
  err_vec <- numeric(length = ncol(c_shell))
  for(i in 1:iter){
    adjmat_i <- matrix(0, nrow=n, ncol=n)
    class <- c_shell[,i]
    for(j in 1:n){
      adjmat_i[j,] <- class==class[j]
    } # compute L2 norm.
    err_vec[i] <- sum((adjmat_i - adjmat)^2)
  }
  
  class_mem <- c_shell[, c(1:ncol(c_shell))[ err_vec==min(err_vec) ] ]
  if(is.matrix(class_mem)) class_mem <- class_mem[,1]
  
  adjmat <- adjmat/iter
  return(list(adjmat=adjmat, class_mem=class_mem))
}