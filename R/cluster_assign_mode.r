#' Function for computing the posterior mode of cluster assignment.
#'
#' This function corrects for label switching and outputs posterior mode cluster assignment. Label switching correction is deterministic - finding the posterior adjacency matrix that is closest to the posterior mode adjacency matrix in the L2 sense.
#' @param d_train Input dataset.
#' @keywords Label Switching
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
    }
    err_vec[i] <- sum((adjmat_i - adjmat)^2)
  }
  
  class_mem <- c_shell[, c(1:ncol(c_shell))[ err_vec==min(err_vec) ] ]
  if(is.matrix(class_mem)) class_mem <- class_mem[,1]
  
  adjmat <- adjmat/iter
  return(list(adjmat=adjmat, class_mem=class_mem))
}