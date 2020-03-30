#' Function for plotting posterior regression line 
#'
#' This function takes as input some continuous covariate \code{x} and output of predictions from either NDPMix, ZDPMix, fDPMix, etc. It outputs a credible band and posterior median line of the predictions across \code{x}
#' 
#' Please see \url{https://stablemarkets.github.io/ChiRPsite/index.html} for examples and detailed model and parameter descriptions.
#' 
#' @import stats
#' 
#' @param x A \code{data.frame} object with outcomes and model covariates/features. All features must be \code{as.numeric} - either continuous or binary with binary variables coded using \code{1} and \code{0}. Categorical features are not supported. We recommend standardizing all continuous features. NA values are not allowed and each row should represent a single subject, longitudinal data is not supported.
#' @param post_draws Optional \code{data.frame} object containing a test set of subjects containing all variables specifed in \code{formula}. The same rules that apply to \code{d_train} apply to \code{d_test}.
#' @param formula Specified in the usual way, e.g. for \code{p=2} covariates, \code{y ~ x1 + x1}. All covariates - continuous and binary - must be \code{as.numeric} , with binary variables coded as \code{1} or \code{0}. We recommend standardizing all continuous features. NA values are not allowed and each row should represent a single subject, longitudinal data is not supported.
#' @param burnin interger specifying number of burn-in MCMC draws. 
#' @return overlays credible band and posterior median line.
#' @keywords Dirichlet Process
#' @examples
#' set.seed(1)
#' 
#' N = 200
#' x<-seq(1,10*pi, length.out = N) # confounder
#' y<-rnorm(n = length(x), sin(.5*x), .07*x )
#' d <- data.frame(x=x, y=y)
#' d$x <- as.numeric(scale(d$x))
#' d$y <- as.numeric(scale(d$y))
#' 
#' d_test = data.frame(x=seq(max(d$x), max(d$x+1 ), .01 ))
#' 
#' 
#' res = fDPMix(d_train = d, d_test = d_test, formula = y ~ x,
#'              iter=100, burnin=50, tau_x = c(.01, .001) )
#' 
#' 
#' 
#' par(mfrow=c(1,1))
#' plot(d$x,d$y, pch=20, xlim=c(min(d$x), max(d$x)+1 ), col='gray')
#' 
#' credible_gradient(d$x, res$train)
#' credible_gradient(d_test$x, res$test, col_gradient = 'pink', col_mean_line = 'darkred')
#' 
#' points(d$x, d$y, pch=20, col='gray') ## re-plot points to bring the to the top.
#' abline(v=max(d$x), lty=2)
#' 
#' @export
credible_gradient = function(x, post_draws, 
                             col_gradient='lightblue', 
                             col_median_line='steelblue'){
  
  func_args<-mget(names(formals()),sys.frame(sys.nframe()))
  error_check_credible_gradient(func_args)
  
  colfunc <- colorRampPalette(c("white", col_gradient))
  ci_perc = seq(.99,.01,-.01)
  colvec = colfunc(length(ci_perc))
  names(colvec) = ci_perc
  
  for(i in ci_perc){
    pci = apply(post_draws, 1, quantile, probs=c( (1-i)/2, (1+i)/2  ) )
    polygon(c(x,rev(x)),c(pci[1,],rev(pci[2,])),
            col = colvec[as.character(i)], border = FALSE)
  }
  
  lines(x, apply(post_draws, 1, median), col=col_median_line , pch=20 )
}