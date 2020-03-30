#' Function for plotting posterior regression line 
#'
#' This function takes as input some continuous covariate \code{x} and predictions from NDPMix, ZDPMix, fDPMix, etc. It outputs a credible band and posterior mean line of the predictions across \code{x}
#' 
#' Please see \url{https://stablemarkets.github.io/ChiRPsite/index.html} for examples and details.
#' 
#' @import stats
#' 
#' @param x A vector of type \code{as.numeric} to plot on the x dimension of the plot.
#' @param post_draws matrix with dimensions length(x) by # posterior draws.
#' @param col_gradient Optional. Defaults to lightblue. Darkest color that gradient goes to.
#' @param col_mean_line Optional. Defaults to steelblue. Color of the posterior mean regression line.
#' @return overlays credible band and posterior mean line.
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
                             col_mean_line='steelblue'){
  
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
  
  lines(x, apply(post_draws, 1, mean), col=col_mean_line , pch=20 )
}