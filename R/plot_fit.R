#' @title Plot A Site-Level Fit For Diagnostic Checking
#' @param fit_list A list resulting from a call to \code{fit_agTrend_ssl}
#' @param title plot title
#' @param plot Logical. Whether or not to produce the plot or just return the \code{ggplot2} object for further modification
#' 
#' @return A \code{ggplot2} object
#' 
#' @import ggplot2
#' @export
#' 
plot_fit <- function(fit_list, title=NULL, plot=TRUE){
  r <- fit_list$summary
  xi <- r$value[grep("xi",names(r$value))]
  d <- fit_list$raw_data 
  tn <- attr(d, "time.name")
  an <- attr(d, "abund.name")
  d$ttt <- d[[tn]]
  d$mu <- r$value[grep("mu",names(r$value))]
  d$mu_sd <- r$sd[grep("mu",names(r$value))]
  d$sigma <- exp(xi[1] + xi[2]*d$mu)
  d$pz <- pnorm(0.5, d$mu, d$sigma)
  
  p <- ggplot()
  probs <- seq(0.1, 0.9, 0.1)
  for(i in 1:length(probs)){
    mu_low=pmax(0,d$mu - qnorm(1-(1-probs[i])/2)*(d$mu_sd + d$sigma))
    mu_up=pmax(0,d$mu + qnorm(1-(1-probs[i])/2)*(d$mu_sd + d$sigma))
    ud <- data.frame(ttt=d$ttt, mu_low, mu_up)
    p <- p + geom_ribbon(aes(x=ttt, ymin=mu_low, ymax=mu_up), data=ud, fill="cadetblue3", alpha=1/8)
  }
  p <- p + geom_point(aes(x=ttt, y=0, alpha=pz), shape=15, size=3, color="cadetblue3", data=d) + 
    scale_alpha_identity(name="Prob. zero obs.") +
    geom_path(aes(x=ttt, y=pmax(0, mu)), color="cadetblue4", alpha=0.2, data=d, lwd=1.5)
  p <- p+ geom_pointrange(
    aes(x=ttt, y=y_mi, ymax=y_mi_up, ymin=y_mi_low), 
    data=d %>% filter(!is.na(y)), color="gray33") + 
    # geom_point(aes(x=time, y=y_mi), data=d %>% filter(!is.na(y)), color="red", pch="-", size=10) +
    geom_point(aes(x=ttt, y=y), data=d %>% filter(!is.na(y)), shape=1) +
    theme_classic() + xlab(toupper(tn)) + ylab(toupper(an))
  if(!is.null(title)) p <- p + ggtitle(title)
  
  if(plot) print(p)
  return(p)
}
