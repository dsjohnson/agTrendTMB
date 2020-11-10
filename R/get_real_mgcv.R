#' @title Add realized SSL survey values to the prediction sample
#' @param x A list object from \code{fit_agtrend_ssl}
#' @param N_sample A matrix of abundance samples from \code{sample_N}
#' @importFrom stats rnorm
#' @export
#' 
get_real_mgcv <- function(fit, sample){
  use_mi <- any(fit$mi_data$w<1)
  d <- fit$data
  for(i in 1:ncol(sample)){
    if(is.na(d$y[i])){
      next
    }else{
      if(d$obl[i]==1 & use_mi){
        sample[,i] <- d$y[i]*exp(rnorm(nrow(sample), 0.03903366, 0.01068773))
      } else{
        sample[,i] <- d$y[i]
      }
    }
  }
  return(sample)
}
