#' @title Add realized SSL survey values to the prediction sample
#' @param x A list object from \code{fit_agtrend_ssl}
#' @param N_sample A matrix of abundance samples from \code{sample_N}
#' @importFrom stats rnorm
#' @export
#' 
get_real <- function(x, N_sample){
  d <- x$raw_data
  use_mi <- nrow(x$q_output$value)>1
  for(i in 1:ncol(N_sample)){
    if(is.na(d$y[i])){
      next
    }else{
      if(d$obl[i]==1 & use_mi){
        N_sample[,i] <- d$y[i]*exp(rnorm(nrow(N_sample), 0.03903366, 0.01068773))
      } else{
        N_sample[,i] <- d$y[i]
      }
    }
  }
  return(N_sample)
}
