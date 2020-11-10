#' @title Sample abundance from fitted agTrendTMB model
#' @param x A list object resulting from \code{fit_agTrend_ssl}
#' @param size Sample size. Defaults to 1000
#' 
#' @import mvnfast truncnorm
#' @export
#' 
sample_N <- function(x, size=1000){
  if(is.null(x)) return(NULL)
  res <- x$summary
  S <- res$cov
  m <- res$value
  o1 <- rev(order(diag(S)))
  revo1 <- order(o1)
  m <- m[o1]
  S <- S[o1,o1]
  eee <- eigen(S)
  V <- eee$vectors%*%diag(sqrt(round(eee$values, 7)))
  spl <- mvnfast::rmvn(size, m, sigma=V, isChol=TRUE) 
  spl <- spl[,revo1]
  tmax <- length(m)-3
  mu_spl <- as.vector(spl[,1:tmax])
  sd_N_spl <- exp(rep(spl[,tmax+1],tmax) + rep(spl[,(tmax+2)],tmax)*mu_spl)
  N_spl <- truncnorm::rtruncnorm(length(mu_spl), mean=mu_spl, sd=sd_N_spl) %>% pmax(.,0) %>%
    round() %>% matrix(.,nrow=size)
  d <- x$raw_data
  tn <- attr(d, "time.name")
  attr(N_spl, "time.name") <- tn
  attr(N_spl, tn) <- d[[tn]]
  attr(N_spl, "abund.name") <- attr(d, "abund.name")
  attr(N_spl, "survey_times") <- attr(d, "survey_times")
  return(N_spl)
}
