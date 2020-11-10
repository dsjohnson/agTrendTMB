#' @title Sample abundance from fitted agTrendTMB model
#' @param x A list object resulting from \code{fit_agTrend_ssl}
#' @param size Sample size. Defaults to 1000
#' @param summary Logical. Whether or not a site-level summary while sampling is occuring.
#' @param
#' 
#' @importFrom mvnfast rmvn 
#' @importFrom mgcv rTweedie
#' @importFrom coda HPDinterval mcmc
#' @export
#' 
sample_N_mgcv <- function(x, size=1000, summary=TRUE, block=1000){
  f <- x$fit$family$family
  
  # abundance sampler
  if(grepl("Tweedie",f)){
    phi <- summary(x$fit)$dispersion
    p <- x$fit$family$getTheta(TRUE)
    sampler <- function(mu){
      mgcv::rTweedie(mu, p, phi)
    }
  }
  if(grepl("Negative Binomial",f)){
    theta <- x$fit$family$getTheta(TRUE)
    sampler <- function(mu){
      rnbinom(length(mu), mu=mu, size=theta)
    }
  }
  if(grepl("poisson",f)){
    sampler <- function(mu){
      rpois(length(mu), lambda=mu)
    }
  }
  
  V <- mgcv::vcov.gam(x$fit, unconditional=TRUE)
  X <- mgcv::predict.gam(x$fit, newdata=x$data, type='lpmatrix')
  m <- coef(x$fit)
  ###
  ###
  N_smp <- NULL
  n_blocks <- floor(size/block)
  rem <- size%%block
  if(n_blocks>0){
    for(j in 1:n_blocks){
      mu_smp <- mvnfast::rmvn(block, m, V)%*%t(X) %>% as.vector() %>% exp()
      tmp <- sampler(mu_smp) %>% matrix(nrow=block)
      N_smp <- rbind(N_smp, tmp)
    }
  }
  if(rem > 0){
    mu_smp <- mvnfast::rmvn(rem, m, V)%*%t(X) %>% as.vector() %>% exp()
    tmp <- sampler(mu_smp) %>% matrix(nrow=rem)
    N_smp <- rbind(N_smp, tmp) 
  }
  ###
  ###
  tn <- attr(x$data, "time.name")
  attr(N_smp, "time.name") <- tn
  attr(N_smp, tn) <- x$data[[tn]]
  attr(N_smp, "abund.name") <- attr(x$data, "abund.name")
  attr(N_smp, "survey_times") <- attr(x$data, "survey_times")
  if(!summary) return(N_smp)
  x$data$Estimate <- exp(X%*%m) %>% round() %>% as.numeric()
  x$data$SE <- apply(N_smp, 2, sd) %>% round()
  x$data <- bind_cols(x$data, data.frame(coda::HPDinterval(coda::mcmc(N_smp))))
  return(tibble(sample=list(N_smp), summary=list(x$data)))
}
