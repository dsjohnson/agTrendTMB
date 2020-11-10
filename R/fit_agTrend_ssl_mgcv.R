#' @title Site-level fitting function
#' 
#' @description Fit separate discrete-normal distribution models to an Steller sea lion counts at a single site
#' 
#' @param .x Data frame for an individual site. Must have been passed through \cite{prep_for_fit} first!
#' @param family Distribution family for GAM model (see \code{?family.mgcv} for guidance).
#' @param obl.corr Logical. Should correction be made for oblique photos.
#' @param debug Logical. Step into fitting function for interactive evaluation
#' @param ... Additional arguments passed to \code{mgcv::gam} for fitting
#' 
#' @details There are 3 models which can be used for within site interpolation:
#' * \code{const} This interpolates using a simple mean and variance for the site
#' * \code{lin} A linear fit with time is used for interpolation
#' * \code{gp} i.e., Gaussian Process. This model used a time trend plus a random walk of order 2 (RW(2)) to interpolate
#' counts. Further, the RW(2) process can be further restrained to minimize 'wiggliness' of the curve when necessary. This
#' prevents overfitting of the curve.
#' @md
#' 
#' For all models a 'discrete-normal' model is use for the counts. This prevents the need to use a link function
#' and the data can be modeled on the identity scale. 
#' 
#' @return A list with the following elements:
#' * summary Parameter and abundance estimate information?
#' * raw_data Original data plus oblique photo correction estimates
#' * q_output quality control diagnostic data
#' @md
#' 
#' @import dplyr 
#' @export

fit_agTrend_ssl_mgcv <- function(.x, family="tweedie", obl.corr=FALSE, debug=FALSE,...){
  if(debug) browser()
  abund.name <- attr(.x, "abund.name")
  time.name <- attr(.x, "time.name")
  t1 <- min(.x[[time.name]])
  t2 <- max(.x[[time.name]])
  if(!family%in%c("tweedie", "poisson","negbin")) stop("family needs to be one of: 'tweedie', 'poisson',or 'negbin'")
  if(obl.corr){
    if(!("obl" %in% colnames(.x))){
      warning(" 'obl' column not in data. No photo format correction will be performed!")
      qpts <- 1
      w <- 1
      corr_idx <- rep(0, nrow(.x))
    } else{
      qpts <- exp(0.03903366 + seq(-3, 3, 0.5)*0.01068773)
      w <- pnorm(c(-Inf, seq(-2.75, 2.75, 0.5), Inf)) %>% diff()
      mi_data <- .x %>% filter(obl==0) %>% mutate(w=1)
      mi_alt <- .x %>% filter(obl==1)
      y_alt <- mi_alt[[abund.name]]
      for(i in 1:length(w)){
        mi_alt[[abund.name]] <- y_alt*qpts[i]
        mi_alt$w <- w[i]
        mi_data <- rbind(mi_data, mi_alt)
      }
    }
  } else{
    mi_data <- .x %>% filter(!is.na(.data[[abund.name]]))
    mi_data$w <- 1
  }
  mi_data <- arrange(mi_data, year, .data[[abund.name]])
  mi_data[[abund.name]] <- round(mi_data[[abund.name]])
  corr_idx <- .x$obl==1
  na <- sum(!is.na(.x[[abund.name]]))
  .x <- .x %>% 
    mutate(
      y = as.integer(.data[[abund.name]]),
      y_mi = corr_idx*exp(0.03903366)*y + (1-corr_idx)*y,
      y_mi_low = corr_idx*exp(0.03903366 - 1.96*0.01068773)*y + (1-corr_idx)*y,
      y_mi_up = corr_idx*exp(0.03903366  + 1.96*0.01068773)*y + (1-corr_idx)*y,
      time = .data[[time.name]]-min(.data[[time.name]]),
    )
  if(model=="const"){
    mod <- as.formula(paste0(abund.name, "~ 1"))
    scale <- 1
  } else if(model=="lin"){
    mod <- as.formula(paste0(abund.name, "~ time"))
    scale <- 0
  } else{
    mod <- as.formula(paste0(abund.name, " ~ s(time, k=", na-1, ", bs='cr')" ))
    scale <- 0
  }
  if(family=="poisson") f <- poisson()
  if(family=="tweedie") f <- tw()
  if(family=="negbin") f <- nb()
  
  a <- mi_data$w
  
  fit <- suppressWarnings(mgcv::gam(mod, data=mi_data, family=f, method="REML", scale=scale, weights=a, select=TRUE))
  fit_list = list(fit=fit, mi_data=mi_data, data=.x)
  return(fit_list)
  
}


###############################
### Tweedie functions
################################
#' @title Tweedie distribution properties in fitted model 
#' @description Obtain the Tweedie distribution properties from a fitted model object
#' @param object A model object from \code{mgcv::gam} fitted with \code{family=tw()}
#' @return a list with Tweedie parameters \code{p}, \code{phi}
#' @import mgcv
#' @export
extract_tweedie <- function(object){
  p <- object$family$getTheta(TRUE)
  phi <- summary(object)$dispersion
  prob_nonzero <- function(mu){
    return(1-exp(-mu^(2-p) / (phi*(2-p))))
  }
  comp_pois <- function(mu){
    return(
      list(
        lambda <- mu^(2-p)/(phi*(2-p)),
        alpha <- (2-p)/(p-1),
        gamma <- phi*(p-1)*mu^(p-1)
      )
    )
  }
  return(list(p=p, phi=phi, prob_nonzero=prob_nonzero, comp_pois=comp_pois))
}


