#' @title Site-level fitting function
#' 
#' @description Fit separate discrete-normal distribution models to an Steller sea lion counts at a single site
#' 
#' @param .x Data frame for an individual site. Must have been passed through \cite{prep_for_fit} first!
#' @param model Model for interpolation: 'const', 'lin', or 'gp' (see Details)
#' @param penalty Wiggliness penalty for 'gp' model
#' @param obl.corr Logical. Should correction be made for oblique photos.
#' @param map.override List. method for overriding parameter parameters set to fixed values. This should rarely need to be used
#' @param silent Logical. Run without messages
#' @param raw.return Logical. Return raw fitting results. This is useful for diagnostic putposes when the models fails to fit
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
#' @import TMB dplyr purrr 
#' @export

fit_agTrend_ssl <- function(
  .x, model, obl.corr=FALSE, penalty=FALSE,map.override=NULL, 
  silent=TRUE, raw.return=FALSE
){
  # browser()
  abund.name <- attr(.x, "abund.name")
  time.name <- attr(.x, "time.name")
  t1 <- min(.x[[time.name]])
  t2 <- max(.x[[time.name]])
  
  if(obl.corr){
    if(!("obl" %in% colnames(.x))){
      warning(" 'obl' column not in data. No photo format correction will be performed!")
      qpts <- 1
      w <- 1
      corr_idx <- rep(0, nrow(.x))
    } else{
      qpts <- exp(0.03903366 + seq(-3, 3, 0.5)*0.01068773)
      w <- pnorm(c(-Inf, seq(-2.75, 2.75, 0.5), Inf)) %>% diff()
      corr_idx <- .x$obl==1
    }
  } else{
    qpts <- 1
    w <- 1
    corr_idx <- rep(0, nrow(.x))
  }
  .x <- .x %>% 
    mutate(
      y = as.integer(.data[[abund.name]]),
      y_mi = corr_idx*exp(0.03903366)*y + (1-corr_idx)*y,
      y_mi_low = corr_idx*exp(0.03903366 - 1.96*0.01068773)*y + (1-corr_idx)*y,
      y_mi_up = corr_idx*exp(0.03903366  + 1.96*0.01068773)*y + (1-corr_idx)*y,
      time = .data[[time.name]]-min(.data[[time.name]]),
    )
  
  fit <- opt <- vector("list", length(qpts))
  
  if(model=="const"){
    X <- model.matrix(~1, data=.x)
    K <- model.matrix(~1, data=.x)
  }
  if(model=="lin"){
    X <- model.matrix(~time, data=.x)
    K <- model.matrix(~1, data=.x)
  }
  if(model=="gp"){
    X <- model.matrix(~time, data=.x)
    P <- X%*%solve(crossprod(X))%*%t(X)
    df <- min(10,floor((sum(!is.na(.x[[abund.name]]))-3)/2))
    knots <- seq(t1, t2, length=df)
    h <- max(diff(knots))
    K <- outer(.x[[time.name]], knots, function(x,y){dnorm(x, y, h)})
    # df <- nrow(.x)-3
    # K <- iar_basis(nrow(.x), 2, df)
    K <- K-P%*%K
    Kaug <- rbind(K, 1*diag(ncol(K)))
  }
  
  ### TMB data
  tmb_data <- list(
    model="disc_norm_gp",
    y = .x$y,
    X = X,
    K = K,
    lambda_tau = -log(1.0e-6)*3,
    lambda_xi = 0.5,
    lambda_beta_0 = max(.x$y, na.rm=TRUE),
    penalty=as.integer(penalty)
  )
  
 
  ### TMB map list
  tmb_map <- list(xi = factor(c(1,NA)))
  if(model!="gp"){
    tmb_map <- c(tmb_map,
                 list(
                   alpha=factor(NA),
                   ln_tau=factor(NA)
                 )
    )
  }
  if(!is.null(map.override)){
    or_fix <- names(map.override)[names(map.override)%in%names(tmb_map)]
    new_fix <- names(map.override)[!names(map.override)%in%names(tmb_map)]
    tmb_map[or_fix] <- map.override[or_fix]
    tmb_map <- c(tmb_map, map.override[new_fix])
  }
  
  # Random desigation
  random <- "alpha"
  
  # Begin MI loop
  safe_optim <- safely(optim)
  safe_sdreport <- safely(TMB::sdreport)
  
  for(j in 1:length(qpts)){
    # browser()
    ### TMB data
    tmb_data$y <- .x$y*(1-corr_idx) + qpts[j]*.x$y*(corr_idx)
    if(model=="const"){
      ff <- lm(tmb_data$y~0+X)
      ln_sigma <- log(summary(ff)$sigma)
      ln_tau <- 0
      alpha <- 0
      beta <- coef(ff)
      xi <- c(log(max(1/6, summary(ff)$sigma)), 0)
    } else if(model=="lin"){
      ff <- lm(tmb_data$y~0+X)
      ln_sigma <- 2+ log(summary(ff)$sigma)
      ln_tau <- 0
      alpha <- 0
      beta <- coef(ff)
      xi <- c(log(max(1/6, summary(ff)$sigma)), 0)
    } else{
      yaug <- c(tmb_data$y, rep(0,ncol(K)))
      Xaug <- rbind(X,matrix(0,ncol(K),2))
      ff <- lm(yaug ~ 0 + Xaug + Kaug)
      ln_sigma <- 2+ log(summary(ff)$sigma)
      ln_tau <- 1
      alpha <- rep(0,ncol(K))
      xi <- c(log(max(1/6, summary(ff)$sigma)), 0)
      alpha <- coef(ff)[-c(1:2)]
      ln_tau <- log(sd(alpha))
      beta <- coef(ff)[1:2]
    }
    
    # if(j==13) browser()
    
    ### TMB parameters
    tmb_par <- list(
      beta=beta,
      alpha=alpha,
      ln_tau=ln_tau,
      xi=xi
    )
    
    foo <- TMB::MakeADFun(
      tmb_data, tmb_par, 
      random=random,
      map = tmb_map, DLL="agTrendTMB_TMBExports",
      silent = silent)
    
    opt[[j]] <- safe_optim(
      foo$par, foo$fn, foo$gr, method="BFGS",
      control = list(maxit=1000, fnscale=1)
    )
    
    # browser()
    
    fit[[j]] <- safe_sdreport(foo, getJointPrecision = TRUE)
    
    # if(!fit[[j]]$result$pdHess & model=="gp"){
    #   # browser()
    #   # tmb_par <- map(tmb_par, ~{.x+rnorm(length(.x), 0, 0.1)})
    #   ff <- lm(tmb_data$y ~ 0 + X)
    #   tmb_par$alpha <- 0.0*tmb_par$alpha 
    #   tmb_par$beta <- coef(ff)[1:2]
    #   tmb_par$xi <- c(log(max(1/6, summary(ff)$sigma)), 0)
    #   foo <- TMB::MakeADFun(
    #     tmb_data, tmb_par, 
    #     random=random,
    #     map = tmb_map, DLL="agTrendTMB_TMBExports",
    #     silent = silent)
    #   opt[[j]] <- safe_optim(
    #     foo$par, foo$fn, foo$gr, method="BFGS",
    #     control = list(maxit=1000, fnscale=2)
    #   )
    #   fit[[j]] <- safe_sdreport(foo, getJointPrecision = TRUE)
    # }
  }  
  
  # browser()
  
  opt_check <- map_lgl(map(opt, ~{.x$error}), is.null)
  fit_check <- map_lgl(map(fit, ~{.x$error}), is.null)
  pdHess_check <- TRUE #map_lgl(fit, ~{.x$result$pdHess}) %>% ifelse(is.na(.), FALSE, .)
  full_check <- opt_check*fit_check*pdHess_check
  ngood <- sum(full_check)
  if(raw.return){
    return(list(opt=opt, fit=fit, opt_check=opt_check, fit_check=fit_check, pdHess_check=pdHess_check, qpts=qpts, w=w))
  }
  if(ngood!=length(qpts)) stop("There were some bad fits! set 'raw.return=TRUE' and check results.")
  
  # browser()
  
  fits <- map(fit, ~{.x$result})
  value_mat <- map(fits, ~{.$value}) %>% do.call("rbind",.)
  cov_list <- map(fits, ~{.$cov}) 
  par.fixed_mat <- map(fits, ~{.$par.fixed}) %>% do.call("rbind",.)
  cov.fixed_list <- map(fits, ~{.$cov.fixed})
  par.random_mat <- map(fits, ~{.$par.random}) %>% do.call("rbind",.)
  diag.cov.random_mat <- map(fits, ~{.$diag.cov.random}) %>% do.call("rbind",.)
  jointPrecision <- map(fits, ~{.$jointPrecision})
  summary <- list(
    value = colSums(value_mat*w),
    cov = Reduce(`+`,Map(`*`, cov_list, w)) + cov.wt(value_mat, w, method="ML")$cov,
    par.fixed = colSums(par.fixed_mat*w),
    cov.fixed = Reduce(`+`,Map(`*`, cov.fixed_list, w)) + cov.wt(par.fixed_mat, w, method="ML")$cov
  )
  summary$sd <- sqrt(diag(summary$cov))
  if(!is.null(par.random_mat)){
    summary <- c(summary,
                 list(
                   par.random = colSums(par.random_mat*w),
                   diag.cov.random = colSums(diag.cov.random_mat*w) + diag(cov.wt(par.random_mat, w, method="ML")$cov)
                 )
    )
  } 
  fit_list <- list(
    summary=summary, 
    raw_data = .x,
    q_output=list(value=value_mat, cov=cov_list, par.fixed=par.fixed_mat, 
                  cov.fixed=cov.fixed_list, par.random=par.random_mat,
                  diag.cov.random=diag.cov.random_mat, 
                  jointPrecision=jointPrecision,
                  opt_check=opt_check, fit_check=fit_check, 
                  weights=w,
                  pdHess_check=pdHess_check
    )
  )
  return(fit_list)
}
