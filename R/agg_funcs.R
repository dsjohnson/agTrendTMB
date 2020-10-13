#' @title Aggregate abundances based on a factor variable
#' @param x A data set containing columns '\code{agg.var}', '\code{N_predict}' 
#' and '\code{N_real}'. The last two are obtained from a calls to 
#' \code{agTrendTMB::sample_N} and \code{agTrendTMB::get_real}.
#' @param agg.var Variable used for aggregation
#' @import purrr dplyr
#' @export
#' 
agg_abund <- function(x, agg.var){
  results <- x %>% group_by(.data[[agg.var]]) %>% nest() %>%
    mutate(N_predict = map(data, ~ reduce(.$N_predict, `+`)))
  if("N_real"%in%colnames(x)){
    results <- results %>% mutate(N_real = map(data, ~ reduce(.$N_real, `+`)))
  }
  results %>% ungroup() %>% select(-data)
}


#' @title Summarize aggregation results
#' @param x An aggregated abundance object from \code{agg_abund} 
#' @param ci.prob Confidence (credible) interval value
#' @importFrom coda HPDinterval mcmc
#' @export
#' 
summary_agg <- function(x, ci.prob=0.95){
  tn <- attr(x$N_predict[[1]], "time.name")
  df <- data.frame(attr(x$N_predict[[1]], tn)) %>% `colnames<-`(tn)
  x$Estimate <- map(x$N_predict, ~{
    out <- data.frame(Estimate=apply(.x, 2, median))
    bind_cols(df, out)
    })
  x$CI <- map(x$N_predict, ~{
    coda::HPDinterval(coda::mcmc(.x), prob=ci.prob) %>% as.data.frame() %>% 
      `colnames<-`(c("CI_predict_lower", "CI_predict_upper"))
    })
  x <- x[,-c(which(colnames(x)=="N_predict"))] 
  if("N_real" %in% colnames(x)){
    x$Estimate_real <- map(x$N_real, ~{data.frame(Estimate_real=apply(.x, 2, median))})
    x$CI_real <- map(x$N_real, ~{
      coda::HPDinterval(coda::mcmc(.x), prob=ci.prob) %>% as.data.frame() %>% 
        `colnames<-`(c("CI_real_lower", "CI_real_upper"))
      })
    x <- x[,-c(which(colnames(x)=="N_real"))] 
    return(unnest(x, cols = c(Estimate, CI, Estimate_real, CI_real)))
  } else{
    return(unnest(x, cols = c(Estimate, CI)))
  }
  
}
  