#' @title Aggregate abundances based on a factor variable
#' @param x A data set containing columns '\code{agg.var}', '\code{sample}' 
#' and '\code{realized}'. The last two are obtained from a calls to 
#' \code{agTrendTMB::sample_N} and \code{agTrendTMB::get_real}.
#' @param agg.var Variable used for aggregation
#' @import purrr dplyr
#' @export
#' 
agg_abund_mgcv <- function(x, agg.var){
  x$surv_times <- map(x$sample, ~{attr(.x, "survey_times")})
  results <- x %>% group_by(.data[[agg.var]]) %>% nest() %>%
    mutate(
      surv_times = map(data, ~{sort(reduce(.x$surv_times, union))}),
      sample = map(data, ~{reduce(.x$sample, `+`)})
    ) %>% mutate(
      sample = map2(sample, surv_times, ~{
        attr(.x, "survey_times") <- .y
        .x
      })
    )
  if("realized"%in%colnames(x)){
    results <- results %>% 
      mutate(realized = map(data, ~ reduce(.x$realized, `+`))) %>% 
      mutate(
        realized = map2(realized, surv_times, ~{
          attr(.x, "survey_times") <- .y
          .x
        })
      )
  }
  results %>% ungroup() %>% select(-data, -surv_times)
}


#' @title Summarize aggregation results
#' @param x An aggregated abundance object from \code{agg_abund} 
#' @param ci.prob Confidence (credible) interval value
#' @importFrom coda HPDinterval mcmc
#' @export
#' 
summary_agg_mgcv <- function(x, ci.prob=0.9){
  tn <- attr(x$sample[[1]], "time.name")
  x[[tn]] <- map(x$sample, ~{attr(.x, "survey_times")})
  surv_times <- x %>% select(1, .data[[tn]]) %>% unnest(cols=.data[[tn]]) %>% 
    mutate(survey=1)
  x <- select(x, -.data[[tn]])
  df <- data.frame(attr(x$sample[[1]], tn)) %>% `colnames<-`(tn)
  x$Estimate <- map(x$sample, ~{
    out <- data.frame(Estimate=apply(.x, 2, median))
    bind_cols(df, out)
  })
  x$CI <- map(x$sample, ~{
    coda::HPDinterval(coda::mcmc(.x), prob=ci.prob) %>% as.data.frame() %>% 
      `colnames<-`(c("CI_predict_lower", "CI_predict_upper"))
  })
  x <- x[,-c(which(colnames(x)=="sample"))] 
  if("realized" %in% colnames(x)){
    x$Estimate_real <- map(x$realized, ~{data.frame(Estimate_real=apply(.x, 2, median))})
    x$CI_real <- map(x$realized, ~{
      coda::HPDinterval(coda::mcmc(.x), prob=ci.prob) %>% as.data.frame() %>% 
        `colnames<-`(c("CI_real_lower", "CI_real_upper"))
    })
    x <- x[,-c(which(colnames(x)=="realized"))] 
    x <- x %>% unnest(cols = c(Estimate, CI, Estimate_real, CI_real))
    
  } else{
    x <- x %>% unnest(cols = c(Estimate, CI))
  }
  x <- x %>% left_join(surv_times, by=colnames(surv_times)[1:2]) %>% 
    mutate(
      survey = ifelse(is.na(survey), 0, survey)
    )
  return(x)
}
