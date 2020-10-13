#' @title Trim sites that do not fit the criteria for model fitting
#' @param x Steller sea lion count data for one DPS and age class, e.g., eDPS nonpups
#' @param site.name Name of column with site names
#' @param abund.name Name of the column with the abundance counts
#' @param time.name Name of the column with the time (e.g., year) data
#' @param model.cuts Number of nonzero surveys for each model type (see Details)
#' @param non.zero Minimum number of nonzero counts for inclusion in model fitting. Must be >2.
#' @param total Minimum total number of sea lions observed at the site over the time span considered
#' @details The \code{model.cuts} is a 2-vector argument that provides the boundaries for each model type, 'const' (constant meand term),
#' 'lin' (linear slope mean term) or 'gp' (Gaussian process mean term). E.g., the default value is
#' \code{model.cuts=c(3,10)} which means that any site with less than or equal to 3 nonzero 
#' survey values with be fit with a 'const' model and if the number of nonzero counts is in (3,10] a
#' 'lin' model is fitted, and so on.
#' @import dplyr purrr
#' @export
prep_for_fit <- function(x, site.name="SITE",abund.name="counts", 
                         time.name="year", model.cuts=c(3,10), non.zero=2, total=10){
  
  time_df <- data.frame(min(x[[time.name]]):max(x[[time.name]])) %>% `colnames<-`(time.name)
  
  x <- x %>% group_by(.data[[site.name]]) %>% nest() %>% ungroup() %>% 
    mutate(
      data = map(data, 
                 ~{
                   full_join(.x, time_df, by=time.name) %>% 
                     arrange(.data[[time.name]]) %>% 
                     mutate(
                       time = .data[[time.name]]-min(.data[[time.name]])#,
                       #obl = ifelse(is.na(obl), 0, obl)
                       )
                 }
      )
    ) 
  
  x <- x %>% mutate(
    num_nonzero = map_int(data, ~{sum(.x[[abund.name]]>0, na.rm=T)}),
    total_abund = map_dbl(data, ~{sum(.x[[abund.name]], na.rm=T)})
    )
  
  model.cuts <- c(0, model.cuts, nrow(x$data[[1]]))

  x <- x %>% filter(num_nonzero>={{non.zero}}, total_abund>={{total}}) %>%
    mutate(
      model = cut(num_nonzero, model.cuts, labels=c("const","lin","gp")) %>% 
        as.character()
    ) %>% select(-num_nonzero, -total_abund)
  
  x <- x %>% mutate(
    data = map(data, ~{
      surv <- filter(.x, !is.na(.x[[abund.name]]))
      surv <- surv[[time.name]] %>% unique() %>% sort()
      attr(.x, "survey_times") <- surv
      attr(.x, "abund.name") <- abund.name
      attr(.x, "time.name") <- time.name
      .x
    })
  )
  
  return(x)
  
}