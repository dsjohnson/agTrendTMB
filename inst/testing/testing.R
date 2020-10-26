library(tidyverse)
# library(agTrendTMB)
library(foreach)
library(doFuture)
registerDoFuture()

load("~/research/projects/sea_lion_analysis/SSL survey power analysis_2019/ssl_counts_v8.RData")
# source("helper_gam.R")

np <- edpsnp %>% #filter(REGION=="SE AK") %>% 
  arrange(SITE, year) %>% left_join(edps_photo) %>% 
  filter(year>=1979, !is.na(counts)) %>% 
  filter(year<=2015)

## Add group info for aggregations
regions <- select(np, SITE, REGION) %>% distinct()
np <- prep_for_fit(np) %>% left_join(regions)
np$Total <- factor("Total")

## Make a 'safe' version of the fitting function
safe_agTrend_ssl <- fit_agTrend_ssl %>% safely()


plan("multisession", workers=6)

np[["fit"]] <- foreach(i=1:nrow(np))%dopar%{
  .x <- np$data[[i]]
  model <- np$model[[i]]
  # if(model!="gp") return(NA)
  ## Fix variance estimate if counts are all equal
  if(var(.x$counts, na.rm=T)==0){
    map.override <- list(xi=factor(c(NA,NA)))
  } else{
    map.override <- NULL
  }
  ###
  output <- safe_agTrend_ssl(.x, model, map.override=map.override, obl.corr=TRUE)
  if(!is.null(output$result)){
    output <- output$result
  } else{
    output <- output$error
  }
  return(output)
}

plan("sequential")

### Make site-level plots
np <- np %>% mutate(
  plots = map2(fit, SITE, plot_fit, plot=FALSE)
)

## Create N samples
np <- np %>% mutate(
  N_predict = map(fit, sample_N, size=5000),
  N_real = map2(fit, N_predict, get_real)
)


## Get aggegated abundance
np_total <- agg_abund(np, "Total")
np_region <- agg_abund(np, "REGION")

### Summary for plotting
N_total_df <- summary_agg(np_total) %>% filter(year>=1985)
N_region_df <- summary_agg(np_region) %>% filter(year>=1985)

p1 <- ggplot(N_total_df) + geom_path(aes(y=Estimate, x=year)) +
  geom_ribbon(aes(x=year, ymax=CI_predict_upper, ymin=CI_predict_lower), alpha=0.2, fill="blue") +
  geom_pointrange(aes(x = year, y=Estimate_real, ymax=CI_real_lower, ymin=CI_real_upper),
                  data=N_total_df) +
  ylab("Count") + xlab("Year") + theme_bw() 
print(p1)

p2 <- ggplot(N_region_df) + geom_path(aes(y=Estimate, x=year)) +
  geom_ribbon(aes(x=year, ymax=CI_predict_upper, ymin=CI_predict_lower), alpha=0.2, fill="blue") +
  geom_pointrange(aes(x = year, y=Estimate_real, ymax=CI_real_lower, ymin=CI_real_upper),
                  data=N_region_df %>% filter(survey==1)) +
  facet_wrap(~REGION, nrow=3, scales="free_y") + theme_bw()
print(p2)
