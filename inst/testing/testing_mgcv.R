library(tidyverse)
# library(agTrendTMB)
devtools::load_all(".")
library(foreach)
library(doFuture)
registerDoFuture()

load("~/research/projects/sea_lion_analysis/SSL survey power analysis_2019/ssl_counts_v8.RData")
# source("helper_gam.R")

np <- edpsnp %>% #filter(REGION=="SE AK") %>% 
  arrange(SITE, year) %>% left_join(edps_photo) %>% 
  filter(year>=1979, !is.na(counts)) %>% 
  filter(year<=2019)

## Add group info for aggregations
regions <- select(np, SITE, REGION) %>% distinct()
np <- prep_for_fit(np) %>% left_join(regions)
np$Total <- factor("Total")

## Make a 'safe' version of the fitting function
# safe_agTrend_ssl <- fit_agTrend_ssl %>% safely()

tictoc::tic()
np[["fit"]] <- foreach(i=c(1:nrow(np))) %do%{
  .x <- np$data[[i]]
  model <- np$model[[i]]
  
  ###
  print(i)
  x <- fit_agTrend_ssl_mgcv(.x, family='tweedie', obl.corr=TRUE)
  return(x)
}
tictoc::toc()

### Make site-level plots
# np <- np %>% mutate(
#   plots = map2(fit, SITE, plot_fit, plot=FALSE)
# )
# map(np$plots, print)

## Create N samples

tictoc::tic()
plan(multisession, workers=6)
# np$N_sample <- foreach(i=c(1:nrow(np))) %dopar% {
foreach(i=1:12) %dopar% {
  out <- sample_N_mgcv(np$fit[[i]], size=3000)
  return(out)
} 
plan(sequential)
tictoc::toc()

np <- unnest(np, cols=N_sample)

np$realized <- foreach(i=c(1:nrow(np))) %do% {
  get_real_mgcv(np$fit[[i]], np$sample[[i]])
}


## Get aggegated abundance
np_total <- agg_abund_mgcv(np, "Total")
np_region <- agg_abund_mgcv(np, "REGION")

### Summary for plotting
N_total_df <- summary_agg_mgcv(np_total) %>% filter(year>=1989)
N_region_df <- summary_agg_mgcv(np_region) %>% filter(year>=1989)

p1 <- ggplot(N_total_df) + geom_path(aes(y=Estimate, x=year)) +
  geom_ribbon(aes(x=year, ymax=CI_predict_upper, ymin=CI_predict_lower), alpha=0.2, fill="blue") +
  geom_pointrange(aes(x = year, y=Estimate_real, ymax=CI_real_lower, ymin=CI_real_upper),
                  data=N_total_df %>% filter(survey==1)) +
  ylab("Count") + xlab("Year") + theme_bw() 
print(p1)

p2 <- ggplot(N_region_df) + geom_path(aes(y=Estimate, x=year)) +
  geom_ribbon(aes(x=year, ymax=CI_predict_upper, ymin=CI_predict_lower), alpha=0.2, fill="blue") +
  geom_pointrange(aes(x = year, y=Estimate_real, ymax=CI_real_lower, ymin=CI_real_upper),
                  data=N_region_df %>% filter(survey==1)) +
  facet_wrap(~REGION, nrow=3, scales="free_y") + theme_bw()
print(p2)
