#### Setup ####
rm(list=ls())

library(plm, quietly=TRUE)
library(tidyverse, quietly=TRUE)
library(glue)
library(gt)
library(lmtest)

# Import standard variables
source("config.R")

# Each day, the data should be recleaned
# source("clean_full.R")

df_long = readr::read_csv(path_full_long, col_types = do.call(
  cols, list(Date = col_date(format = "%Y-%m-%d"))))

regressors = c("airPassengersArrived", "touristArrivals", "broadbandAccess",
               "dischargeRateDiabetes", "dischargeRateRespiratory",
               "dischargeRateHypertension", "dischargeRateCancer",
               "dischargeRateChd", "dischargeRatePneumonia", "dischargeRateTB",
               "availableBeds", "maritimePassengersDisembarked",
               "riskOfPovertyOrSocialExclusion", "railTravelers", "medianAge")

#### Transform variables into proportions ####
# TODO: Transform variables to proportions (in clean_full.R or here?)
make_prop = function(x, na.rm = FALSE) { x / sum(x, na.rm = na.rm) }

regressors_prop = c("airPassengersArrived", "touristArrivals",
                    "dischargeRateDiabetes", "dischargeRateRespiratory",
                    "dischargeRateHypertension", "dischargeRateCancer",
                    "dischargeRateChd", "dischargeRatePneumonia",
                    "dischargeRateTB", "availableBeds",
                    "maritimePassengersDisembarked", "railTravelers")

df_long = df_long %>% 
  group_by(Date) %>% 
  mutate_at(regressors_prop, make_prop) %>% 
  ungroup

# TODO: Run models
# Create formula to be the product of the regressors with the lagged incidence
# and susceptible rates

# TODO: Add week/weekend effect
df_long = df_long %>%
  mutate(weekNumber = lubridate::week(df_long$Date)) %>%
  mutate(weekend = lubridate::wday(df_long$Date, label = TRUE)
         %in% c("Sat", "Sun") %>% as.integer %>% as.factor)

X_regressors = c("weekend", "weekNumber")
base_vars = c("Date", "Code", "incidenceRate", "susceptibleRate", X_regressors)

df_long_pd = pdata.frame(df_long %>%
                           select(all_of(c(base_vars, regressors))) %>%
                           arrange(Code),
                         index=c("Code","Date"), drop.index=TRUE,
                         row.names=TRUE)

lags = 1:6 # Incubation period
lsdv_results = vector("list")

resids_col = vector()
lag_col = vector()
dw_col = vector()
dw_p_col = vector()
t_col = vector()
t_p_col = vector()

for (lag in lags){
  print("--------------------")
  print(glue("Running models for incubation period {lag}"))
  
  # Construct formulae
  base_fm = paste("incidenceRate ~",
                  
                  # The eurostat regressors are multiplied by the lagged Inc and S
                  paste0(glue("dplyr::lag(incidenceRate, {lag})", 
                              ":dplyr::lag(susceptibleRate, {lag}):")) %>%
                    paste0(glue("dplyr::lag({regressors}, {lag})")) %>%
                    paste(collapse="+"), "+", 
                  
                  # These include the weekend and weekNumber effect
                  paste(X_regressors, collapse="+"))
  
  fm = as.formula(base_fm)
  fm_lsdv = base_fm %>%
    paste("+factor(Code)") %>%
    as.formula()
  
  # Run models
  # plmo_pooled = plm(fm, data = df_long_pd, model = "pooling")
  # plmo_fe = plm(fm, data = df_long_pd, model = "within")
  
  lsdv = lm(fm_lsdv, data=df_long)
  lsdv_results[[as.character(lag)]] = lsdv
  dw_results = dwtest(fm_lsdv, data=df_long)
  
  residual = df_long$incidenceRate[-1:-lag] - lsdv$fitted.values
  
  resids_col = c(resids_col, residual)
  lag_col = c(lag_col, rep(lag, length(residual)))
  dw_col = c(dw_col, rep(dw_results$statistic, length(residual)))
  dw_p_col = c(dw_p_col, rep(dw_results$p.value, length(residual)))
  
  residuals = c(residual, rep(NA, lag))
  t_results = t.test(residuals)
  
  t_col = c(t_col, rep(dw_results$statistic, length(residual)))
  t_p_col = c(t_p_col, rep(t_results$p.value, length(residual)))
}

# Test
lag_col_counts = table(lag_col)
index = vector()
for (len in lag_col_counts){
  index = c(index, 1:len)
}

tbl = tibble(index, residuals = resids_col, lag=lag_col, dw_col, dw_p_col,
             t_col, t_p_col)
tbl$lag = factor(tbl$lag, levels=unique(tbl$lag),
                 labels=paste0("Lag: ", lags,
                               "\n Durbin-Watson statistic: ",
                               signif(unique(tbl$dw_col), 5),
                               " (p=", signif(unique(tbl$dw_p_col), 5), ")",
                               "\n t-statistic for mean 0: ",
                               signif(unique(tbl$t_col), 5),
                               " (p=", signif(unique(tbl$t_p_col), 5), ")"))

ggplot(tbl, aes(x=index, y=residuals)) +
  facet_wrap(~lag) +
  theme(plot.title = element_text(face = "bold")) +
  geom_point(alpha=0.6, color='firebrick')
ggsave("residuals_plot_lag_1_6.png", path=output_path)


# Create HTML table of LSDV results
rownames_tbl = c("(Intercept)", "weekend1", "weekNumber", "BAS",
                 "BZ", "CAL", "CAM", "EMR", "FVG", "LAZ", "LIG",
                 "LOM", "MAR", "MOL", "PIE", "PUG", "SAR", "SIC",
                 "TN", "TOS", "UMB", "VDA", "VEN",
                 "airPassengersArrived", "touristArrivals",
                 "broadbandAccess", "dischargeRateDiabetes",
                 "dischargeRateRespiratory",
                 "dischargeRateHypertension", "dischargeRateCancer",
                 "dischargeRateChd", "dischargeRatePneumonia",
                 "dischargeRateTB", "availableBeds",
                 "maritimePassengersDisembarked",
                 "riskOfPovertyOrSocialExclusion", "railTravelers", "medianAge")

coefs_tbl = tibble(variable=rownames_tbl)
lags = vector()

for (item in lsdv_results){
  stars = vector()
  
  lag = item$coefficients %>% names() %>% tail(1) %>% str_extract("\\d{1,2}")
  lags = c(lags, lag)
  
  coefs = item$coefficients
  names(coefs) = NULL
  
  tvals = summary(item)$coefficients[, "t value"]
  names(tvals) = NULL
  
  pvals = summary(item)$coefficients[, 4]
  names(pvals) = NULL
  
  for (pval in pvals){
    if (pval < 0.001){
      stars = c(stars, "***")
    } else if (pval < 0.01) {
      stars = c(stars, "**")
    } else if (pval < 0.05) {
      stars = c(stars, "*")
    } else {
      stars = c(stars, "")
    }
  }
  
  coefs_tbl = coefs_tbl %>%
    add_column(!!glue("({lag})") := glue("{signif(coefs, digits=5)}{stars}\n",
                                         "({signif(tvals, digits=5)})"))
}

# Save to HTML table -> path= does not work.
# coefs_tbl %>% gt() %>% gtsave("table_lsdv.html", path = output_path)
coefs_tbl %>% gt() %>% gtsave("table_lsdv.html")

# TODO: Model validation, e.g. walk-forward approach
