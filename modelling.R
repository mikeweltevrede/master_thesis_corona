#### Setup ####
rm(list=ls())

library(plm, quietly=TRUE)
library(tidyverse, quietly=TRUE)
library(glue)
library(gt)

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

lags = 1:14 # Incubation period
lsdv_results = vector("list")

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
  plmo_pooled = plm(fm, data = df_long_pd, model = "pooling")
  print(summary(plmo_pooled))
  
  plmo_fe = plm(fm, data = df_long_pd, model = "within")
  print(summary(plmo_fe))
  
  lsdv = lm(fm_lsdv, data=df_long)
  lsdv_results[[as.character(lag)]] = lsdv
  print(summary(lsdv))
}

# Create HTML table
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
    add_column(!!glue("({lag})") := glue("{coefs}{stars}\n({tvals})"))
}

coefs_tbl %>% gt() %>% gtsave("table_lsdv.html")

# TODO: Model validation, e.g. walk-forward approach
