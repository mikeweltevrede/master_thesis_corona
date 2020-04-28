#### Setup ####
rm(list=ls())

library(plm, quietly=TRUE)
library(tidyverse, quietly=TRUE)
library(glue)

# Import standard variables
source("config.R")

# Each day, the data should be recleaned
# source("clean_full.R")

df_long = readr::read_csv(path_full_long, col_types = do.call(
  cols, list(Date = col_date(format = "%Y-%m-%d"))))

regressors = c("airPassengersArrived", "touristArrivals", "broadbandAccess",
               "deathRateDiabetes", "deathRateInfluenza", "deathRateChd",
               "deathRateCancer", "deathRatePneumonia", "availableBeds",
               "maritimePassengersDisembarked",
               "riskOfPovertyOrSocialExclusion", "railTravelers")

#### Transform variables into proportions ####
# TODO: Transform variables to proportions (in clean_full.R or here?)
make_prop = function(x, na.rm = FALSE) { x / sum(x, na.rm = na.rm) }

regressors_prop = c("airPassengersArrived", "touristArrivals",
                    "deathRateDiabetes", "deathRateInfluenza", "deathRateChd",
                    "deathRateCancer", "deathRatePneumonia", "availableBeds",
                    "maritimePassengersDisembarked")

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

lag = 1

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

df_long_pd = pdata.frame(df_long %>%
                           select(all_of(c(base_vars, regressors))) %>%
                           arrange(Code),
                         index=c("Code","Date"), drop.index=TRUE,
                         row.names=TRUE)

plmo_pooled = plm(fm, data = df_long_pd, model = "pooling")
plmo_fe = plm(fm, data = df_long_pd, model = "within")
lsdv = lm(fm_lsdv, data=df_long)

pooltest(plmo_pooled, plmo_fe) # p-value = 5.463e-10 => H1: unstability => FE
summary(plmo_pooled)
summary(plmo_fe)
summary(lsdv)

# Check NAs
df_long %>% is.na %>% apply(2, which)

# TODO: Model validation, e.g. walk-forward approach
