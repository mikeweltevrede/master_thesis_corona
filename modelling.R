#### Setup ####
rm(list=ls())

library(plm, quietly=TRUE)
library(readxl, quietly=TRUE)
library(tidyverse, quietly=TRUE)

# Import standard variables
source("config.R")

# Run `clean_full.R` to create the file at `path_full_long` imported on the
# next line.
df_long = readr::read_csv(path_full_long, col_types = do.call(
  cols, list(Date = col_date(format = "%Y-%m-%d"))))

regressors = c("airPassengersArrived", "airPassengersDeparted",
               "touristArrivals", "broadbandAccess", "deathRateDiabetes",
               "deathRateInfluenza", "deathRateChd", "deathRateCancer",
               "deathRatePneumonia", "availableBeds",
               "maritimePassengersDisembarked", "maritimePassengersEmbarked",
               "riskOfPovertyOrSocialExclusion", "railTravelers")

#### Transform variables into proportions ####
# TODO: Transform variables to proportions (in clean_full.R or here?)
make_prop = function(x, na.rm = FALSE) { x / sum(x, na.rm = na.rm) }


# TODO: Run models
# Create formula to be the product of the regressors with the lagged incidence
# and susceptible rates
fm = paste("incidenceRate ~",
           "dplyr::lag(incidenceRate, 1):dplyr::lag(susceptibleRate, 1):" %>%
  paste0(regressors) %>%
  paste(collapse="+")) %>%
  as.formula()

plmo_pooled = plm(fm, data = df_long, model = "pooling")
plmo_fe = plm(fm, data = df_long, model = "within")
pooltest(plmo_pooled, plmo_fe) # p-value < 2.2e-16 => H1: unstability => FE
