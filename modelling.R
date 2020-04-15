#### Setup ####
rm(list=ls())

library(readxl, quietly=TRUE)
library(tidyverse, quietly=TRUE)
library(tidyquant, quietly=TRUE)
library(reticulate, quietly=TRUE)
library(forecast, quietly=TRUE)

# Import standard variables
source("config.R")

# Read in metadata
df_meta = read_excel(path_wiki, sheet="Metadata")

# Read in the data containing per region aggregated statistics
df_region_aggregated = read_excel(path_wiki, sheet="RegionPerDate") %>%
  mutate(Date = as.Date(Date, format = "%Y-%m-%d")) %>%
  drop_na()

# We need to clean the Wikipedia data before being able to process it in R, also
# to include new dates. For this, run `clean_wide_full.R` to create the file
# `path_full_wide` imported on the next line.
df_wide = read_csv(path_full_wide, col_types = do.call(
  cols, list(Date=col_date(format="%Y-%m-%d"))))

#### One region ####
region = "LOM"
region_full = df_meta$Region[which(df_meta$Code == region)]

# We do not have the number of active cases per region and this is difficult to
# compute since there is no data on the number of recoveries per region. For
# now, we will assume that the number of recoveries is 0.

# We later assume that the number of recoveries per region is proportional to
# the number of confirmed cases (e.g. if Lombardy has 10% of their total
# confirmed cases on March 1, then the number of recoveries on March 1 is set to
# be 10% of the total recoveries in Lombardy).
region_cols = colnames(df_wide)[
  sapply(colnames(df_wide), function(x){grepl(region, x, fixed=TRUE)})]

df_region = df_wide %>%
  select(c("Date", region_cols)) 
df_region = df_region %>%
  bind_cols(
    # Add total confirmed cases for this region (cumulative sum)
    df_region %>%
      select(ends_with("Confirmed")) %>%
      cumsum() %>%
      rename("Total_Confirmed" =  paste0(region, "_Confirmed"))
    ) %>%
  bind_cols(
    # Add number of active cases, defined as the total number of cases minus the
    # number of deaths, i.e. we assume no recoveries
    .$Total_Confirmed - select(., contains("_Deaths")) %>%
      as_tibble() %>%
      rename("Active" =  paste0(region, "_Deaths"))
    ) %>%
  bind_cols(
    # Add the growth rate
    df_region %>%
      select(paste0(region, "_Confirmed")) / .$Active
    ) %>%
  rename("Growth_Rate" = paste0(region, "_Confirmed1"))

#### Simple forecasting models ####
ts_region = ts(df_region[2:ncol(df_region)])
ts_variable = ts_region[, paste0(region, "_ICU")]
na_tsvar = is.na(ts_variable)
ts_variable = ts_variable[!na_tsvar]

Box.test(ts_variable, lag=14, type="Lj") # p-value = 2.543e-07
ggAcf(ts_variable, na.action = na.pass)
ar_obj = ar(ts_variable, na.action = na.pass)

# Forecast period
h = 5

# Naive
fc = naive(ts_variable, h=h)
checkresiduals(fc)
# Not Normal, so be careful with prediction intervals
# Not white noise!
autoplot(fc, series="Data") +
  autolayer(fitted(fc), series="Fitted")

# Simple Exponential Smoothing
fc = ses(ts_variable, h=h)
checkresiduals(fc) # Not white noise!
autoplot(fc, series="Data") +
  autolayer(fitted(fc), series="Fitted")

# Holt's trend method
fc = holt(ts_variable, h=h)
checkresiduals(fc)
autoplot(fc, series="Data") +
  autolayer(fitted(fc), series="Fitted")

# Dampened trend method
fc = holt(ts_variable, h=h, damped=TRUE)
checkresiduals(fc)
autoplot(fc, series="Data") +
  autolayer(fitted(fc), series="Fitted")

# Holt-Winters' methods are not applicable because we do not have seasonality

# General model: ETS(M,Ad,N): Multiplicative errors, dampened trend, no
# seasonality
ets(ts_variable)
fc = ts_variable %>%
  ets() %>%
  forecast(h=h)
checkresiduals(fc)
autoplot(fc, series="Data") +
  autolayer(fitted(fc), series="Fitted")

# ARIMA model - ARIMA(0,2,2); so we take 2 differences to create stationarity
# and include 2 lagged errors
BoxCox.lambda(ts_variable) # 1.387894 (LOM_ICU)
auto.arima(ts_variable)
fc = ts_variable %>%
  auto.arima() %>%
  forecast(h=h)
checkresiduals(fc)
autoplot(fc, series="Data") +
  autolayer(fitted(fc), series="Fitted")

#### Advanced models ####
# Dynamic regression
df_restrictions = read_excel(path_wiki, sheet="NationwideRestrictions") %>%
  mutate(Date = as.Date(Date, format = "%Y-%m-%d")) %>%
  drop_na()
df_restrictions = df_restrictions[!na_tsvar,]

ts_restrictions = ts(df_restrictions[2:ncol(df_restrictions)])
auto.arima(ts_variable, xreg=ts_restrictions[, "SchoolsClosed"])

# xreg = -1.2372
# -> LOM_ICU changes by -1.2372 point as SchoolsClosed changes to 1
