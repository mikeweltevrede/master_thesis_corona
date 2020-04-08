rm(list=ls())

library(tidyverse, verbose=FALSE)
library(tidyquant, verbose=FALSE)
library(reticulate, verbose=FALSE)
library(forecast, verbose=FALSE)

data_path = "data"

# Read in metadata
df_meta = readxl::read_excel(paste0(data_path, "/italy_wikipedia.xlsx"),
                             sheet="Metadata")

# We need to clean the Wikipedia data before being able to process it in R.
# Uncomment the next line to do so (you may need to install Miniconda).
# py_run_file("clean_wide.py")

# Read in the cleaned Wikipedia data
df_wide = readr::read_csv(
  paste0(data_path, "/italy_wikipedia_cleaned.csv"),
  col_types = do.call(cols, list(Date=col_date(format="%Y-%m-%d"))))

# We first are only interested in one region
region = "LOM"
region_full = df_meta$Region[which(df_meta$Code == region)]

region_cols = colnames(df_wide)[
  sapply(colnames(df_wide), function(x){grepl(region, x, fixed=TRUE)})]

# We do not have the number of active cases per region and this is difficult to
# compute since there is no data on the number of recoveries per region. For
# now, we will assume that the number of recoveries is 0.

# We can then later assume that the number of recoveries per region is
# proportional to the number of confirmed cases (e.g. if Lombardy has 40% of the
# confirmed cases on that day (or with a certain lag), it gets 40% of the
# nationwide recoveries).

df_region = df_wide %>%
  select(c("Date", region_cols)) %>%
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
  rename("Growth_Rate" = paste0(region, "_Confirmed1")) %>%
  mutate(Growth_Rate_MA3 = rollapply(Growth_Rate, 3, mean, align='right',
                                     fill=NaN))

#### Simple time series models ####
sma_obj = SMA(df_region$Growth_Rate, n=3) # Simple Moving Average
ema_obj = EMA(df_region$Growth_Rate, n=3) # Exponential Moving Average (note: unstable in the short-term)

ts_region = ts(df_region[2:ncol(df_region)])
growth_rate = ts_region[, "Growth_Rate"]
na_growth = is.na(growth_rate)
growth_rate = growth_rate[!na_growth]

Box.test(growth_rate, lag=14, type="Lj") # p-value = 2.543e-07
ggAcf(df_region$Growth_Rate, na.action = na.pass)
ar_obj = ar(df_region$Growth_Rate, na.action = na.pass)

gr_naive = naive(growth_rate, h=3)
autoplot(gr_naive)

# Naive
fc = naive(growth_rate, h=3)
checkresiduals(fc) # Not Normal; also lack of data; so be careful with prediction intervals
autoplot(fc, series="Data") +
  autolayer(fitted(fc), series="Fitted")

# Simple Exponential Smoothing
fc = ses(growth_rate, h=3)
checkresiduals(fc)
autoplot(fc, series="Data") +
  autolayer(fitted(fc), series="Fitted")

# Holt's trend method
fc = holt(growth_rate, h=3)
checkresiduals(fc)
autoplot(fc, series="Data") +
  autolayer(fitted(fc), series="Fitted")

# Dampened trend method
fc = holt(growth_rate, h=3, damped=TRUE)
checkresiduals(fc)
autoplot(fc, series="Data") +
  autolayer(fitted(fc), series="Fitted")

# Holt-Winters' methods are not applicable because we do not have seasonality

# General model: ETS(M,Ad,N): Multiplicative errors, dampened trend, no seasonality
ets(growth_rate)
fc = growth_rate %>%
  ets() %>%
  forecast(h=3)
checkresiduals(fc)
autoplot(fc, series="Data") +
  autolayer(fitted(fc), series="Fitted")

# ARIMA model - ARIMA(0,2,2); so we take 2 differences to create stationarity and include 2 lagged errors
BoxCox.lambda(growth_rate) # -0.5473554
auto.arima(growth_rate)
fc = growth_rate %>%
  auto.arima() %>%
  forecast(h=3)
checkresiduals(fc)
autoplot(fc, series="Data") +
  autolayer(fitted(fc), series="Fitted")

#### Advanced models ####
# Dynamic regression
df_restrictions = readxl::read_excel(paste0(data_path, "/italy_wikipedia.xlsx"),
                                     sheet="NationwideRestrictions")
df_restrictions = drop_na(df_restrictions)
df_restrictions = df_restrictions[!na_growth,]

ts_restrictions = ts(df_restrictions[2:ncol(df_restrictions)])
auto.arima(growth_rate, xreg=ts_restrictions[, "SchoolsClosed"])

# xreg = -0.0098 -> growth_rate changes by -0.98 percentage point as
# SchoolsClosed changes by 1 percentage point (how to interpret?)

#### Import Eurostat files ####

# To merge all Eurostat zip files, uncomment the next line if the file does not
# yet exist.
py_run_file("eurostat_reader.py")

df_eurostat = readr::read_csv(
  paste0(data_path, "/merged_eurostat.csv"),
  col_types = do.call(cols, list(region=col_character())))

# Only keep rows where the `region` is an Italian region, not a direction or the
# entire country.
df_eurostat = df_eurostat[sapply(df_eurostat$region,
                                 function(x){x %in% df_meta$Region}), ]

# Select the region we are interested in and replicate the data to have the same
# amount of rows as the time series data
df_eurostat_region = df_eurostat[which(df_eurostat$region == region_full), ] %>%
  select(-c("region")) %>%
  slice(rep(1:n(), each = nrow(df_wide)))


