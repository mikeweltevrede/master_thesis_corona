#### Setup ####
rm(list=ls())

library(readxl, quietly=TRUE)
library(tidyverse, quietly=TRUE)
library(tidyquant, quietly=TRUE)
library(reticulate, quietly=TRUE)
library(forecast, quietly=TRUE)

data_path = "data"

# Activate the Conda environment. If it does not exist yet, create it with the
# required packages.
env_name = "r-thesis_corona"

tryCatch(use_condaenv(env_name),
         error=function(e){
           requirements = scan(file="requirements.txt", what=character(),
                               quiet=TRUE) 
           conda_create(env_name, packages=requirements)
           })

# Read in metadata
df_meta = read_excel(paste0(data_path, "/italy_wikipedia.xlsx"),
                     sheet="Metadata")

# Read in the data containing per region aggregated statistics
df_region_aggregated = read_excel(
  paste0(data_path, "/italy_wikipedia.xlsx"), sheet="RegionPerDate") %>%
  mutate(Date = as.Date(Date, format = "%Y-%m-%d")) %>%
  drop_na()

# We need to clean the Wikipedia data before being able to process it in R, also
# to include new dates. Uncomment the next line to do so (you may need to
# install Miniconda as a Python interpreter).
# py_run_file("clean_wide.py")

# Read in the cleaned Wikipedia data
df_wide = read_csv(paste0(data_path, "/italy_wikipedia_cleaned.csv"),
                   col_types = do.call(cols,
                                       list(Date=col_date(format="%Y-%m-%d"))))

# We now add missing dates to the data for equal spacing. # The following dates
# will be filled in:
all_dates = seq.Date(min(df_wide$Date), max(df_wide$Date), by="day")
missing_dates = all_dates[which(!all_dates %in% df_wide$Date)][-1]
sprintf("The following %d dates are missing and will be added: %s",
        length(missing_dates), paste(missing_dates, collapse=", "))

# For the non-aggregated variables, we enter 0. For the totals, we use the
# previous value.
cols_replace_0 = c(colnames(df_wide)[
  2:(which(colnames(df_wide) == "Confirmed_New")-1)], "Confirmed_New",
  "Deaths_New")
named_cols_replace_0 = set_names(as.list(rep(0, length(cols_replace_0))),
                                cols_replace_0)

# The columns to backfill are the remaining columns minus "Date".
cols_fill = colnames(df_wide)[
  which(!colnames(df_wide) %in% cols_replace_0)][-1]

df_wide = df_wide %>%
  complete(Date = seq.Date(min(Date), max(Date), by="day")) %>%
  replace_na(named_cols_replace_0) %>%
  fill(all_of(cols_fill))

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

#### Simple time series models ####
ts_region = ts(df_region[2:ncol(df_region)])
ts_variable = ts_region[, paste0(region, "_ICU")]
na_growth = is.na(ts_variable)
ts_variable = ts_variable[!na_growth]

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
df_restrictions = read_excel(paste0(data_path, "/italy_wikipedia.xlsx"),
                             sheet="NationwideRestrictions") %>%
  mutate(Date = as.Date(Date, format = "%Y-%m-%d"))
df_restrictions = drop_na(df_restrictions)
df_restrictions = df_restrictions[!na_growth,]

ts_restrictions = ts(df_restrictions[2:ncol(df_restrictions)])
auto.arima(ts_variable, xreg=ts_restrictions[, "SchoolsClosed"])

# xreg = -1.2372
# -> LOM_ICU changes by -1.2372 point as SchoolsClosed changes to 1

#### Import Eurostat files ####

# To merge all Eurostat zip files, uncomment the next line if the file does not
# yet exist.
# py_run_file("eurostat_reader.py")

df_eurostat = read_csv(paste0(data_path, "/merged_eurostat.csv"),
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
