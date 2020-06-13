#### Setup ####
# Import standard variables
source("config.R")

library(readxl)
library(tidyverse)
library(lubridate)
library(reticulate)
library(zoo)

# Read in metadata
df_meta = readxl::read_xlsx(path_wiki, sheet = "Metadata")

# We need to clean the Wikipedia data before being able to process it in R, also
# to include new dates. Run the next line to do so (you may need to install
# Miniconda as a Python interpreter).
reticulate::py_run_file("clean_wide.py")

# Read in the cleaned Wikipedia data
df_wide = readr::read_csv(path_cleaned_wide, col_types = do.call(
  cols, list(Date=col_date(format="%Y-%m-%d"))))

#### Handle missing values ####
# We now add missing dates to the data for equal spacing.
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

# Complete the dates and execute the filling procedure described above
df_wide = df_wide %>%
  complete(Date = all_dates) %>%
  replace_na(named_cols_replace_0) %>%
  fill(all_of(cols_fill))

#### Create population variables ####
# We are currently only interested in regional (non-aggregated) data
df_wide = df_wide %>%
  select(Date:SAR_Deaths)

# To merge all Eurostat zip files, uncomment the next line if the file does not
# yet exist or if new data gets added.
# reticulate::py_run_file("eurostat_reader.py")

# We know the amount of people on January 1, 2019 as defined in df_eurostat. We
# only keep rows where the `region` is an Italian region, not a direction/NUTS-1
# region or the entire country by right joining with df_meta, since df_meta only
# contains NUTS-2 regions.

df_eurostat = readr::read_csv(path_full_eurostat, col_types = do.call(
  cols, list(region=col_character()))) %>%
  right_join(df_meta %>% select(Region, Code), by=c("region"="Region"))

# From https://www.worldometers.info/world-population/italy-population/, we find
# that the yearly growth rate for Italy in 2019 was -0.13% and for 2020 it was
# estimated to be -0.15% (not taking the coronacrisis into account). We assume
# that these rates are constant for all regions, for lack of a better metric.
growth_rate_pop_2019 = -0.0013
growth_rate_pop_2020 = -0.0015

# Find day number in 2020 of the day before the first date in the data. We take
# the day before because then we can factor in the amount of reported deaths due
# to COVID-19.
date_diff = as.integer(df_wide$Date[1] - as.Date("2020-01-01", "%Y-%m-%d")) - 1

for (regio in df_eurostat$Code) {
  # Only select the data for this region
  df_region = df_wide %>%
    select(starts_with(regio))
  
  # Initialize the baseline values at time t=1
  pop_region = df_eurostat %>%
    filter(Code == regio) %>%
    .[["populationNumbers"]]
  
  # We divide by 366 on the next line because 2020 is a leap year
  pop_base = pop_region * (1+growth_rate_pop_2019) *
    (1+growth_rate_pop_2020)^(date_diff/366) - 
    df_region[[paste0(regio, "_Deaths")]][1] 
  total_pop = pop_base - df_region[[paste0(regio, "_Deaths")]][1]
  suscept = pop_base - df_region[[paste0(regio, "_Confirmed")]][1] -
    df_region[[paste0(regio, "_Deaths")]][1]
  
  # The variables below are assumed to happen on day t at 5pm, namely the time
  # at which new data is reported. That is, the new cases, deaths, and
  # recoveries from that day contribute to the calculations rather than those of
  # the day before. This is to be consistent with the definitions from the
  # reports.
  for (t in 2:nrow(df_wide)) {
    
    # This column assumes the growth rates above and no extraordinary deaths due
    # to COVID-19. That is, what would be expected to happen if SARS-CoV-2 never
    # struck?
    pop_base = c(pop_base, (1+growth_rate_pop_2020)^(1/366)*tail(pop_base, 1))
    
    # The total population at time t is defined as the amount of people alive at
    # that time: total_pop(t) = total_pop(t-1) - deaths(t-1). In the code below,
    # because we have no other information, we assume that the deaths are only
    # those according to the coronavirus reports. Deaths statistics for Italy
    # include coronavirus victims who died in hospital, as well as those who
    # died outside of hospitals and were tested before or after dying. Indeed,
    # post-mortem tests are routinely carried out, and there is no distinction
    # between people who died "with" vs "because of" coronavirus, including
    # deaths of patients with pre-existing conditions. It should be said,
    # however, that in specific towns in regions where the healthcare system has
    # been overwhelmed by the pandemic (i.e. Lombardy), official death
    # statistics may have missed a portion of deaths outside hospitals.
    ## This is important because that means that it is likely that deaths due to
    ## comorbidities of the coronavirus were recorded under COVID-19 deaths.
    total_pop = c(total_pop,
                  tail(total_pop, 1) - df_region[[paste0(regio, "_Deaths")]][t])
    
    # Let:
    # I(t) = New confirmed cases at time t
    # S(t) = The total susceptible population at time t
    # The susceptible population at time t (beginning of the day) is then defined as:
    ## S(t) = S(t-1) - I(t-1)
    # This is because the people that die or recover have previously been tested
    # positive and are therefore already included in the confirmed cases (under
    # the assumption that only corona patients die). We also assume that other
    # causes of death are negligible. This can be relaxed later.
    suscept = c(suscept,
                tail(suscept, 1) - df_region[[paste0(regio, "_Confirmed")]][t])
  }
  
  df_wide = df_wide %>%
    add_column(!!paste0(regio, "_populationBaseline") := pop_base) %>%
    add_column(!!paste0(regio, "_totalPopulation") := total_pop) %>%
    add_column(!!paste0(regio, "_susceptiblePopulation") := suscept)
}

# Join the data with the extra data (containing the amount of active ICU
# patients, recoveries, tested people, and positively tested people)
df_wide = read_xlsx(path_wiki, sheet="Extra") %>%
  mutate(Date = as.Date(Date, format = "%Y-%m-%d")) %>%
  complete(Date = all_dates) %>%
  drop_na %>%
  full_join(df_wide, by="Date") %>%
  arrange(Date) # Sort by Date

# Note that, for example, cumsum(ABR_Confirmed) != ABR_TestedPositive. However,
# we do assume that the missing values, i.e. those before March 2, are correctly
# imputed by the cumsum of the confirmed cases, as long as the final element in
# the cumsum is less than the first non-missing element of TestedPositive.
# Unfortunately, we cannot (reasonably) impute _ICU. If the first known value
# for _Recovered and _Tested is 0, then we could potentially backpropagate this.

for (regio in df_eurostat$Code){
  # For completeness sake, even though the NAs (should) align, we find them per
  # region. It does not take much computing time
  which_NA = df_wide %>%
    select(!!paste0(regio, "_TestedPositive")) %>%
    is.na %>%
    which
  
  csum = df_wide %>%
    select(!!paste0(regio, "_Confirmed")) %>%
    cumsum %>%
    .[which_NA, ]
  
  first_elt = df_wide %>%
    select(!!paste0(regio, "_TestedPositive")) %>%
    na.omit %>%
    unlist %>%
    first
  
  # Check if the final element in the cumsum is lower than
  if (tail(csum, 1) <= first_elt) {
    df_wide[which_NA, paste0(regio, "_TestedPositive")] = csum
  } else {
    print(paste0("For region ", regio, ", the final element of csum (",
                 tail(csum, 1), ") was not lower than the first NA element (",
                 first_elt, "). We cannot impute values for this region."))
  }
  
  # In Adda's notation, region_Confirmed[t] is the incidence. We now convert
  # these to get the rates Inc and S. Do note that these will be extremely small
  df_wide = df_wide %>%
    mutate(!!paste0(regio, "_susceptibleRate") :=
             .data[[paste0(regio, "_susceptiblePopulation")]] /
             .data[[paste0(regio, "_totalPopulation")]]) %>%
    mutate(!!paste0(regio, "_incidenceRate") :=
             .data[[paste0(regio, "_Confirmed")]] /
             .data[[paste0(regio, "_totalPopulation")]])
}

# Sort the columns, keeping Date as the first column.
df_wide = df_wide[, c("Date", sort(colnames(df_wide)[-1]))]

# Convert the data to long format.
df_long_full = df_wide %>%
  pivot_longer(cols = -Date, names_to = c("Code", ".value"), names_sep = "_")

#### Process Eurostat regressors (import on line 60) ####
# Define the variables we are interested in
eurostat_variables = c("airPassengersArrived", "airPassengersDeparted",
                       "touristArrivals", "broadbandAccess",
                       "deathRateDiabetes", "deathRateInfluenza",
                       "deathRateChd", "deathRateCancer", 
                       "deathRatePneumonia", "availableBeds",
                       "maritimePassengersDisembarked",
                       "maritimePassengersEmbarked",
                       "riskOfPovertyOrSocialExclusion")
eurostat_variables = colnames(df_eurostat)

# Only keep rows where the `region` is an Italian region, not a direction or the
# entire country.
df_eurostat = df_eurostat[sapply(df_eurostat$region,
                                 function(x){x %in% df_meta$Region}), ] %>%
  select(c("region", "populationDensity", all_of(eurostat_variables)))

# Impute columns based on nearest neighbors tourists or population
na_cols = sapply(df_eurostat, function(x){x %>% is.na() %>% any()}) %>%
  .[.] %>%
  names

num_neighbors = 1 # Amount of neighbors to take into account when imputing

for (na_col in na_cols) {
  if (grepl("passenger", na_col, fixed=TRUE)) {
    diff_col = "touristArrivals"
  } else {
    diff_col = "populationDensity"
  }
  
  na_check = df_eurostat[[na_col]] %>% is.na()
  
  # Loop over regions/rows that are NA
  for (r in which(na_check)) {
    imputation = vector()
    
    df_eurostat_full = df_eurostat[!na_check, ]
    diffs = abs(df_eurostat_full[[diff_col]] - df_eurostat[[diff_col]][r])
    mins = sort(diffs)[1:num_neighbors]
    
    for (m in mins) {
      ratio = m / df_eurostat[[diff_col]][r]
      imputation = c(imputation,
                     df_eurostat_full[[na_col]][which(diffs == m)] *
                       ratio)
    }
    
    # We take the mean of the scaled 
    df_eurostat[r, na_col] = mean(imputation)
  }
}

# We need to expand the time-constant variables to the same scale as df_long
iddat = expand.grid(Date = unique(df_long_full$Date),
                    Code = unique(df_long_full$Code))
iddat = iddat[order(iddat$Date, iddat$Code), ]
rownames(iddat) <- NULL
df_eurostat = left_join(iddat, df_eurostat, by = "Code") %>% as_tibble

df_long_full = df_long_full %>%
  left_join(df_eurostat, by = c("Date", "Code"))

#### Process Google Mobility Report ####
df_gmr = readr::read_csv(path_mobility_report_official,
                         col_types = do.call(cols, list(
                           sub_region_2 = col_character(),
                           date = col_date(format = "%Y-%m-%d")))) %>%
  filter(country_region_code == "IT") %>%
  
  # Drop unused columns
  select(-country_region_code, -country_region, -sub_region_2) %>%
  
  # Drop rows with NAs in column `sub_region_1` (region name)
  drop_na(any_of("sub_region_1")) %>%
  
  # Clean column names
  rename_at(vars(ends_with("_percent_change_from_baseline")),
            funs(str_replace(., "_percent_change_from_baseline", "")))

colnames(df_gmr) = colnames(df_gmr) %>%
  str_replace("sub_region_1", "Region") %>%
  snakecase::to_upper_camel_case()

# Translate the regions to their Italian equivalent
df_gmr$Region = df_gmr$Region %>%
  str_replace("Aosta", "Valle d'Aosta/Vallée d'Aoste") %>%
  str_replace("Apulia", "Puglia") %>%
  str_replace("Lombardy", "Lombardia") %>%
  str_replace("Piedmont", "Piemonte") %>%
  str_replace("Sardinia", "Sardegna") %>%
  str_replace("Sicily", "Sicilia") %>%
  
  # Consider how to divide up in TN and BZ; for now: just the same (see below)
  str_replace("Trentino-South Tyrol", "Provincia Autonoma di Trento") %>%
  str_replace("Tuscany", "Toscana")

# Assume that the changes in Trentino-South Tyrol are the same for Trento and
# Bolzano
temp = df_gmr %>%
  filter(Region == "Provincia Autonoma di Trento")
temp$Region = "Provincia Autonoma di Bolzano/Bozen"

df_gmr = df_gmr %>%
  bind_rows(temp)

rm(temp) # Remove temp from the workspace

# Transform into percentages (decimal form)
df_gmr = df_gmr %>%
  mutate(across(-c(Region, Date), function(x) {x/100 + 1}))

# Add regional codes
df_gmr = df_meta %>%
  select(Region, Code) %>%
  right_join(df_gmr , by="Region") %>%
  select(-Region)

# We now are only interested in a decrease in the rail travellers, so we only
# select TransitStations.
df_gmr = df_gmr %>% select(Code, Date, TransitStations)

# For railroad transport, we can multiply by the baseline value. The latest data
# from Eurostat is from 2015. We assume that this value has changed in the same
# way as the population growth rate.
df_rail = readr::read_csv(path_railroad,
                          col_types = do.call(cols, list(
                            C_LOAD = col_character()))) %>%
  filter(TIME == max(TIME))

regions = df_gmr$Code %>% unique
baselines = vector()
date_diff = as.integer(df_gmr$Date[1] - as.Date("2020-01-01", "%Y-%m-%d")) - 1

for (region_code in regions) {
  region_name = df_meta[df_meta$Code == region_code, ][["Region"]]
  
  baseline = df_rail %>%
    # Add the amount of passengers that travelled FROM the region
    filter(C_LOAD == region_name) %>%
    select(-c("TIME", "C_LOAD")) %>%
    sum(na.rm = TRUE) +
    
    # Add the amount of passengers that travelled TO the region
    df_rail[[region_name]] %>%
    sum(na.rm = TRUE) -
    
    # We double counted the amount of passengers that travelled WITHIN the
    # region so we subtract this
    df_rail %>%
    filter(C_LOAD == region_name) %>%
    .[[region_name]]
  
  # We convert the amount of rail travellers in accordance with the population
  # growth rate. We divide by 366 on the next line because 2020 is a leap year
  baseline = baseline * (1+growth_rate_pop_2019) *
    (1+growth_rate_pop_2020)^(date_diff/366) - 
    head(df_wide[[glue("{region_code}_Deaths")]], 1)
  
  # The given numbers in df_rail are per year so we need daily numbers. This
  # depends on whether it is a leap year
  if (df_rail$TIME[1] %>% leap_year) {
    baseline = baseline / 366
  } else {
    baseline = baseline / 365
  }
  
  # TODO: Assume constant behaviour for the dates in df_long_full missing in 
  #### here ####
  dates_before = seq.Date(df_long_full %>%
                           filter(Code == region_code) %>%
                           .[["Date"]] %>%
                           min,
                         df_gmr %>%
                           filter(Code == region_code) %>%
                           .[["Date"]] %>%
                           min-1, by="day")
  dates_after = seq.Date(df_gmr %>%
                           filter(Code == region_code) %>%
                           .[["Date"]] %>%
                           max+1,
                         df_long_full %>%
                           filter(Code == region_code) %>%
                           .[["Date"]] %>%
                           max, by="day")
  
  df_gmr = df_gmr %>%
    bind_rows(tibble(Code = rep(region_code, length(dates_before)),
                     Date = dates_before,
                     TransitStations = rep(df_gmr %>%
                                             filter(Code == region_code) %>%
                                             select(TransitStations) %>%
                                             head(1) %>%
                                             unlist(use.names=FALSE),
                                           length(dates_before)))) %>%
    bind_rows(tibble(Code = rep(region_code, length(dates_after)),
                     Date = dates_after,
                     TransitStations = rep(df_gmr %>%
                                             filter(Code == region_code) %>%
                                             select(TransitStations) %>%
                                             tail(1) %>%
                                             unlist(use.names=FALSE),
                                           length(dates_after))))
  
  all_dates = seq.Date(min(df_long_full$Date), max(df_long_full$Date), by="day")
  missing_dates = all_dates[which(!df_long_full$Date %in%
                                    (df_gmr %>% filter(Code == region_code) %>%
                                    .[["Date"]]))][-1]
  
  if (length(missing_dates) > 0){
    # TODO: Impute missing values by the mean of the surrounding values
    df_gmr = df_gmr %>%
      bind_rows(tibble(
        Code = rep(region_code, length(missing_dates)),
        Date = missing_dates,
        TransitStations = rep(NA, length(missing_dates)))) 
  }
  
  # Sort now that new dates have been added
  df_gmr = df_gmr %>% arrange(Code, Date)
  
  # Impute possible NAs with the mean of the surrounding values
  df_gmr$TransitStations = df_gmr$TransitStations %>%
    (function(x){(na.locf(x) + rev(na.locf(rev(x))))/2})
  
  # Add baseline to a column baselines to be added to df_gmr
  num_rows = df_gmr %>% filter(Code == region_code) %>% nrow
  baselines = c(baselines, rep(baseline, num_rows))
}

df_gmr = df_gmr %>%
  transmute(Date = Date, Code = Code,
            railTravelers = round(TransitStations*eval(baselines)))

# Join the expanded eurostat data and the railway data with the long data.
df_long_full = df_long_full %>%
  left_join(df_gmr, by = c("Date", "Code"))

#### Final processing ####
# Pivot the long data to wide data
df_wide = df_long_full %>%
  pivot_wider(names_from = Code,
              values_from = all_of(colnames(df_long_full)[-c(1,2)]))

# Colnames are now of the form "variable_regionCode". We want to have them of
# the form "regionCode_variable".
colnames(df_wide) = colnames(df_wide) %>%
  sapply(function(s) {
    s %>% str_split("_", simplify = TRUE) %>% rev() %>% paste(collapse="_")
    }, USE.NAMES=FALSE)

#### Export to file ####
readr::write_csv(df_wide, path_full_wide)
readr::write_csv(df_long_full, path_full_long)
