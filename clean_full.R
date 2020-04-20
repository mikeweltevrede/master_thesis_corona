#### Setup ####
rm(list=ls())

library(readxl, quietly=TRUE)
library(tidyverse, quietly=TRUE)
library(reticulate, quietly=TRUE)

# Import standard variables
source("config.R")

# Read in metadata
df_meta = readxl::read_excel(path_wiki, sheet = "Metadata")

# We need to clean the Wikipedia data before being able to process it in R, also
# to include new dates. Run the next line to do so (you may need to
# install Miniconda as a Python interpreter).
reticulate::py_run_file("clean_wide.py")

# Read in the cleaned Wikipedia data
df_wide = readr::read_csv(path_cleaned_wide, col_types = do.call(
  cols, list(Date=col_date(format="%Y-%m-%d"))))

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

# We are currently only interested in regional (non-aggregated) data
df_wide = df_wide %>%
  select(Date:SAR_Deaths)

# Add the population numbers per region. We know the amount of people on January
# 1, 2019 as defined in df_eurostat. We only keep rows where the `region` is an
# Italian region, not a direction/NUTS-1 region or the entire country.
df_eurostat = readr::read_csv(path_full_eurostat, col_types = do.call(
  cols, list(region=col_character()))) %>%
  select(c("region", "population_numbers")) %>%
  right_join(df_meta %>% select(c("Region", "Code")), by=c("region"="Region"))

# From https://www.worldometers.info/world-population/italy-population/, we see
# that the yearly growth rate for Italy in 2019 was -0.13% and for 2020 it is
# estimated to be -0.15% (not taking the coronacrisis into account). We assume
# that these rates are constant for all regions, for lack of a better metric.
growth_rate_pop_2019 = -0.0013
growth_rate_pop_2020 = -0.0015

date_diff = as.integer(df_wide$Date[1] - as.Date("2020-01-01", "%Y-%m-%d"))

for (regio in df_eurostat$Code) {
  df_region = df_wide %>%
    select(starts_with(regio))
  
  # Initialize the baseline values at time t=1
  pop_region = df_eurostat %>%
    filter(Code == regio) %>%
    .[["population_numbers"]]
  pop_base = pop_region * (1+growth_rate_pop_2019) *
    (1+growth_rate_pop_2020)^(date_diff/366) # 366 because 2020 is a leap year
  total_pop = pop_base
  suscept = pop_base - df_region[[paste0(regio, "_Confirmed")]][1]
  
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
df_wide_full = read_excel(path_wiki, sheet="Extra") %>%
  mutate(Date = as.Date(Date, format = "%Y-%m-%d")) %>%
  drop_na() %>%
  full_join(df_wide, by="Date") %>%
  arrange(Date) # Sort by Date

# Note that, for example, cumsum(ABR_Confirmed) != ABR_TestedPositive. However,
# we do assume that the missing values, i.e. those before March 2, are correctly
# imputed by the cumsum of the confirmed cases, as long as the final element in
# the cumsum is less than the first non-missing element of TestedPositive.
# Unfortunately, we cannot (reasonably) impute _ICU, _Recovered, and _Tested.
for (regio in df_eurostat$Code){
  # For completeness sake, even though the NAs (should) align, we find them per
  # region. It does not take much computing time
  which_NA = df_wide_full %>%
    select(!!paste0(regio, "_TestedPositive")) %>%
    is.na() %>%
    which()
  
  csum = df_wide_full %>%
    select(!!paste0(regio, "_Confirmed")) %>%
    cumsum() %>%
    .[which_NA, ]
  
  first_elt = df_wide_full %>%
    select(!!paste0(regio, "_TestedPositive")) %>%
    na.omit() %>%
    unlist() %>%
    first()
  
  # Check if the final element in the cumsum is lower than
  if (tail(csum, 1) <= first_elt) {
    df_wide_full[which_NA, paste0(regio, "_TestedPositive")] = csum
  } else {
    print(paste0("For region ", regio, ", the final element of csum (",
                 tail(csum, 1), ") was not lower than the first NA element (",
                 first_elt, "). We cannot impute values for this region."))
  }
  
  # In Adda's notation, region_Confirmed[t] is the incidence. We now convert
  # these to get the rates Inc and S. Do note that these will be extremely small
  df_wide_full = df_wide_full %>%
    mutate(!!paste0(regio, "_susceptibleRate") :=
             .data[[paste0(regio, "_susceptiblePopulation")]] /
             .data[[paste0(regio, "_totalPopulation")]]) %>%
    mutate(!!paste0(regio, "_incidenceRate") :=
             .data[[paste0(regio, "_Confirmed")]] /
             .data[[paste0(regio, "_totalPopulation")]])
}

# Sort the columns, keeping Date as the first column.
df_wide_full = df_wide_full[, c("Date", sort(colnames(df_wide_full)[-1]))]

# Convert the data to long format.
df_long_full = df_wide_full %>%
  pivot_longer(cols = -Date, names_to = c("Code", ".value"), names_sep = "_")

# Save the tibbles to a file
readr::write_csv(df_wide_full, path_full_wide)
readr::write_csv(df_long_full, path_full_long)
