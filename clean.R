#### Setup ####
# Import standard variables
source("config.R")

library(glue)
library(lubridate)
library(readxl)
library(reticulate)
library(tidyverse)
library(zoo)

# Import metadata
df_meta = readxl::read_xlsx(path_metadata, sheet = "Metadata")

# Column names are now of the form "variable_regionCode". We want to have them
# in the form "regionCode_variable".
pivot_to_df_wide = function(df_long) {
  # Turn long data into wide data
  df_wide = df_long %>%
    pivot_wider(names_from = code,
                values_from = all_of(colnames(df_long)[-c(1,2)]))
  
  # Column names are now of the form "variable_regionCode". We want to have them
  # in the form "regionCode_variable".
  colnames(df_wide) = colnames(df_wide) %>%
    sapply(function(s) {
      s %>%
        str_split("_", simplify = TRUE) %>%
        rev %>%
        paste(collapse="_")
    }, USE.NAMES=FALSE)
  
  # Order the columns alphabetically, keeping date at the start
  df_wide = df_wide[, c("date", sort(colnames(df_wide[-1])))]
  
  return(df_wide)
}

#### Import the data ####
# From the official Github account; make sure to pull the repository's latest
# version yourself!

# Check if the larger data sheet already exists
if (file.exists(new_data_path)) {
  # Then append new data to the existing data file
  
  # All dates for which reports are posted
  all_dates = list.files(data_path_gh) %>%
    stringr::str_extract("\\d{8}") 
  all_dates = all_dates[!is.na(all_dates)]
  
  # All dates for which the data has already been added
  completed_dates = scan(file=completed_dates_path, what=character(),
                         quiet=TRUE)
  
  # All dates for which reports are posted but for which the data has not yet
  # been added
  missing_dates = setdiff(completed_dates, all_dates)
  
  if (length(missing_dates) > 0){
    sprintf("The following %d dates are missing and will be added: %s",
            length(missing_dates), paste(missing_dates, collapse=", "))
    
    for (date in missing_dates){
      # Per date, import the data
      df = read_csv(glue("{data_path_gh}/dpc-covid19-ita-regioni-{date}.csv"),
                    col_types=do.call(
                      cols_only, # Only retain the columns specified below
                      list(
                        data=col_date(format="%Y-%m-%dT%H:%M:%S"), # Date
                        denominazione_regione=col_character(), # Region name
                        deceduti=col_integer(), # Total number of deceased
                        dimessi_guariti=col_integer(), # Total number recovered
                        totale_casi=col_integer(), # Total positive tests
                        tamponi=col_integer() # Total number of tests executed
                      ))) %>%
        rename(date = data,
               regionGH = denominazione_regione,
               deaths = deceduti,
               recovered = dimessi_guariti,
               infectives = totale_casi,
               tested = tamponi
        ) %>%
        
        # We want to use the region names in df_meta and not the one in the
        # official data
        left_join(df_meta %>% select(code, regionGH), by="regionGH") %>%
        select(-regionGH) %>%
        
        # Reorder the columns
        select(date, code, everything()) %>%
        
        # Append to the existing file
        write.table(file = new_data_path, sep = ",", append = TRUE, quote = FALSE,
                    col.names = FALSE, row.names = FALSE)
      
      # Append this date to the completed dates file
      write(date, file=completed_dates_path, append=TRUE)
    }
  }
} else{
  # Then process the entire file
  df = read_csv(glue("{data_path_gh}/dpc-covid19-ita-regioni.csv"),
                col_types=do.call(
                  cols_only, # Only retain the columns specified below
                  list(
                    data=col_date(format="%Y-%m-%dT%H:%M:%S"), # Date
                    denominazione_regione=col_character(), # region name
                    deceduti=col_integer(), # Total number of deceased
                    dimessi_guariti=col_integer(), # Total number recovered
                    totale_casi=col_integer(), # Total positive tests
                    tamponi=col_integer() # Total number of tests executed
                  ))) %>%
    rename(date = data,
           regionGH = denominazione_regione,
           deaths = deceduti,
           recovered = dimessi_guariti,
           infectives = totale_casi,
           tested = tamponi
    ) %>%
    
    # We want to use the region names in df_meta and not the one in the
    # official data
    left_join(df_meta %>% select(code, regionGH), by="regionGH") %>%
    select(-regionGH) %>%
    
    # Reorder the columns
    select(date, code, everything()) %>%
    
    # Save to a csv
    readr::write_csv(new_data_path)
  
  # Write the processed dates to the completed dates file - note: this
  # overwrites the existing file
  df$date %>%
    unique %>%
    str_replace_all("-", "") %>%
    write(file=completed_dates_path)
}

# df is only a temporary file and can be removed
tryCatch(rm(df), warning = function(cond) {})

#### Import the entire data ####
df_long = read_csv(new_data_path) %>%
  arrange(date, code)

# On June 12, Campania reported -229 confirmed cases, whereas the number of new
# cases in the week before that date only ranges from 0 to 5. Similarly, on June
# 19, Sicilia reported -394 confirmed cases and Trentino reported 387 new
# confirmed cases on June 24. This is likely to correct for the number of wrong
# confirmed cases reported in the past. We have no way to reliably determine how
# these were distributed and using our usual propagation method would lead to 0
# confirmed cases from May 13 till June 12 for Campania and Sicily,
# respectively.  Therefore, we drop these regions from our analysis.
df_long = df_long %>%
  filter(code != "CAM" & code != "SIC" & code != "TN")

# Construct group of removed population
df_long = df_long %>%
  mutate(removed = recovered + deaths)

# Take the first difference of the columns except for date
df_wide = pivot_to_df_wide(df_long) %>%
  mutate_at(vars(!matches("date")), function(x){c(NA,diff(x))}) %>%
  drop_na

# df_long currently contains totals. To keep these, we adapt the names
colnames(df_long) = c("date", "code", paste0(colnames(df_long %>%
                                                        select(!c(date, code))),
                                             "Total"))

df_long = df_wide %>%
  pivot_longer(cols = -date, names_to = c("code", ".value"), names_sep = "_") %>%
  left_join(df_long, by=c("date", "code")) %>% 
  group_by(code) %>% 
  mutate(
    activeInfectives =
      infectives - dplyr::lag(recovered, 1) - dplyr::lag(deaths, 1)
    ) %>%
  drop_na() %>%
  mutate(
    activeInfectives = cumsum(activeInfectives)
  ) %>%
  ungroup()

# The number of tests executed cannot be lower than the number of people tested
# positive. If this is the case, we set the number of tests equal to the number
# of positive tests.
## We need to do this later for the first differences too. That is, if from one
## day to the next, more people are tested positive than tests are executed,
## this is illogical and should be corrected.
index = df_long$tested < df_long$infectives
df_long[index, "tested"] = df_long[index, "infectives"]

df_wide = pivot_to_df_wide(df_long)

readr::write_csv(df_wide, new_data_path_wide)
readr::write_csv(df_long, new_data_path_long_cleaned)

# We need to clean the data before being able to process it in R, also
# to include new dates. Run the next line to do so (you may need to install
# Miniconda as a Python interpreter). This processes the negative numbers in
# propagating the value backwards until no negative numbers remain.
reticulate::py_run_file("clean_wide.py")

# Read in the cleaned data
df_wide = readr::read_csv(new_data_path_wide_cleaned, col_types = do.call(
  cols, list(date=col_date(format="%Y-%m-%d"))))

#### Handle missing values ####
# We now add missing dates to the data for equal spacing.
all_dates = seq.Date(min(df_wide$date), max(df_wide$date), by="day")
missing_dates = all_dates[which(!all_dates %in% df_wide$date)][-1]

if (length(missing_dates) > 0){
  sprintf("The following %d dates are missing and will be added: %s",
          length(missing_dates), paste(missing_dates, collapse=", "))

 # Complete the dates and execute the filling procedure described above
  df_wide = df_wide %>%
    complete(date = all_dates) %>%
    replace(is.na(.), 0)
}

#### Create population variables ####
# To merge all Eurostat zip files, uncomment the next line if the file does not
# yet exist or if new data gets added.
# reticulate::py_run_file("eurostat_reader.py")

# We know the amount of people on January 1, 2019 as defined in df_eurostat. We
# only keep rows where the `region` is an Italian region, not a direction/NUTS-1
# region or the entire country by right joining with df_meta, since df_meta only
# contains NUTS-2 regions.
df_eurostat = readr::read_csv(path_full_eurostat, col_types = do.call(
  cols, list(region=col_character()))) %>%
  right_join(df_meta %>% select(region, code), by="region")

# From https://www.worldometers.info/world-population/italy-population/, we find
# that the yearly growth rate for Italy in 2019 was -0.13% and for 2020 it was
# estimated to be -0.15% (not taking the coronacrisis into account). We assume
# that these rates are constant for all regions, for lack of a better metric.
growth_rate_pop_2019 = -0.0013
growth_rate_pop_2020 = -0.0015

# Find day number in 2020 of the day before the first date in the data. We take
# the day before because then we can factor in the amount of reported deaths due
# to COVID-19.
date_diff = as.integer(df_wide$date[1] - as.Date("2020-01-01", "%Y-%m-%d")) - 1

# Add population variables for the base and total population
for (regio in unique(df_long$code)) {
  # Only select the data for this region
  df_region = df_wide %>%
    select(starts_with(regio))
  
  # Initialize the baseline values at time t=1
  pop_region = df_eurostat %>%
    filter(code == regio) %>%
    .[["populationNumbers"]]
  
  # We divide by 366 on the next line because 2020 is a leap year
  pop_base = round(pop_region * (1+growth_rate_pop_2019) *
                     (1+growth_rate_pop_2020)^(date_diff/366) - 
                     df_region[[glue("{regio}_deaths")]][1])
  total_pop = pop_base
  
  # The population is constant over time; only deaths by COVID-19 are counted
  populationN = pop_base
  
  # The variables below are assumed to happen on day t at 5pm, namely the time
  # at which new data is reported. That is, the new cases, deaths, and
  # recoveries from that day contribute to the calculations rather than those of
  # the day before. This is to be consistent with the definitions from the
  # reports.
  for (t in 2:nrow(df_wide)) {
    
    # This column assumes the growth rates above and no extraordinary deaths due
    # to COVID-19. That is, what would be expected to happen if SARS-CoV-2 never
    # struck?
    pop_base = c(pop_base,
                 round((1+growth_rate_pop_2020)^(1/366)*tail(pop_base, 1)))
    
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
                  tail(total_pop, 1) - df_region[[glue("{regio}_deaths")]][t])
  }

  df_wide = df_wide %>%
    mutate(
      "{regio}_populationN" := populationN,
      "{regio}_populationBaseline" := pop_base,
      "{regio}_totalPopulation" := total_pop
    )
}

df_long = df_wide %>%
  pivot_longer(cols = -date, names_to = c("code", ".value"), names_sep = "_") %>%
  mutate(populationAlive = populationN - deaths)

#### Add undocumented infections ####
# We add undocumented infections because the susceptible population depends on
# the number of total infectives.
source("undocumented_infections.R")

# We assume some values for fmin and the gammas for the cubic form. For the
# quadratic form, we try several values of gamma below.
fmin = 0.1
gamma1_cubic = 0.6
gamma2_cubic = 0.8

df_long = df_long %>% 
  mutate(Smax = populationN - removedTotal) %>%
  rowwise() %>%
  # Compute f_t
  mutate(proportionDocumentedLinear =
           undocumented_infections(Smax, testedTotal,
                                   form = "Linear", fmin = fmin),
         proportionDocumentedDownwardsVertex =
           undocumented_infections(Smax, testedTotal,
                                   form = "DownwardsVertex", fmin = fmin),
         proportionDocumentedUpwardsVertex =
           undocumented_infections(Smax, testedTotal,
                                   form = "UpwardsVertex", fmin = fmin),
         proportionDocumentedCubic =
           undocumented_infections(Smax, testedTotal, form = "Cubic",
                                   gamma = gamma1_cubic, gamma2 = gamma2_cubic,
                                   fmin = fmin)
         ) %>%
  
  # Compute I_t = DI_t / f_t
  mutate(activeInfectivesLinear =
           round(activeInfectives / proportionDocumentedLinear),
         
         # TODO: For some reason, downwards vertex does not work, but gamma = 0.75 does work for quadratic
         activeInfectivesDownwardsVertex =
           round(activeInfectives / proportionDocumentedDownwardsVertex),
         activeInfectivesUpwardsVertex =
           round(activeInfectives / proportionDocumentedUpwardsVertex),
         activeInfectivesCubic =
           round(activeInfectives / proportionDocumentedCubic)
         ) %>%
  # NAs are formed when no tests have been executed. Instead of counting these
  # as NAs, TC_t=0 implies f_t=fmin. So, we replace the NAs with fmin.
  replace(is.na(.), fmin)

# For the Quadratic form, we try several values of gamma
gammas_quadratic = c(0.6, 0.65, 0.7, 0.75)  
names(gammas_quadratic) = c("Six", "SixtyFive", "", "SeventyFive") # 0.7 is default

for (i in 1:length(gammas_quadratic)) {
  name_gamma = names(gammas_quadratic)[i]
  df_long = df_long %>% 
    rowwise() %>%
    mutate("proportionDocumentedQuadratic{name_gamma}" :=
             undocumented_infections(Smax, testedTotal, form = "Quadratic",
                                     gamma = gammas_quadratic[[i]],
                                     fmin = fmin)
           ) %>%
    # Compute I_t = DI_t / f_t
    mutate("activeInfectivesQuadratic{name_gamma}" :=
             round(activeInfectives /
                     .data[[glue("proportionDocumentedQuadratic{name_gamma}")]])
    ) %>% 
    # NAs are formed when no tests have been executed. Instead of counting these
    # as NAs, TC_t=0 implies f_t=fmin. So, we replace the NAs with fmin.
    replace(is.na(.), fmin)
}

#### Compute the susceptible population ####
# Let:
# S(t) = The total susceptible population at time t
# N(t) = The total population at time t
# D(t) = New deaths between time t-1 and t
# I(t) = New confirmed cases between time t-1 and t
# R(t) = New recoveries between time t-1 and t

# All data is reported at 17:00 (5pm) on day t. The susceptible population
# at time t (so at 17:00) is then defined as:
# S(t) = N(t) - D(t) - I(t) - R(t)
for (form in c("Linear", "DownwardsVertex", "UpwardsVertex", "Cubic", "")) {
  
  infective_variable = glue("activeInfectives{form}")
  df_long = df_long %>%
    group_by(code) %>%
    mutate(
      "susceptiblePopulation{form}" :=
        populationN - removedTotal - .data[[infective_variable]],
      "susceptibleRate{form}" :=
        .data[[glue("susceptiblePopulation{form}")]] / .data[["populationAlive"]],
      "activeInfectivesTotal{form}" := cumsum(.data[[infective_variable]]),
      "infectivesRateTotal{form}" :=
        .data[[glue("activeInfectivesTotal{form}")]] / .data[["populationAlive"]],
      "infectivesRate{form}" :=
        .data[[infective_variable]] / .data[["populationAlive"]]
      ) %>%
    ungroup()
}

form = "Quadratic"
for (i in 1:length(gammas_quadratic)) {
  name_gamma = names(gammas_quadratic)[i]
  
  infective_variable = glue("activeInfectives{form}{name_gamma}")
  df_long = df_long %>%
    group_by(code) %>%
    mutate(
      "susceptiblePopulation{form}{name_gamma}" :=
        populationN - removedTotal - .data[[infective_variable]],
      "susceptibleRate{form}{name_gamma}" :=
        .data[[glue("susceptiblePopulation{form}{name_gamma}")]] /
        .data[["populationAlive"]],
      "activeInfectivesTotal{form}{name_gamma}" :=
        cumsum(.data[[infective_variable]]),
      "infectivesRateTotal{form}{name_gamma}" :=
        .data[[glue("activeInfectivesTotal{form}{name_gamma}")]] / 
        .data[["populationAlive"]],
      "infectivesRate{form}" :=
        .data[[infective_variable]] / .data[["populationAlive"]]
    ) %>%
    ungroup()
}
  
df_wide = pivot_to_df_wide(df_long)

#### Process Eurostat regressors (import on line 213) ####
# Define the variables we are interested in
eurostat_variables = c("touristArrivals", "broadbandAccess",
                       "dischargeRateDiabetes", "dischargeRateHypertension",
                       "dischargeRateCancer", "dischargeRateChd",
                       "dischargeRateTB", "availableBeds",
                       "riskOfPovertyOrSocialExclusion", "medianAge",
                       "populationDensity")
eurostat_variables = colnames(df_eurostat)

# Only keep rows where the `region` is an Italian region, not a direction or the
# entire country.
df_eurostat = df_eurostat[sapply(df_eurostat$region,
                                 function(x){x %in% df_meta$region}), ] %>%
  select(c("region", "code", all_of(eurostat_variables)))

# Impute columns based on nearest neighbors tourists or population
na_cols = sapply(df_eurostat, function(x){x %>% is.na() %>% any()}) %>%
  .[.] %>%
  names

num_neighbors = 1 # Amount of neighbors to take into account when imputing

for (na_col in na_cols) {
  if (grepl("passenger", na_col, fixed=TRUE)) {
    # If the variable contains "passenger", we use "touristArrivals" as the
    # matching variable
    diff_col = "touristArrivals"
  } else {
    # Else, we use "populationDensity"
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
iddat = expand.grid(date = unique(df_long$date),
                    code = unique(df_long$code))
iddat = iddat[order(iddat$date, iddat$code), ]
rownames(iddat) <- NULL
df_eurostat = left_join(iddat, df_eurostat, by = "code") %>% as_tibble

df_long = df_long %>%
  left_join(df_eurostat, by = c("date", "code"))

#### Process Google Mobility Report ####
# Grocery & pharmacy: Mobility trends for places like grocery markets, food
#   warehouses, farmers markets, specialty food shops, drug stores, and
#   pharmacies.
# Parks: Mobility trends for places like local parks, national parks, public
#   beaches, marinas, dog parks, plazas, and public gardens.
# Transit stations: Mobility trends for places like public transport hubs such
#   as subway, bus, and train stations.
# Retail & recreation: Mobility trends for places like restaurants, cafes,
#   shopping centers, theme parks, museums, libraries, and movie theaters.
# Residential: Mobility trends for places of residence.
# Workplaces: Mobility trends for places of work.

if (file.exists(path_mobility_report_cleaned)){
  df_gmr = readr::read_csv(path_mobility_report_cleaned,
                           col_types = do.call(
                             cols_only, list(
                               country_region_code=col_character(),
                               sub_region_1=col_character(),
                               sub_region_2=col_character(),
                               date=col_date(format = "%Y-%m-%d"),
                               transit_stations_percent_change_from_baseline=
                                 col_double())))
} else {
  df_gmr = readr::read_csv(path_mobility_report_official,
                           col_types = do.call(
                             cols_only, list(
                               country_region_code=col_character(),
                               sub_region_1=col_character(),
                               sub_region_2=col_character(),
                               date=col_date(format = "%Y-%m-%d"),
                               transit_stations_percent_change_from_baseline=
                                 col_double()))) %>%
    filter(country_region_code == "IT") %>%
    
    # Google recently added the NUTS 3 regions. We want to only have the
    # region-wide (NUTS 2 region) data
    filter(is.na(sub_region_2))

  write_csv(df_gmr, path_mobility_report_cleaned)
}

df_gmr = df_gmr %>%
  # Drop unused columns
  select(-country_region_code) %>%
  
  # Drop rows with NAs in column `sub_region_1` (region name)
  drop_na(any_of("sub_region_1")) %>%
    
  # Clean column names
  rename_at(vars(ends_with("_percent_change_from_baseline")),
            function(x){str_replace(x, "_percent_change_from_baseline", "")})

colnames(df_gmr) = colnames(df_gmr) %>%
  str_replace("sub_region_1", "region") %>%
  snakecase::to_lower_camel_case()

# Translate the regions to their Italian equivalent
df_gmr$region = df_gmr$region %>%
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
  filter(region == "Provincia Autonoma di Trento")
temp$region = "Provincia Autonoma di Bolzano/Bozen"

df_gmr = df_gmr %>%
  bind_rows(temp)

rm(temp) # Remove temp from the workspace

# Replace region name by regional codes
df_gmr = df_meta %>%
  select(region, code) %>%
  right_join(df_gmr, by="region") %>%
  select(-region) %>%
  
  # Make sure that the regions concur with df_long
  filter(code %in% unique(df_long$code)) %>%
  
  # We now are only interested in a decrease in the rail travelers, so we only
  # select transitStations. Note that this is not only for train stations but
  # similar places like public transport hubs.
  # TODO: Use all variables?
  select(code, date, transitStations)

# Transform into percentages (decimal form)
df_gmr = df_gmr %>%
  mutate(across(-c(code, date), function(x) {x/100 + 1}))

# For railroad transport, we can multiply by the baseline value. The latest data
# from Eurostat is from 2015. We assume that this value has changed in the same
# way as the population growth rate.
df_rail = readr::read_csv(path_railroad,
                          col_types = do.call(cols, list(
                            C_LOAD = col_character()))) %>%
  filter(TIME == max(TIME))

# Since the amount of rail travellers was *in* 2015, this was reported at the
# end of the year. As such, we need to start in 2016
growth_rate_pop_2016 = 0.0014
growth_rate_pop_2017 = 0.0002
growth_rate_pop_2018 = -0.0008

baselines = vector()
date_diff = as.integer(df_gmr$date[1] - as.Date("2020-01-01", "%Y-%m-%d")) - 1
all_dates = seq.Date(min(df_long$date), max(df_long$date), by="day")

for (region_code in unique(df_gmr$code)) {
  if (is.na(region_code)){
    next
  }
  
  region_name = df_meta[df_meta$code == region_code, ][["region"]]
  
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
  baseline = baseline * (1+growth_rate_pop_2016) * (1+growth_rate_pop_2017) *
    (1+growth_rate_pop_2018) * (1+growth_rate_pop_2019) *
    (1+growth_rate_pop_2020)^(date_diff/366) - 
    head(df_wide[[glue("{region_code}_deaths")]], 1)
  
  # The given numbers in df_rail are per year so we need daily numbers. This
  # depends on whether it is a leap year
  if (df_rail$TIME[1] %>% leap_year) {
    baseline = baseline / 366
  } else {
    baseline = baseline / 365
  }
  
  # Assume constant behaviour before and after the limiting dates, if applicable
  min_date_long = df_long %>%
    filter(code == region_code) %>%
    .[["date"]] %>%
    min
  max_date_long = df_long %>%
    filter(code == region_code) %>%
    .[["date"]] %>%
    max
  min_date_gmr = df_gmr %>%
    filter(code == region_code) %>%
    .[["date"]] %>%
    min
  max_date_gmr = df_gmr %>%
    filter(code == region_code) %>%
    .[["date"]] %>%
    max
  
  if (min_date_long < min_date_gmr){
    # Then we do need to draw out the values before the minimum date
  
    dates_before = seq.Date(min_date_long, min_date_gmr-1, by="day")
    df_gmr = df_gmr %>%
      bind_rows(tibble(code = rep(region_code, length(dates_before)),
                       date = dates_before,
                       transitStations = rep(df_gmr %>%
                                               filter(code == region_code) %>%
                                               select(transitStations) %>%
                                               head(1) %>%
                                               unlist(use.names=FALSE),
                                             length(dates_before))))
  }
  
  if (max_date_long > max_date_gmr){
    # Then we do need to draw out the values after the maximum date
    dates_after = seq.Date(max_date_gmr+1, max_date_long, by="day")
    
    df_gmr = df_gmr %>%
      bind_rows(tibble(code = rep(region_code, length(dates_after)),
                       date = dates_after,
                       transitStations = rep(df_gmr %>%
                                               filter(code == region_code) %>%
                                               select(transitStations) %>%
                                               tail(1) %>%
                                               unlist(use.names=FALSE),
                                             length(dates_after))))
  }
  
  missing_dates = all_dates[which(!unique(df_long$date) %in%
                                    (df_gmr %>% filter(code == region_code) %>%
                                       .[["date"]]))][-1]
  
  if (length(missing_dates) > 0) {
    df_gmr = df_gmr %>%
      bind_rows(tibble(
        code = rep(region_code, length(missing_dates)),
        date = missing_dates,
        transitStations = rep(NA, length(missing_dates)))) 
  }
  
  # Sort now that new dates have been added
  df_gmr = df_gmr %>% arrange(code, date)
  
  # Impute possible NAs with the mean of the surrounding values
  df_gmr$transitStations = df_gmr$transitStations %>%
    (function(x){(zoo::na.locf(x) + rev(zoo::na.locf(rev(x))))/2})
  
  # Add baseline to a column baselines to be added to df_gmr
  num_rows = df_gmr %>% filter(code == region_code) %>% nrow
  baselines = c(baselines, rep(baseline, num_rows))
}

df_gmr = df_gmr %>%
  transmute(date = date, code = code,
            transitStations = round(transitStations*eval(baselines) / 366))

# Join the railway data with the long data
df_long = df_long %>%
  left_join(df_gmr, by = c("date", "code")) %>%
  mutate(transitStations = transitStations / totalPopulation)

#### Final processing ####
df_long = df_long %>%
  mutate(populationDensity = totalPopulation/area)

# The number of tests executed cannot be lower than the number of people tested
# positive. We have now taken first differences and this still holds. That is,
# if from one day to the next, more people are tested positive than tests are
# executed, this is illogical and should be corrected. If this is the case, we
# set the number of tests equal to the number of positive tests. 
index = df_long$tested < df_long$infectives
df_long[index, "tested"] = df_long[index, "infectives"]

# Pivot the long data to wide data
df_wide = pivot_to_df_wide(df_long)

# Add nationwide variables by summing the individual regions' variables as well
# as weekend and lockdown dummies
lockdown_start = "2020-03-10"
lockdown_end = "2020-05-18"

totalPopulationNational = df_wide %>%
  select(ends_with("totalPopulation")) %>%
  rowSums
populationAliveNational = df_wide %>%
  select(ends_with("populationAlive")) %>%
  rowSums
areaNational = df_wide %>%
  select(ends_with("area")) %>%
  rowSums
df_wide = df_wide %>%
  mutate(totalPopulationNational = totalPopulationNational,
         populationAliveNational = populationAliveNational,
         populationDensityNational = populationAliveNational/areaNational,
         weekend =
           lubridate::wday(date, label = TRUE) %in% c("Sat", "Sun") %>%
           as.integer %>% as.factor,
         lockdown =
           ifelse(date > as.Date(lockdown_start, format = "%Y-%m-%d") &
                    date < as.Date(lockdown_end, format = "%Y-%m-%d"),
                  1, 0) %>%
           as.factor)

for (form in c("Linear", "DownwardsVertex", "UpwardsVertex", "Cubic", "")) {
  
  infective_variable = glue("activeInfectives{form}")
  
  susceptiblePopulationNational = df_wide %>%
    select(ends_with(glue("susceptiblePopulation{form}"))) %>%
    rowSums
  activeInfectivesNational = df_wide %>% 
    select(ends_with(glue("_{infective_variable}"))) %>% 
    rowSums
  
  df_wide = df_wide %>%
    mutate("susceptiblePopulationNational{form}" :=
             susceptiblePopulationNational,
           "susceptibleRateNational{form}" :=
             susceptiblePopulationNational / populationAliveNational,
           "activeInfectivesNational{form}" := activeInfectivesNational,
           "activeInfectivesTotalNational{form}" :=
             cumsum(activeInfectivesNational),
           "activeInfectivesRateNational{form}" :=
             !!sym(glue("activeInfectivesNational{form}")) /
             populationAliveNational,
           "activeInfectivesRateTotalNational{form}" :=
             !!sym(glue("activeInfectivesTotalNational{form}")) /
             populationAliveNational)
}

form = "Quadratic"
for (i in 1:length(gammas_quadratic)) {
  name_gamma = names(gammas_quadratic)[i]
  
  infective_variable = glue("activeInfectives{form}{name_gamma}")
  
  susceptiblePopulationNational = df_wide %>%
    select(ends_with(glue("susceptiblePopulation{form}{name_gamma}"))) %>%
    rowSums
  infectivesNational = df_wide %>% 
    select(ends_with(glue("_{infective_variable}"))) %>% 
    rowSums
  
  df_wide = df_wide %>%
    mutate("susceptiblePopulationNational{form}{name_gamma}" :=
             susceptiblePopulationNational,
           "susceptibleRateNational{form}{name_gamma}" :=
             susceptiblePopulationNational / populationAliveNational,
           "activeInfectivesNational{form}{name_gamma}" :=
             infectivesNational,
           "activeInfectivesTotalNational{form}{name_gamma}" :=
             cumsum(infectivesNational),
           "activeInfectivesRateNational{form}{name_gamma}" :=
             !!sym(glue("activeInfectivesNational{form}{name_gamma}")) /
             populationAliveNational,
           "activeInfectivesRateTotalNational{form}{name_gamma}" :=
             !!sym(glue("activeInfectivesTotalNational{form}{name_gamma}")) /
             populationAliveNational)
}

df_long = df_long %>%
  mutate(weekend =
           lubridate::wday(date, label = TRUE) %in% c("Sat", "Sun") %>%
           as.integer %>% as.factor,
         lockdown =
           ifelse(date > as.Date(lockdown_start, format = "%Y-%m-%d") &
                    date < as.Date(lockdown_end, format = "%Y-%m-%d"),
                  1, 0) %>%
           as.factor)

#### Export to file ####
readr::write_csv(df_wide, path_full_wide)
readr::write_csv(df_long, path_full_long)
