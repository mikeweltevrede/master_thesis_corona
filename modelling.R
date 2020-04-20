#### Setup ####
rm(list=ls())

library(readxl, quietly=TRUE)
library(tidyverse, quietly=TRUE)
library(splm, quietly=TRUE)
library(plm, quietly=TRUE)

# Import standard variables
source("config.R")

# Read in metadata
df_meta = readxl::read_excel(path_wiki, sheet = "Metadata")

# Run `clean_full.R` to create the file at `path_full_long` imported on the
# next line.
df_long = readr::read_csv(path_full_long, col_types = do.call(
  cols, list(Date = col_date(format = "%Y-%m-%d"))))

# Weighting matrix: distance between the largest cities - currently unused
# if (file.exists(path_distances)) {
#   load(path_distances)
# } else {
#   df_distance = readxl::read_excel(path_wiki, sheet = "Distances") %>%
#     left_join(df_meta[, c("Code", "LargestCity")],
#               by=c("Distance (km)"="LargestCity")) %>%
#     column_to_rownames("Code") %>%
#     select(-"Distance (km)") %>%
#     data.matrix()
#   colnames(df_distance) = rownames(df_distance)
#   save(df_distance, file = path_distances)
# }

# Regressors: eurostat data, which are assumed to be time-constant
# To merge all Eurostat zip files, uncomment the next line if the file does not
# yet exist or if new data gets added.
# py_run_file("eurostat_reader.py")

df_eurostat = readr::read_csv(path_full_eurostat, col_types = do.call(
  cols, list(region=col_character())))

eurostat_variables = c("air_passengers_arrived", "air_passengers_departed",
                       "tourist_arrivals", "broadband_access",
                       "death_rate_diabetes", "death_rate_influenza",
                       "death_rate_chd", "death_rate_cancer", 
                       "death_rate_pneumonia", "available_beds",
                       "maritime_passengers_disembarked",
                       "maritime_passengers_embarked",
                       "risk_of_poverty_or_social_exclusion")

# Only keep rows where the `region` is an Italian region, not a direction or the
# entire country.
df_eurostat = df_eurostat[sapply(df_eurostat$region,
                                 function(x){x %in% df_meta$Region}), ] %>%
  select(c("region", "population_density", all_of(eurostat_variables))) %>%
  right_join(df_meta %>% select(c("Region", "Code")), by=c("region"="Region"))

# Impute columns based on nearest neighbors tourists or population
na_cols = sapply(df_eurostat, function(x){x %>% is.na() %>% any()}) %>%
  .[.] %>%
  names()

num_neighbors = 1 # Amount of neighbors to take into account when imputing

for (na_col in na_cols) {
  if (grepl("passenger", na_col, fixed=TRUE)) {
    diff_col = "tourist_arrivals"
  } else {
    diff_col = "population_density"
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
iddat = expand.grid(Date = unique(df_long$Date), Code = unique(df_long$Code))
iddat = iddat[order(iddat$Date, iddat$Code), ]
rownames(iddat) <- NULL
df_eurostat_panel = left_join(iddat, df_eurostat, by = "Code") %>%
  as_tibble()

# Import railway travellers data. This is interpolated from the Google Mobility
# Report (by eye and hand).
df_rail_travel = readr::read_csv(path_interpolated_rail, col_types = do.call(
  cols, list(Date = col_date(format = "%Y-%m-%d")))) %>%
  pivot_longer(cols = -Date, names_to = "Code", values_to = "rail_travelers")

# Join the expanded eurostat data and the railway data with the long data.
df_long = df_long %>%
  left_join(df_eurostat_panel, by = c("Date", "Code")) %>%
  left_join(df_rail_travel, by = c("Date", "Code"))

# TODO: Run models
# TODO: Add railway too
# Create formula to be the product of the regressors with the lagged incidence
# and susceptible rates
fm = paste("incidenceRate ~",
           "lag(incidenceRate, -1):lag(susceptibleRate, -1):" %>%
  paste0(eurostat_variables) %>%
  paste(collapse="+")) %>%
  as.formula()

plmo_pooled = plm(fm, data = df_long, model = "pooling")
plmo_fe = plm(fm, data = df_long, model = "within")
pooltest(plmo_pooled, plmo_fe) # p-value < 2.2e-16 => H1: unstability => FE
