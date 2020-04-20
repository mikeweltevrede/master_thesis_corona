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

# Weighting matrix: distance between the largest cities
if (file.exists(path_distances)) {
  load(path_distances)
} else {
  df_distance = readxl::read_excel(path_wiki, sheet = "Distances") %>%
    left_join(df_meta[, c("Code", "LargestCity")],
              by=c("Distance (km)"="LargestCity")) %>%
    column_to_rownames("Code") %>%
    select(-"Distance (km)") %>%
    data.matrix()
  colnames(df_distance) = rownames(df_distance)
  save(df_distance, file = path_distances)
}

# Regressors: eurostat data, which are assumed to be time-constant
# To merge all Eurostat zip files, uncomment the next line if the file does not
# yet exist or if new data gets added.
# py_run_file("eurostat_reader.py")

df_eurostat = readr::read_csv(path_full_eurostat, col_types = do.call(
  cols, list(region=col_character())))

# Only keep rows where the `region` is an Italian region, not a direction or the
# entire country.
df_eurostat = df_eurostat[sapply(df_eurostat$region,
                                 function(x){x %in% df_meta$Region}), ]

eurostat_variables = c("air_passengers_arrived", "air_passengers_departed",
                       "tourist_arrivals", "broadband_access",
                       "death_rate_diabetes", "death_rate_influenza",
                       "death_rate_chd", "death_rate_cancer", 
                       "death_rate_pneumonia", "available_beds",
                       "maritime_passengers_disembarked",
                       "maritime_passengers_embarked",
                       "risk_of_poverty_or_social_exclusion")
X = df_eurostat[, c("region", "population_density", eurostat_variables)]

# Impute columns based on nearest neighbors tourists or population
na_cols = sapply(X, function(x){x %>% is.na() %>% any()}) %>%
  .[.] %>%
  names()

for (na_col in na_cols) {
  if (grepl("passenger", na_col, fixed=TRUE)) {
    diff_col = "tourist_arrivals"
  } else {
    diff_col = "population_density"
  }
  
  na_check = X[[na_col]] %>% is.na()
  
  # Loop over regions/rows that are NA
  for (r in which(na_check)) {
    imputation = vector()
    
    X_full = X[!na_check, ]
    x = abs(X_full[[diff_col]] - X[[diff_col]][w])
    mins = sort(x)[1:2]
    
    for (m in mins) {
      ratio = m / X[[diff_col]][w]
      imputation = c(imputation,
                     X_full[[na_col]][which(x == m)] *
                       ratio)
    }
    
    X[r, na_col] = mean(imputation)
  }
}

# TODO: Process the Eurostat data so that we can use these as regressors.
df_rail_travel = readr::read_csv(path_interpolated_rail, col_types = do.call(
  cols, list(Date = col_date(format = "%Y-%m-%d"))))

# Select the region we are interested in and replicate the data to have the same
# amount of rows as the time series data
# region = "LOM"
# region_full = df_meta$Region[which(df_meta$Code == region)]
# df_eurostat_region = df_eurostat[which(df_eurostat$region == region_full), ] %>%
#   select(-c("region")) %>%
#   slice(rep(1:n(), each = nrow(df_wide)))