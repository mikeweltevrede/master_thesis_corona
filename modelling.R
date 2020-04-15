#### Setup ####
rm(list=ls())

library(readxl, quietly=TRUE)
library(tidyverse, quietly=TRUE)
library(splm, quietly=TRUE)

# Import standard variables
source("config.R")

# Read in metadata
df_meta = read_excel(path_wiki, sheet = "Metadata")

# Run `clean_full.R` to create the file at `path_full_long` imported on the
# next line.
df_long = read_csv(path_full_long, col_types = do.call(
  cols, list(Date = col_date(format = "%Y-%m-%d"))))

# Weighting matrix: distance between the largest cities
df_distance = read_excel(path_wiki, sheet = "Distances") %>%
  left_join(df_meta[, c("Code", "LargestCity")],
            by=c("Distance (km)"="LargestCity")) %>%
  column_to_rownames("Code") %>%
  select(-"Distance (km)") %>%
  data.matrix()
colnames(df_distance) = rownames(df_distance)

# Regressors: eurostat data, which are assumed to be time-constant
# To merge all Eurostat zip files, uncomment the next line if the file does not
# yet exist.
# py_run_file("eurostat_reader.py")

df_eurostat = read_csv(path_full_eurostat, col_types = do.call(
  cols, list(region=col_character())))

# Only keep rows where the `region` is an Italian region, not a direction or the
# entire country.
df_eurostat = df_eurostat[sapply(df_eurostat$region,
                                 function(x){x %in% df_meta$Region}), ]

# TODO: Process the Eurostat data so that we can use these as regressors.
df_rail_travel = read_csv(path_interpolated_rail, col_types = do.call(
  cols, list(Date = col_date(format = "%Y-%m-%d"))))

# Select the region we are interested in and replicate the data to have the same
# amount of rows as the time series data
# region = "LOM"
# region_full = df_meta$Region[which(df_meta$Code == region)]
# df_eurostat_region = df_eurostat[which(df_eurostat$region == region_full), ] %>%
#   select(-c("region")) %>%
#   slice(rep(1:n(), each = nrow(df_wide)))