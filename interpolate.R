#### Setup ####
rm(list = ls())

library(readxl)
library(tidyverse)
library(lubridate)

source("config.R")

dfs_mobility = path_mobility_report %>% 
  excel_sheets() %>% 
  set_names() %>% 
  map(read_excel, path = path_mobility_report)

df_meta = readxl::read_xlsx(path_wiki, sheet = "Metadata")

#### Interpolate proportions ####
interpolate = function(data, metadata, max_date = NULL, only=NULL) {
  
  if (is.null(max_date)) {
    data = data %>%
      mutate(Date = as.Date(Date, format = "%d/%m/%Y")) %>%
      complete(Date = seq.Date(min(Date), max(Date), by = "day"))
  } else {
    max_date = as.Date(max_date, format = "%d/%m/%Y")
    data = data %>%
      mutate(Date = as.Date(Date, format = "%d/%m/%Y")) %>%
      complete(Date = seq.Date(min(Date), max_date, by = "day"))
  }
  
  if (is_null(only)){
    # Take all columns except for "Date"; else, take the ones specified in only
    only = colnames(data)[colnames(data) != "Date"]
  }
  
  # We now run the interpolation for each column
  for (region_code in only) {
    region_name = metadata[df_meta$Code == region_code, ][["Region"]]
    missing_dates = data[data[, region_code] %>% is.na(), ][["Date"]]
    
    for (date in missing_dates) {
      
      # Note that we have not specified a value for the final date in the data
      # since this is unknown. We assume that the value from March 29 stays
      # constant, rule=2 argument.
      data[data$Date == date, region_code] = approx(
        data$Date, data[[region_code]], xout = date,
        ties = "ordered", rule = 2)$y
    }
  }
  
  return(data)
}

# Use the function; save each interpolated sheet as a list element
dfs_interpolated = vector("list")

# We get the maximum date available in our cleaned wikipedia data
max_date = readr::read_csv(path_cleaned_wide, col_types = do.call(
  cols, list(Date=col_date(format="%Y-%m-%d")))) %>%
  .[["Date"]] %>%
  max()

for (data in names(dfs_mobility)[names(dfs_mobility) != "Overall"]) {
  dfs_interpolated[[data]] = interpolate(dfs_mobility[[data]], df_meta,
                                         max_date)
}

# For railroad transport, we can multiply by the baseline value
df_rail = read_csv(path_railroad) %>%
  filter(TIME == max(TIME))

regions = colnames(dfs_interpolated[["Transit stations"]])[
  colnames(dfs_interpolated[["Transit stations"]]) != "Date"]

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
  
  # The given numbers in df_rail are per year so we need daily numbers. This
  # depends on whether it is a leap year.
  if (df_rail$TIME %>% max() %>% leap_year()) {
    baseline = baseline / 366
  } else {
    baseline = baseline / 365
  }
  
  dfs_interpolated[["Transit stations"]][[region_code]] =
    dfs_interpolated[["Transit stations"]][[region_code]] * baseline
}

readr::write_csv(dfs_interpolated[["Transit stations"]], path_interpolated_rail)
