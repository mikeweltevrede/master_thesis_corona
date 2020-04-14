rm(list = ls())

library(readxl)
library(tidyverse)

file_path = "data/google_mobility_report.xlsx"

df = file_path %>% 
  excel_sheets() %>% 
  set_names() %>% 
  map(read_excel, path = file_path)

df_meta = read_excel("data/italy_wikipedia.xlsx",
                     sheet = "Metadata")

interpolate = function(data, metadata, only=NULL) {
  
  data = data %>%
    mutate(Date = as.Date(Date, format = "%d/%m/%Y")) %>%
    complete(Date = seq.Date(min(Date), max(Date), by = "day"))
  
  if (is_null(only)){
    # Take all columns except for "Date"; else, take the ones specified in only
    only = colnames(data)[colnames(data) != "Date"]
  }
  
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

# Use the function
df_interpolated = vector("list")

for (data in names(df)[-1]) {
  df_interpolated[[data]] = interpolate(df[[data]], df_meta, only="LOM")
}

df_interpolated

# For railway, we can multiply by the baseline value
df_rail = read_csv("data/eurostat/interregion_railroad_travel.csv") %>%
  filter(TIME == max(TIME))

regions = colnames(df_interpolated[["Transit stations"]])[
  colnames(df_interpolated[["Transit stations"]]) != "Date"]

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
  
  df_interpolated[["Transit stations"]][[region_code]] =
    df_interpolated[["Transit stations"]][[region_code]] * baseline
}


