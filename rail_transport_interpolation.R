rm(list = ls())

library(lubridate)

#### Read in the data #####
# df_changes contains percentages of passengers on that day according to the
# Google Mobility Report of March 29, 2020.
df_changes = read_excel("data/google_transport_interpolation.xlsx") %>%
  mutate(Date = as.Date(Date, format = "%d/%m/%Y")) %>%
  complete(Date = seq.Date(min(Date), max(Date), by = "day"))

# We assume that the amount of passengers in 2015 is representative of the
# amount of passengers in 2020.
df_rail = read_csv("data/eurostat/interregion_railroad_travel.csv") %>%
  filter(TIME == max(TIME))

df_meta = read_excel("data/italy_wikipedia.xlsx",
                     sheet = "Metadata")

# Initialise tibble of passengers
df_passengers = df_changes

# Initialise empty list for baseline passengers
baseline = vector("list")

for (region_code in colnames(df_passengers)[-1]) {
  region_name = df_meta[df_meta$Code == region_code, ][["Region"]]

  # Define the baseline as the amount of passengers that arrive at and depart
  # from that region. The Google Mobility report defines this as the median
  # amount of passengers.
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
  
  missing_dates = df_passengers[df_passengers[, region_code] %>% is.na(),
                                "Date"]
  
  for (date in unlist(missing_dates)) {
    
    # Note that we have not specified a value for the final date in the data
    # since this is unknown. We assume that the value from March 29 stays
    # constant, rule=2 argument.
    
    df_passengers[df_passengers$Date == date, region_code] = approx(
      df_passengers$Date, df_passengers[[region_code]], xout = date,
      ties = "ordered", rule = 2)$y
  }
  
  # Having computed all percentages, we can multiply these with the baseline
  df_passengers[[region_code]] = df_passengers[[region_code]] * baseline
  
}

tail(df_passengers, 10)

# Note that this may be unrealistic. In total, we still have a high total amount
# of passengers in all of Italy at the final date, namely 298735.1:
df_passengers %>% select(-Date) %>% rowSums %>% tail(1)

# Two comments on this:
# 1. There are still some people who do have to use public transport to go to
# work, for instance essential workers like doctors and nurses.
# 2. We may also choose to set the percentage to be 0 at a certain date (e.g.
# the date of today) and to then let the amount of passengers decay linearly.
# For now, we neglect this.
