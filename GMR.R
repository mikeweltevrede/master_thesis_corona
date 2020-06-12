#### Setup ####
# Import standard variables
source("config.R")

library(tidyverse)
library(snakecase)

# Read in metadata
df_meta = readxl::read_xlsx(path_wiki, sheet = "Metadata")

df_gmr = read_csv(path_mobility_report_official,
                  col_types = do.call(
                    cols,
                    list(
                      country_region_code = col_character(),
                      country_region = col_character(),
                      sub_region_1 = col_character(),
                      sub_region_2 = col_character(),
                      date = col_date(format = "%Y-%m-%d"),
                      retail_and_recreation_percent_change_from_baseline = col_double(),
                      grocery_and_pharmacy_percent_change_from_baseline = col_double(),
                      parks_percent_change_from_baseline = col_double(),
                      transit_stations_percent_change_from_baseline = col_double(),
                      workplaces_percent_change_from_baseline = col_double(),
                      residential_percent_change_from_baseline = col_double()))) %>%
  filter(country_region_code == "IT") %>%
  
  # Drop unused columns
  select(-country_region_code, -country_region, -sub_region_2) %>%
  
  # Drop rows with NAs in column `sub_region_1` (region name)
  drop_na(any_of("sub_region_1")) %>%
  
  # Clean column names
  rename_at(vars(ends_with("_percent_change_from_baseline")),
            funs(str_replace(., "_percent_change_from_baseline", "")))

colnames(df_gmr) = colnames(df_gmr) %>%
  snakecase::to_upper_camel_case() %>%
  str_replace("SubRegion1", "Region")
  
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
# from Eurostat is from 2015. We assume that this value has not changed much.
df_rail = read_csv(path_railroad) %>%
  filter(TIME == max(TIME))

regions = df_gmr$Code %>% unique()
baselines = vector()

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
  
  num_rows = df_gmr[df_gmr$Code == region_code, ] %>% nrow
  baselines = c(baselines, rep(baseline, num_rows))
}

df_gmr = df_gmr %>%
  transmute(Date = Date, Code = Code,
            railTravelers = round(TransitStations*eval(baselines)))
