#### Setup ####
rm(list=ls())

library(readxl)
library(tidyverse)
library(reticulate)
library(GGally)

# Import standard variables
source("config.R")

# Read in metadata
df_meta = readxl::read_xlsx(path_wiki, sheet = "Metadata")

# To merge all Eurostat zip files, uncomment the next line if the file does not
# yet exist or if new data gets added.
# reticulate::py_run_file("eurostat_reader.py")

dischargeRates = readr::read_csv(path_full_eurostat, col_types = do.call(
  cols, list(region=col_character()))) %>%
  right_join(df_meta %>%
               select(c("Region", "Code")), by=c("region"="Region")) %>%
  select(starts_with("dischargeRate")) %>%
  rename_with(function(x){str_replace(x, "dischargeRate", "")}) %>%
  select(-Respiratory) # TODO: Remove respiratory diseases from Eurostat

GGally::ggcorr(dischargeRates, label=TRUE, hjust = 0.75, size = 4.5,
               layout.exp = 1, label_round = 2)
ggsave("correlations_discharge_rates.png", path=output_path, dpi=300)
