#### Setup ####
rm(list=ls())

library(readxl)
library(tidyverse)
library(reticulate)
library(GGally)

# Import standard variables
source("config.R")

# Read in metadata
df_meta = readxl::read_xlsx(path_metadata, sheet = "Metadata")

#### Make correlation matrix of discharge rates  ####
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
ggsave("correlations_discharge_rates.pdf", path=output_path)

#### Predict values ####
# # Predict values -> Odd results with negative values and a completely different
# # curve
# fit_my_model = function(t, data) {
#   model = lm(fm, data=data[1:t, ])
#   return(predict(model, newdata=data[1:(t+1), ])[[t+1]])
# }
# 
# start = 20
# end = nrow(df_wide)-1
# infectives = sapply(start:end, fit_my_model, data=df_wide)
# plot(df_wide$infectivesTotal[start:end], type="l")
# lines(infectives, col="red")
