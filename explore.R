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

df_eurostat = readr::read_csv(path_full_eurostat, col_types = do.call(
  cols, list(region=col_character()))) %>%
  right_join(df_meta %>% select(c("Region", "Code")), by=c("region"="Region"))

dischargeRates = df_eurostat %>% select(starts_with("dischargeRate"))

corrs = correlate(dischargeRates)
corrs %>% rplot(print_cor = TRUE)
corrs %>% xtable %>% print(booktabs=TRUE)