#### Setup ####
rm(list=ls())

library(readxl, quietly=TRUE)
library(tidyverse, quietly=TRUE)
library(lubridate, quietly=TRUE)
library(reticulate, quietly=TRUE)
library(corrr)

# Import standard variables
source("config.R")

# Read in metadata
df_meta = readxl::read_xlsx(path_wiki, sheet = "Metadata")

# To merge all Eurostat zip files, uncomment the next line if the file does not
# yet exist or if new data gets added.
# reticulate::py_run_file("eurostat_reader.py")

# We know the amount of people on January 1, 2019 as defined in df_eurostat. We
# only keep rows where the `region` is an Italian region, not a direction/NUTS-1
# region or the entire country by right joining with df_meta, since df_meta only
# contains NUTS-2 regions.

df_eurostat = readr::read_csv(path_full_eurostat, col_types = do.call(
  cols, list(region=col_character()))) %>%
  right_join(df_meta %>% select(c("Region", "Code")), by=c("region"="Region"))

dischargeRates = df_eurostat %>% select(starts_with("dischargeRate"))

corrs = correlate(dischargeRates)
corrs %>% rplot(print_cor = TRUE)
corrs %>% xtable %>% print(booktabs=TRUE)