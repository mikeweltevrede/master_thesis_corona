rm(list=ls())
# Environment variables
env_name = "r-thesis_corona"
path_requirements = "requirements.txt"

# Activate the Conda environment. If it does not exist yet, create it with the
# required packages.
library(reticulate)
if (env_name %in% conda_list()$name){
  use_condaenv(env_name)
} else {
  conda_create(env_name)
  conda_install(env_name, packages=scan(file=path_requirements,
                                        what=character(), quiet=TRUE))
  use_condaenv(env_name)
}

# Path variables for data
data_path = "data"
output_path = "output"
path_wiki = paste0(data_path, "/italy_wikipedia.xlsx")
path_cleaned_wide = paste0(data_path, "/italy_wikipedia_cleaned.csv")
path_full_wide = paste0(data_path, "/italy_wikipedia_full_wide.csv")
path_full_long = paste0(data_path, "/italy_wikipedia_full_long.csv")
path_full_eurostat = paste0(data_path, "/merged_eurostat.csv")
path_mobility_report = paste0(data_path, "/google_mobility_report.xlsx")
path_railroad = paste0(data_path, "/eurostat/interregion_railroad_travel.csv")
path_distances = paste0(data_path, "/distances.RData")
