library(reticulate)
library(glue)

# Environment variables
env_name = "r-thesis_corona"
path_requirements = "requirements.txt"

# Activate the Conda environment. If it does not exist yet, create it with the
# required packages.
tryCatch(
  use_condaenv(env_name),
  error = function(e){
    conda_create(env_name)
    conda_install(env_name, packages=scan(file=path_requirements,
                                          what=character(), quiet=TRUE))
    use_condaenv(env_name)
  }
)

# Path variables for data
data_path = "data"

# For data import
data_path_gh = glue("{data_path}/COVID-19/dati-regioni")
completed_dates_path = glue("{data_path}/completed_dates.txt")
new_data_path = glue("{data_path}/data_long.csv")
new_data_path_wide = glue("{data_path}/data_wide.csv")
new_data_path_wide_cleaned = glue("{data_path}/data_wide_cleaned.csv")
new_data_path_long_cleaned = glue("{data_path}/data_long_cleaned.csv")
path_full_long = glue("{data_path}/data_long_full.csv")
path_full_wide = glue("{data_path}/data_wide_full.csv")

path_metadata = glue("{data_path}/metadata.xlsx")
path_full_eurostat = glue("{data_path}/merged_eurostat.csv")
path_mobility_report_official = glue("{data_path}/Global_Mobility_Report.csv")
path_mobility_report_cleaned = glue("{data_path}/Global_Mobility_Report_Cleaned.csv")
path_railroad = glue("{data_path}/eurostat/interregion_railroad_travel.csv")

# For images and such
output_path = "output"
