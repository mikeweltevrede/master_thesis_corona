# Environment variables
env_name = "r-thesis_corona"
path_requirements = "requirements.txt"

# Activate the Conda environment. If it does not exist yet, create it with the
# required packages.
tryCatch(use_condaenv(env_name),
         error=function(e){
           requirements = scan(file=path_requirements, what=character(),
                               quiet=TRUE) 
           conda_create(env_name, packages=requirements)
         })

# Path variables for data
data_path = "data"
path_wiki = paste0(data_path, "/italy_wikipedia.xlsx")
path_cleaned_wide = paste0(data_path, "/italy_wikipedia_cleaned.csv")
path_full_wide = paste0(data_path, "/italy_wikipedia_full_wide.csv")
path_full_long = paste0(data_path, "/italy_wikipedia_full_long.csv")
path_full_eurostat = paste0(data_path, "/merged_eurostat.csv")
path_mobility_report = paste0(data_path, "/google_mobility_report.xlsx")
path_railroad = paste0(data_path, "/eurostat/interregion_railroad_travel.csv")
path_interpolated_rail = paste0(data_path, "/interpolated_railroad_travel.csv")
