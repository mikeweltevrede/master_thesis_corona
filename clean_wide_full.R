#### Setup ####
rm(list=ls())

library(readxl, quietly=TRUE)
library(tidyverse, quietly=TRUE)
library(reticulate, quietly=TRUE)

data_path = "data"

# Activate the Conda environment. If it does not exist yet, create it with the
# required packages.
env_name = "r-thesis_corona"

tryCatch(use_condaenv(env_name),
         error=function(e){
           requirements = scan(file="requirements.txt", what=character(),
                               quiet=TRUE) 
           conda_create(env_name, packages=requirements)
         })

# We need to clean the Wikipedia data before being able to process it in R, also
# to include new dates. Run the next line to do so (you may need to
# install Miniconda as a Python interpreter).
py_run_file("clean_wide.py")

# Read in the cleaned Wikipedia data
df_wide = read_csv(paste0(data_path, "/italy_wikipedia_cleaned.csv"),
                   col_types = do.call(
                     cols, list(Date=col_date(format="%Y-%m-%d"))))

# We now add missing dates to the data for equal spacing.
all_dates = seq.Date(min(df_wide$Date), max(df_wide$Date), by="day")
missing_dates = all_dates[which(!all_dates %in% df_wide$Date)][-1]
sprintf("The following %d dates are missing and will be added: %s",
        length(missing_dates), paste(missing_dates, collapse=", "))

# For the non-aggregated variables, we enter 0. For the totals, we use the
# previous value.
cols_replace_0 = c(colnames(df_wide)[
  2:(which(colnames(df_wide) == "Confirmed_New")-1)], "Confirmed_New",
  "Deaths_New")
named_cols_replace_0 = set_names(as.list(rep(0, length(cols_replace_0))),
                                 cols_replace_0)

# The columns to backfill are the remaining columns minus "Date".
cols_fill = colnames(df_wide)[
  which(!colnames(df_wide) %in% cols_replace_0)][-1]

df_wide = df_wide %>%
  complete(Date = all_dates) %>%
  replace_na(named_cols_replace_0) %>%
  fill(all_of(cols_fill))

# Save the tibble to a file
write_csv(df_wide, paste0(data_path, "/italy_wikipedia_full_wide.csv"))

# Convert the data to wide and join it with the extra data (containing the
# amount of active ICU patients, recoveries, tested people, and positively
# tested people)
df_long = df_wide %>%
  select(Date:SAR_Deaths) %>%
  pivot_longer(cols = -Date,
               names_to = c("group", ".value"),
               names_sep = "_")

df_extra_long = read_excel(
  paste0(data_path, "/italy_wikipedia.xlsx"), sheet="Extra") %>%
  mutate(Date = as.Date(Date, format = "%Y-%m-%d")) %>%
  drop_na() %>%
  pivot_longer(cols = -Date,
               names_to = c("group", ".value"),
               names_sep = "_")

df_long_full = full_join(df_long, df_extra_long, by=c("Date", "group"))

write_csv(df_long_full, paste0(data_path, "/italy_wikipedia_full_long.csv"))
