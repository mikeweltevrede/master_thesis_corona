#### Setup ####
rm(list=ls())

library(readxl, quietly=TRUE)
library(tidyverse, quietly=TRUE)
library(reticulate, quietly=TRUE)

# Import standard variables
source("config.R")

# We need to clean the Wikipedia data before being able to process it in R, also
# to include new dates. Run the next line to do so (you may need to
# install Miniconda as a Python interpreter).
py_run_file("clean_wide.py")

# Read in the cleaned Wikipedia data
df_wide = read_csv(path_cleaned_wide, col_types = do.call(
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

# Complete the dates and execute the filling procedure described above
df_wide = df_wide %>%
  complete(Date = all_dates) %>%
  replace_na(named_cols_replace_0) %>%
  fill(all_of(cols_fill))

# Convert the data to long and join it with the extra data (containing the
# amount of active ICU patients, recoveries, tested people, and positively
# tested people)
df_wide = df_wide %>%
  select(Date:SAR_Deaths)

df_long = df_wide %>%
  pivot_longer(cols = -Date,
               names_to = c("group", ".value"),
               names_sep = "_")

df_extra = read_excel(path_wiki, sheet="Extra") %>%
  mutate(Date = as.Date(Date, format = "%Y-%m-%d")) %>%
  drop_na()

df_extra_long = df_extra %>% pivot_longer(cols = -Date,
                                          names_to = c("group", ".value"),
                                          names_sep = "_")

df_wide_full = full_join(df_wide, df_extra, by="Date")
df_long_full = full_join(df_long, df_extra_long, by=c("Date", "group"))

# Note that, for example, cumsum(ABR_Confirmed) != ABR_TestedPositive. However,
# we do assume that the missing values, i.e. those before March 2, are correctly
# imputed by the cumsum of the confirmed cases, as long as the final element in
# the cumsum is less than the first non-missing element of TestedPositive.
# Unfortunately, we cannot impute _Tested.
df_meta = read_excel(path_wiki, sheet = "Metadata")

for (region in df_meta$Code){
  # For completeness sake, even though the NAs (should) align, we find them per
  # region
  which_NA = df_wide_full %>%
    select(!!paste0(region, "_TestedPositive")) %>%
    is.na() %>%
    which()
  
  csum = df_wide_full %>%
    select(!!paste0(region, "_Confirmed")) %>%
    cumsum() %>%
    .[which_NA, ]
  
  first_elt = df_wide_full %>%
    select(!!paste0(region, "_TestedPositive")) %>%
    na.omit() %>%
    first()
  
  # Check if the final element in the cumsum is lower than
  if (tail(csum, 1) <= first_elt) {
    df_wide_full[which_NA, paste0(region, "_TestedPositive")] = csum
  } else {
    print(paste0("For region ", region, ", the final element of csum (",
                 tail(csum, 1), ") was not lower than the first NA element (",
                 first_elt, "). We cannot impute values for this region."))
  }
}

# Save the tibbles to a file
write_csv(df_wide_full, path_full_wide)
write_csv(df_long_full, path_full_long)
