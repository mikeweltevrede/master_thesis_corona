eurostat_reader = function(zip_path){
  require(tidyverse)
  require(magrittr)
  
  df = data.frame()
  files = unzip(zip_path, list=TRUE)
  
  for (j in 1:nrow(files)) {
    # Only keep files that have Data appended
    if (grepl("Data", files[j, "Name"])){
      df = rbind(df, readr::read_csv(unz(zip_path, files[j, "Name"]),
                                     col_names=TRUE))
    }
  }
  
  ## Currently hardcoded: remove certain constant columns
  cols_to_drop = c("UNIT", "C_RESID", "AGE")
  if (length(intersect(colnames(df), cols_to_drop)) > 0){
    df = df %>%
      select(-intersect(colnames(df), cols_to_drop))
  }
  
  # Clean - because missing values are indicated with a colon instead of a blank
  # space, the Value column is of type character.
  df = df %>%
    mutate(Value=replace(Value, Value==":", NA)) %>%
    mutate(Value = str_replace(Value, " ", "")) %>% # Numbers contain spaces or commas, inducing NAs by conversion
    mutate(Value = str_replace(Value, ",", "")) %<>% # Idem
    mutate(Value = as.double(Value))
  
  # Find the column containing the relevant values to `spread` on
  value_col = colnames(df)[!colnames(df) %in% c("TIME", "GEO", "Value")]
  
  # If there are 2 or more matches, stop the code: we won't know how to `spread`
  ## If we have no matches, then we can simply keep `Value` as our data column
  stopifnot((length(value_col)==1 | length(value_col)==0))
  
  # Either spread on the columns or rename the Value column to be informative
  if (length(value_col) > 0) {
    df = spread(df, value_col, Value)
  } else {
    df = rename(df, !!basename(zip_path) := Value)
  }
  
  # Drop all columns that consist only of missing values
  df = select_if(df, function(x){!all(is.na(x))})
  
  return(df)
}

# Find all filenames in `data_path` ending in `.zip` - One could consider making
# a separate folder "data/eurostat" and changing `data_path` accordingly
data_path = "data"
filelist = paste0("data/", list.files(data_path, pattern=".zip"))

# Initialise full data.frame with the first file
df_full = eurostat_reader(filelist[1])

# DO an outer join with the next data
for (file in filelist[2:length(filelist)]){
  
  df = eurostat_reader(file)
  
  # TODO: Check that the columns in df do not already appear in df_full
  
  df_full = dplyr::full_join(df_full, df, by = c("TIME", "GEO"))
}

df_full