#### Unused code - keep for future ####
# Weighting matrix: distance between the largest cities - currently unused
if (file.exists(path_distances)) {
  load(path_distances)
} else {
  df_distance = readxl::read_xlsx(path_wiki, sheet = "Distances") %>%
    left_join(df_meta[, c("Code", "LargestCity")],
              by=c("Distance (km)"="LargestCity")) %>%
    column_to_rownames("Code") %>%
    select(-"Distance (km)") %>%
    data.matrix()
  colnames(df_distance) = rownames(df_distance)
  save(df_distance, file = path_distances)
}

