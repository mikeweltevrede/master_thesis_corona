source("config.R")

library(incidence)
library(EpiEstim)
library(ggplot2)

df_long = readr::read_csv(path_full_long, col_types = do.call(
  cols, list(date = col_date(format = "%Y-%m-%d"))))

region = "LOM"

df_incidence = df_long %>%
  filter(code == !!region) %>%
  transmute(dates = date, I = infectives)

# SI: serial interval (time interval between symptoms onset in a case and in
# their infector)
# Du et al: https://wwwnc.cdc.gov/eid/article/26/6/20-0357_article
mean_si = 3.96 # Du et al
std_si = 4.75 # Du et al
res = estimate_R(df_incidence, method="parametric_si",
                 config = make_config(list(mean_si = mean_si, std_si = std_si)))
plot(res, legend = FALSE)
