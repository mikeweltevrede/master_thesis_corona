#### Setup ####
rm(list=ls())

library(readxl)
library(tidyverse)
library(reticulate)
library(GGally)

# Import standard variables
source("config.R")

# Read in metadata
df_meta = readxl::read_xlsx(path_metadata, sheet = "Metadata")

#### Make correlation matrix of discharge rates  ####
# To merge all Eurostat zip files, uncomment the next line if the file does not
# yet exist or if new data gets added.
# reticulate::py_run_file("eurostat_reader.py")

dischargeRates = readr::read_csv(path_full_eurostat, col_types = do.call(
  cols, list(region=col_character()))) %>%
  right_join(df_meta %>%
               select(c("Region", "Code")), by=c("region"="Region")) %>%
  select(starts_with("dischargeRate")) %>%
  rename_with(function(x){str_replace(x, "dischargeRate", "")}) %>%
  select(-Respiratory) # TODO: Remove respiratory diseases from Eurostat

GGally::ggcorr(dischargeRates, label=TRUE, hjust = 0.75, size = 4.5,
               layout.exp = 1, label_round = 2)
ggsave("correlations_discharge_rates.pdf", path=output_path)

#### Predict values ####
# # Predict values -> Odd results with negative values and a completely different
# # curve
# fit_my_model = function(t, data) {
#   model = lm(fm, data=data[1:t, ])
#   return(predict(model, newdata=data[1:(t+1), ])[[t+1]])
# }
# 
# start = 20
# end = nrow(df_wide)-1
# infectives = sapply(start:end, fit_my_model, data=df_wide)
# plot(df_wide$infectivesTotal[start:end], type="l")
# lines(infectives, col="red")

#### Plot undocumented infections factor over time ####
df_long = readr::read_csv(path_full_long, col_types = do.call(
  cols, list(date = col_date(format = "%Y-%m-%d")))) %>%
  mutate(proportionTested = totalTested / totalPopulation)

testedColor = "#0072B2"
propColor = "#D55E00"

# Then we can make the plots
# First, with the total number of tests executed...
# For a secondary axis, we need a scaling factor; the 2 is hardcoded to make the
# axis range from 0 to 1
coeff = mean(df_long$totalTested / df_long$proportionDocumentedQuadratic) * 2

ggplot(df_long, aes(x=date)) +
  geom_line(aes(y=totalTested), size=0.8, color=testedColor) + 
  facet_wrap(vars(code)) +
  geom_line(aes(y=proportionDocumentedQuadratic * coeff), size=0.8,
            color = propColor) + 
  xlab("") +
  scale_y_continuous(name = "Total number of tests executed\n",
                     sec.axis = sec_axis(
                       ~./coeff, name="Proportion of infectives documented\n")) +
  theme(
    axis.text.y = element_text(color = testedColor),
    axis.text.y.right = element_text(color = propColor),
    axis.title.y = element_text(color = testedColor, size=12),
    axis.title.y.right = element_text(color = propColor, size=12),
    panel.spacing = unit(0.8, "lines"))

ggsave("tamponi_vs_ft.pdf", path=output_path)

# Then for the quotient of tests to the total population
coeff = mean(df_long$proportionTested / df_long$proportionDocumentedQuadratic) * 1.2

ggplot(df_long, aes(x=date)) +
  geom_line(aes(y=proportionTested), size=0.8, color=testedColor) + 
  facet_wrap(vars(code)) +
  geom_line(aes(y=proportionDocumentedQuadratic * coeff), size=0.8,
            color = propColor) + 
  xlab("") +
  scale_y_continuous(name = "Tests executed over total population\n",
                     sec.axis = sec_axis(
                       ~./coeff, name="Proportion of infectives documented\n")) +
  theme(
        axis.text.y = element_text(color = testedColor),
        axis.text.y.right = element_text(color = propColor),
        axis.title.y = element_text(color = testedColor, size=12),
        axis.title.y.right = element_text(color = propColor, size=12),
        panel.spacing = unit(0.8, "lines"))

ggsave("tamponiprop_vs_ft.pdf", path=output_path)
