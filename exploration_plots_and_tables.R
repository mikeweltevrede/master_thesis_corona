#### Setup ####
rm(list=ls())

library(GGally)
library(ggplot2)
library(glue)
library(latex2exp)
library(readxl)
library(reticulate)
library(tidyverse)

# Import standard variables
source("config.R")

# Read in metadata
df_meta = readxl::read_xlsx(path_metadata, sheet = "Metadata")

df_long = readr::read_csv(path_full_long, col_types = do.call(
  cols, list(date = col_date(format = "%Y-%m-%d")))) %>%
  left_join(df_meta %>% select(code, regionGH, direction), by="code") %>%
  mutate(weekday = date %>% lubridate::wday(label = TRUE) %>% as.factor)

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
df_long = df_long %>% 
  mutate(proportionTested = totalTested / totalPopulation,
         proportionUndocumentedQuadratic = 1 - proportionDocumentedQuadratic)

testedColor = "#0072B2"
propColor = "#D55E00"

# First, with the total number of tests executed...
# For a secondary axis, we need a scaling factor; the 2 is hardcoded to make the
# axis range from 0 to 1
coeff = mean(df_long$totalTested / df_long$proportionUndocumentedQuadratic) * 6

ggplot(df_long, aes(x=date)) +
  geom_line(aes(y=totalTested), size=0.8, color=testedColor) + 
  facet_wrap(vars(code)) +
  geom_line(aes(y=proportionUndocumentedQuadratic * coeff), size=0.8,
            color = propColor) +
  xlab("") +
  scale_y_continuous(name = "Total number of tests executed\n",
                     sec.axis = sec_axis(
                       ~./coeff, name="Proportion of infectives undocumented\n")) +
  theme(
    axis.text.y = element_text(color = testedColor),
    axis.text.y.right = element_text(color = propColor),
    axis.title.y = element_text(color = testedColor, size=12),
    axis.title.y.right = element_text(color = propColor, size=12),
    panel.spacing = unit(0.8, "lines"))

ggsave("tamponi_vs_ft.pdf", path=output_path)

# Then for the quotient of tests to the total population
coeff = mean(df_long$proportionTested / df_long$proportionUndocumentedQuadratic) * 3.6

ggplot(df_long, aes(x=date)) +
  geom_line(aes(y=proportionTested), size=0.8, color=testedColor) + 
  facet_wrap(vars(code)) +
  geom_line(aes(y=proportionUndocumentedQuadratic * coeff), size=0.8,
            color = propColor) + 
  xlab("") +
  scale_y_continuous(name = "Tests executed over total population\n",
                     sec.axis = sec_axis(
                       ~./coeff, name="Proportion of infectives undocumented\n")) +
  theme(
        axis.text.y = element_text(color = testedColor),
        axis.text.y.right = element_text(color = propColor),
        axis.title.y = element_text(color = testedColor, size=12),
        axis.title.y.right = element_text(color = propColor, size=12),
        panel.spacing = unit(0.8, "lines"))

ggsave("tamponiprop_vs_ft.pdf", path=output_path)

#### Numbers of documented infections table ####
data = df_long %>%
  filter(code %in% c("CAL", "LOM", "VEN")) %>%
  group_by(code) %>%
  mutate(totalInfectives = cumsum(infectives),
         totalInfectivesQuadratic = cumsum(infectivesQuadratic)) %>%
  ungroup() %>%
  select(date, code, totalInfectives, proportionDocumentedQuadratic,
         totalInfectivesQuadratic)

for (date in c("2020-04-01", "2020-06-01", "2020-08-01")) {
  data %>%
    filter(date == !!date) %>%
    print()
}

#### Plot heterogeneity ####
# Can we pool the data or not? If the means are similar enough, then the
# regions are similar and there is no region-specific effect apparent
df_long %>%
  group_by(regionGH) %>%
  summarize(mean = mean(infectivesRate),
            sd = sd(infectivesRate),
            se = sd / sqrt(n()),
            direction = direction,
            .groups = "drop") %>%
  distinct %>%
  arrange(direction, regionGH) %>%
  ggplot(aes(x=regionGH, y=mean, group=direction)) +
  xlab("Regions (grouped by NUTS 1 region)") +
  ylab(TeX("Mean infectives rate ($\\pm 2 \\times$ standard error)")) +
  geom_point(aes(colour = factor(direction)), stat = "identity",
             position = position_dodge(width = 1), size=3) +
  geom_path(aes(colour = factor(direction)), alpha=0.3) +
  geom_errorbar(aes(ymin=mean-2*se, ymax=mean+2*se,
                    color=direction), width=.3) +
  facet_grid(.~ direction, scales = "free", switch = "x", space = "free_x") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        strip.placement = "outside", legend.position = "none") +
  scale_colour_manual(values=c("#0072B2", # Dark blue
                               "#D55E00", # Orange-brown
                               "#CC79A7", # Pink
                               "#009E73", # Green
                               "#56B4E9", # Light blue
                               "#E69F00")) # Yellow

ggsave("heterogeneity_over_regions.pdf", path=output_path)

#### Infective rate per NUTS 1 region ####
df_long %>%
  group_by(direction, date) %>%
  mutate(infectiveRate = sum(infectives) / sum(totalPopulation),
         direction = as.factor(direction)) %>%
  ungroup() %>%
  select(date, direction, infectiveRate) %>%
  distinct() %>%
  ggplot(aes(x=date, y=infectiveRate, fill=direction)) +
  geom_bar(stat="identity") +
  facet_wrap(vars(direction)) +
  xlab("") +
  ylab("") +
  theme(legend.position="none") +
  scale_fill_manual(values=c("#0072B2", # Dark blue
                             "#D55E00", # Orange-brown
                             "#CC79A7", # Pink
                             "#009E73", # Green
                             "#56B4E9", # Light blue
                             "#E69F00")) # Yellow

ggsave("infective_rate_per_NUTS1.pdf", path=output_path)

#### Infective rate per NUTS 2 region for Nord-Est ####
for (direc in unique(df_long$direction)) {
  df_long %>%
    filter(direction == !!direc) %>%
    mutate(infectiveRate = infectives / totalPopulation,
           regionGH = as.factor(regionGH)) %>%
    select(date, regionGH, infectiveRate) %>%
    distinct() %>%
    ggplot(aes(x=date, y=infectiveRate, fill=regionGH)) +
    geom_bar(stat="identity") +
    facet_wrap(vars(regionGH)) +
    xlab("") +
    ylab("") +
    theme(legend.position="none") +
    scale_fill_manual(values=c("#0072B2", # Dark blue
                               "#D55E00", # Orange-brown
                               "#CC79A7", # Pink
                               "#009E73", # Green
                               "#56B4E9", # Light blue
                               "#E69F00")) # Yellow
  
  ggsave(glue("infective_rates_{direc}.pdf"), path=output_path)
}

#### Plot per day ####
for (direc in unique(df_long$direction)) {
  df_long %>%
    filter(direction == !!direc) %>%
    group_by(weekday, regionGH) %>%
    summarize(mean = mean(infectivesRate),
              sd = sd(infectivesRate),
              se = sd / sqrt(n()),
              weekday = weekday,
              regionGH = regionGH,
              .groups = "drop") %>%
    distinct %>%
    arrange(weekday) %>%
    ggplot(aes(x=weekday, y=mean)) +
    xlab("Day of the week") +
    ylab(TeX("Mean infectives rate ($\\pm 2 \\times$ standard error)")) +
    geom_point(color="#0072B2", stat = "identity",
               position = position_dodge(width = 1), size=3) +
    geom_errorbar(aes(ymin=mean-2*se, ymax=mean+2*se,
                      color="#0072B2"), width=.3) +
    facet_wrap(vars(regionGH)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          strip.placement = "outside", legend.position = "none") +
    scale_colour_manual(values=c("#0072B2")) # Dark Blue
  
  ggsave(glue("infective_rates_weekday_{direc}.pdf"), path=output_path)
}
