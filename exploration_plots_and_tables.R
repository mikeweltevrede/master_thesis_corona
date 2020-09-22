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

df_wide = readr::read_csv(path_full_wide, col_types = do.call(
  cols, list(date = col_date(format = "%Y-%m-%d"))))

# Shade the lockdown
tbl_rects = tibble(
  xstart = head(df_wide$date[which(df_wide$lockdown == 1)], 1),
  xend = tail(df_wide$date[which(df_wide$lockdown == 1)], 1)
)

#### Plot undocumented infections factor over time ####
df_long = df_long %>% 
  mutate(proportionTested = testedTotal / totalPopulation,
         proportionUndocumentedQuadratic = 1 - proportionDocumentedQuadratic)

testedColor = "#0072B2"
propColor = "#D55E00"

# First, with the total number of tests executed...
# For a secondary axis, we need a scaling factor; the 2 is hardcoded to make the
# axis range from 0 to 1
coeff = mean(df_long$testedTotal / df_long$proportionUndocumentedQuadratic) * 6

g = ggplot(df_long) +
  geom_line(aes(x=date, y=testedTotal), size=0.8, color=testedColor) + 
  geom_rect(data = tbl_rects, aes(xmin = xstart, xmax = xend, ymin = -Inf, # Shadow over plot
                                  ymax = Inf), alpha = 0.2) +
  facet_wrap(vars(regionGH),
             nrow=unique(df_long$regionGH) %>% length %>% sqrt %>% floor) +
  geom_line(aes(x=date, y=proportionUndocumentedQuadratic * coeff), size=0.8,
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

ggsave("tamponi_vs_ft.pdf", plot = g, path=output_path, width = 10.8,
       height = 6.62, units = "in")

# Then for the quotient of tests to the total population
coeff = mean(df_long$proportionTested / df_long$proportionUndocumentedQuadratic) * 3.6

g = ggplot(df_long) +
  geom_line(aes(x=date, y=proportionTested), size=0.8, color=testedColor) + 
  geom_rect(data = tbl_rects, aes(xmin = xstart, xmax = xend, ymin = -Inf, # Shadow over plot
                                  ymax = Inf), alpha = 0.2) +
  facet_wrap(vars(regionGH),
             nrow=unique(df_long$regionGH) %>% length %>% sqrt %>% floor) +
  geom_line(aes(x=date, y=proportionUndocumentedQuadratic * coeff), size=0.8,
            color = propColor) + 
  xlab("") +
  scale_y_continuous(name = "Tests executed over total population\n",
                     sec.axis = sec_axis(
                       ~./coeff, name="Proportion of infectives undocumented\n")) +
  theme(
    axis.title = element_text(size=16),
    axis.text = element_text(size=14),
    axis.text.y = element_text(color = testedColor),
    axis.text.y.right = element_text(color = propColor),
    axis.title.y = element_text(color = testedColor, size=12),
    axis.title.y.right = element_text(color = propColor, size=12),
    panel.spacing = unit(0.8, "lines"))

ggsave("tamponiprop_vs_ft.pdf", plot = g, path=output_path, width = 10.8,
       height = 6.62, units = "in")

#### Plot comparing undocumented infectives (absolute) ####
for (region in unique(df_long$code)) {
  df_temp = df_long %>%
    filter(code == !!region)
  
  g = df_temp %>%
    ggplot() +
    geom_line(aes(x=date, y=activeInfectivesLinear, color="#0072B2"), size=0.8) + 
    geom_line(aes(x=date, y=activeInfectivesCubic, color="#D55E00"), size=0.8) +
    geom_line(aes(x=date, y=activeInfectivesQuadratic, color="#CC79A7"), size=0.8) +
    geom_line(aes(x=date, y=activeInfectivesQuadraticSix, color="#009E73"), size=0.8) +
    geom_line(aes(x=date, y=activeInfectivesQuadraticSixtyFive, color="#56B4E9"), size=0.8) +
    geom_line(aes(x=date, y=activeInfectives, color="#E69F00"), size=0.8) +
    geom_rect(data = tbl_rects, aes(xmin = xstart, xmax = xend, ymin = -Inf, # Shadow over plot
                                    ymax = Inf), alpha = 0.2) +
    xlab("") +
    ylab("Number of infectives") +
    scale_colour_manual(
      name = "Functional form",
      labels = c("Linear",
                 unname(TeX("Cubic ($\\gamma_1 = 0.6$, $\\gamma_2 = 0.8$)")),
                 unname(TeX("Quadratic $(\\gamma = 0.7)$")),
                 unname(TeX("Quadratic $(\\gamma = 0.6)$")),
                 unname(TeX("Quadratic $(\\gamma = 0.65)$")),
                 "Baseline"),
      values=c("#0072B2", # Dark blue
               "#D55E00", # Orange-brown
               "#CC79A7", # Pink
               "#009E73", # Green
               "#56B4E9", # Light blue
               "#E69F00")) + # Yellow
  theme(
    axis.title = element_text(size=16),
    axis.text = element_text(size=14),
    legend.text.align = 0)
  
  ggsave(glue("undocumented_comparison_{region}.pdf"), plot = g,
         path=output_path, width = 10.8, height = 6.62, units = "in")
}

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
g = df_long %>%
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
  ylab(TeX("Mean incidence rate ($\\pm 2 \\times$ standard error)")) +
  geom_point(aes(colour = factor(direction)), stat = "identity",
             position = position_dodge(width = 1), size=3) +
  geom_path(aes(colour = factor(direction)), alpha=0.3) +
  geom_errorbar(aes(ymin=mean-2*se, ymax=mean+2*se,
                    color=direction), width=.3) +
  facet_grid(.~ direction, scales = "free", switch = "x", space = "free_x") + 
  theme(
    axis.title = element_text(size=16),
    axis.text = element_text(size=14),
    axis.text.x = element_text(angle = 90, hjust = 1),
    strip.placement = "outside", legend.position = "none") +
  scale_colour_manual(values=c("#0072B2", # Dark blue
                               "#D55E00", # Orange-brown
                               "#CC79A7", # Pink
                               "#009E73", # Green
                               "#56B4E9", # Light blue
                               "#E69F00")) # Yellow

ggsave("heterogeneity_over_regions.pdf", plot = g,
       path=output_path, width = 10.8, height = 6.62, units = "in")

#### Infective rate per NUTS 1 region ####
g = df_long %>%
  group_by(direction, date) %>%
  mutate(infectiveRate = sum(activeInfectives) / sum(populationAlive),
         direction = as.factor(direction)) %>%
  ungroup() %>%
  select(date, direction, infectiveRate) %>%
  distinct() %>%
  ggplot() +
  # geom_rect(data = tbl_rects, aes(xmin = xstart, xmax = xend, ymin = -Inf, # Shadow below plot
  #                                 ymax = Inf), alpha = 0.2) +
  geom_bar(aes(x=date, y=infectiveRate, fill=direction), stat="identity",
           width=1) +
  facet_wrap(vars(direction)) +
  geom_rect(data = tbl_rects, aes(xmin = xstart, xmax = xend, ymin = -Inf, # Shadow over plot
                                  ymax = Inf), alpha = 0.2) +
  xlab("") +
  ylab("") +
  theme(
    axis.title = element_text(size=16),
    axis.text = element_text(size=14),
    legend.position="none") +
  scale_fill_manual(values=c("#0072B2", # Dark blue
                             "#D55E00", # Orange-brown
                             "#CC79A7", # Pink
                             "#009E73", # Green
                             "#56B4E9", # Light blue
                             "#E69F00")) # Yellow

ggsave("infective_rate_per_NUTS1.pdf", plot = g, path=output_path, width = 10.8,
       height = 6.62, units = "in")

#### Infective rate per NUTS 2 region for Nord-Est ####
for (direc in unique(df_long$direction)) {
  g = df_long %>%
    filter(direction == !!direc) %>%
    mutate(infectiveRate = activeInfectives / populationAlive,
           regionGH = as.factor(regionGH)) %>%
    select(date, regionGH, infectiveRate) %>%
    distinct() %>%
    ggplot() +
    # geom_rect(data = tbl_rects, aes(xmin = xstart, xmax = xend, ymin = -Inf, # Shadow below plot
    #                                 ymax = Inf), alpha = 0.2) +
    geom_bar(aes(x=date, y=infectiveRate, fill=regionGH), stat="identity",
             width=1) +
    facet_wrap(vars(regionGH)) +
    geom_rect(data = tbl_rects, aes(xmin = xstart, xmax = xend, ymin = -Inf, # Shadow over plot
                                    ymax = Inf), alpha = 0.2) +
    xlab("") +
    ylab("") +
    theme(
      axis.title = element_text(size=16),
      axis.text = element_text(size=14),
      legend.position="none") +
    scale_fill_manual(values=c("#0072B2", # Dark blue
                               "#D55E00", # Orange-brown
                               "#CC79A7", # Pink
                               "#009E73", # Green
                               "#56B4E9", # Light blue
                               "#E69F00")) # Yellow
  
  ggsave(glue("infective_rates_{direc}.pdf"), plot = g, path=output_path,
         width = 10.8, height = 6.62, units = "in")
}

#### Plot per day ####
for (direc in unique(df_long$direction)) {
  g = df_long %>%
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
    theme(
      axis.title = element_text(size=16),
      axis.text = element_text(size=14),
      axis.text.x = element_text(angle = 90, hjust = 1),
      strip.placement = "outside", legend.position = "none") +
    scale_colour_manual(values=c("#0072B2")) # Dark Blue
  
  ggsave(glue("infective_rates_weekday_{direc}.pdf"), plot = g,
         path=output_path, width = 10.8, height = 6.62, units = "in")
}
