library(ggplot2)
library(glue)
library(latex2exp)
library(tidyverse)

source("config.R")

df_meta = readxl::read_xlsx(path_metadata, sheet = "Metadata")
df_long = readr::read_csv(path_full_long, col_types = do.call(
  cols, list(date = col_date(format = "%Y-%m-%d")))) %>%
  left_join(df_meta %>% select(code, regionGH, direction), by="code") %>%
  mutate(weekday =
           lubridate::wday(date, label = TRUE) %>%
           as.factor)

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
