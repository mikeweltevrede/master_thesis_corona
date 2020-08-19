library(ggplot2)
library(glue)
library(latex2exp)
library(tidyverse)

source("config.R")

df_meta = readxl::read_xlsx(path_metadata, sheet = "Metadata")
df_long = readr::read_csv(path_full_long, col_types = do.call(
  cols, list(date = col_date(format = "%Y-%m-%d")))) %>%
  left_join(df_meta %>% select(code, regionGH, direction), by="code")

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

#### New ####
for (direc in unique(df_long$direction)){
  num_regions = df_long %>%
    filter(direction == direc) %>%
    .[["regionGH"]] %>%
    unique %>%
    length
  g = ggplot(df_long %>% filter(direction == direc),
             aes(date, infectivesRate, colour = regionGH, group = regionGH)) +
    geom_line() +
    ggtitle(direc) +
    xlab("") +
    ylab(TeX("Infectives rate")) +
    scale_colour_manual(values=c("#0072B2", # Dark blue
                                 "#D55E00", # Orange-brown
                                 "#CC79A7", # Pink
                                 "#009E73", # Green
                                 "#56B4E9", # Light blue
                                 "#E69F00")) # Yellow
  print(g)
  ggsave(glue("infective_rate_{direc}.pdf"), path=output_path,
         width = 10.8, height = 7.47, units = "in")
}
