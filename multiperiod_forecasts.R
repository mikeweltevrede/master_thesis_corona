rm(list=ls())

#### Within-region spread model ####
# Delta i_rt = beta_within*S_rt-tau*Delta i_rt-tau + delta*X_rt + nu_rt
# i_rt is the absolute number of cases!

#### Set-up ####
# Import standard variables and activate Python environment
source("config.R")

# Import packages
library(cowplot)
library(glue)
library(gridExtra)
library(gtools)
library(latex2exp)
library(plm)
library(snakecase)
library(tidyverse)
library(xtable)

# Import data
df_meta = readxl::read_xlsx(path_metadata, sheet = "Metadata")
df_long = readr::read_csv(path_full_long, col_types = do.call(
  cols, list(date = col_date(format = "%Y-%m-%d")))) %>%
  left_join(df_meta %>% select(code, regionGH, direction), by="code")
df_wide = readr::read_csv(path_full_wide, col_types = do.call(
  cols, list(date = col_date(format = "%Y-%m-%d"))))

# Get all region abbreviations
regions = unique(df_long$code)

# Select regressors
M_regressors = c("weekend")

#### Decide the parameters ####
# You need to adapt, if desired, the following parameters:
# tau: int, the maximum incubation period
# rolling: boolean, whether to apply a rolling window
# window_size: int, if using a rolling window, how large?
# form: str, the form of undocumented infections to model with (if any)

# Maximum incubation period
tau = 14

# Do we want to use a rolling window_size, i.e. only use the most recent
# `window_size` observations?
rolling = TRUE
window_size = 30

if (window_size == 100){
  window_flag = ""
} else {
  window_flag = glue("window{window_size}")
}

window_size = window_size + tau

if (rolling) {
  rolling_flag = "_rolling"
  
  df_long_trimmed = df_long %>% 
    group_by(code) %>% 
    slice(tail(row_number(), window_size)) %>%
    ungroup()
  
} else {
  rolling_flag = ""
}

# Determine if we want to model undocumented infectives and, if so, by which
# method. Note that infective_variable is the number of new cases, i.e. \Delta i
form = "" %>%
  to_upper_camel_case
name_gamma = ""

if (!(name_gamma %in% c("Six", "SixtyFive", "", "SeventyFive"))) {
  sprintf(paste("The variable `name_gamma` is %s but it should be one",
                "of %s. Setting gamma to 0.7, the default."),
          name_gamma,
          paste(c("Six", "SixtyFive", "", "SeventyFive"),
                collapse=", "))
  
  name_gamma = ""
}

if (form %in% c("Linear", "Quadratic", "DownwardsVertex", "UpwardsVertex",
                "Cubic")){
  print(glue("#### Running models while modelling undocumented infections with",
             " the {form}{name_gamma} functional form! ####"))
  
  infective_variable = glue("activeInfectives{form}{name_gamma}")
  undoc_flag = glue("_Undoc{form}{name_gamma}")
  undoc_type = glue("{form}{name_gamma}")
  
} else if (form == ""){
  # Then do not use the undocumented infections modelling
  print("#### Running models WITHOUT modelling undocumented infections! ####")
  infective_variable = "activeInfectives"
  undoc_flag = ""
  undoc_type = ""
  
} else {
  sprintf(paste("The variable `form` is %s but it should be one of %s.",
                "Choosing 'infectives' as the infectives variable, so *not*",
                "modelling undocumented infections."),
          form,
          paste(c("Linear", "Quadratic", "DownwardsVertex",
                  "UpwardsVertex", "Cubic", ""), collapse=", ")) %>%
    print
  
  infective_variable = "activeInfectives"
  undoc_flag = ""
  undoc_type = ""
}

df_long = df_long %>%
  mutate(
    susceptibleLag = dplyr::lag(!!sym(glue("susceptibleRate{form}")), tau),
    infectivesLag = dplyr::lag(!!sym(glue("{infective_variable}")), tau)
    )

#### Forecast ####
fm = glue("{infective_variable} ~ -1 + ",
          "dplyr::lag({infective_variable}, {tau}):",
          "dplyr::lag(susceptibleRate{undoc_type}, {tau}) +",
          paste(M_regressors, collapse="+")) %>%
  as.formula

#### NEW ####
# Table to shade the lockdown
lockdown_start = head(df_wide$date[which(df_wide$lockdown == 1)], 1)
lockdown_end = tail(df_wide$date[which(df_wide$lockdown == 1)], 1)

tbl_rects = tibble(
  xstart = lockdown_start,
  xend = lockdown_end
)

for (direc in unique(df_meta$direction)){
  df_temp = df_long %>%
    filter(direction == !!direc)
  
  regions = unique(df_temp$code)
  plots = vector("list")
  
  for (region in regions){
    region_full = filter(df_meta, code == !!region)$regionGH[1]
    
    df_region = df_temp %>%
      filter(code == !!region) %>%
      transmute(
        date,
        weekend,
        I = !!sym(glue("{infective_variable}")),
        s = !!sym(glue("susceptibleRate{form}")),
        S = !!sym(glue("susceptiblePopulation{form}")),
        N = populationAlive,
        lagI = dplyr::lag(I, !!tau),
        lag_s = dplyr::lag(s, !!tau),
        lagS = dplyr::lag(S, 1),
        susceptibleRate_pred = s
      ) %>%
      drop_na() %>%
      mutate(I_pred = NA,
             I_pred = as.numeric(I_pred),
             group = NA,
             group = as.character(group))
    
    # Given data until April 23 (30+tau data points), what can we forecast?
    # - The tau days after the end: till May 7
    # - Use the forecasts to predict beyond May 7
    model = lm(I ~ -1 + lag_s:lagI + weekend, data=df_region[1:window_size, ])
    beta = coef(model)[["lag_s:lagI"]]
    delta = coef(model)[["weekend"]]
    
    for (row in (window_size+1):(window_size+1+tau)) {
      data_row = df_region[row, ]
      I_pred = round(beta*data_row$lagI*data_row$lag_s + delta*data_row$weekend)
      S_new = df_region[row-1, "S"] - I_pred
      
      df_region[row, "I_pred"] = I_pred
      df_region[row, "S"] = S_new
      df_region[row, "s"] = S_new / data_row$N
      df_region[row, "group"] = "Forecast"
    }
    
    for (row in (window_size+2+tau):(window_size+2+2*tau)) {
      data_row = df_region[row, ]
      I_pred = round(beta*data_row$lagI*data_row$lag_s + delta*data_row$weekend)
      S_new = df_region[row-1, "S"] - I_pred
      
      df_region[row, "I_pred"] = I_pred
      df_region[row, "S"] = S_new
      df_region[row, "s"] = S_new / data_row$N
      df_region[row, "group"] = "Two-step"
    }
    
    g_temp = df_region %>%
      filter(date <= date[(window_size+2+2*tau)]) %>%
      ggplot() +
      geom_line(aes(x=date, y=I, colour="True", group=1), size=1) +
      geom_line(aes(x=date, y=I_pred, colour=group, group=1), size=1) +
      xlab("") +
      ylab("Infectives\n") +
      ggtitle(region_full)
    
    if (lockdown_end >= min(df_region$date)) {
      g_temp = g_temp +
        geom_rect(data = tbl_rects, aes(xmin = xstart, xmax = xend, ymin = -Inf, # Shadow over plot
                                        ymax = Inf), alpha = 0.2)
    }
    
    g_temp = g_temp +  
      scale_colour_manual(name = "",
                          breaks=c("True", "Forecast", "Two-step"),
                          values=c("Forecast" = "#0072B2", # Dark blue
                                   "Two-step" = "#D55E00", # Orange-brown
                                   "True" = "#CC79A7")) + # Pink
      theme(
        axis.title = element_text(size=14),
        axis.text = element_text(size=12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)
      )
    
    if (region == regions[1]) {
      # Extract legend
      # https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
      g_legend = function(a.gplot){
        tmp = ggplot_gtable(ggplot_build(a.gplot))
        leg = which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
        legend = tmp$grobs[[leg]]
        return(legend)}
      
      mylegend = g_legend(g_temp)
    }
    
    g_temp = g_temp +
      theme(legend.position = "none")
    
    plots[[region]] = g_temp
  }
  
  g = do.call("grid.arrange", c(plots, ncol=ceiling(sqrt(length(plots))))) %>%
    plot_grid(mylegend, ncol = 1, rel_heights = c(1, .2))
  
  ggsave(glue("model_within_lag{tau}_forecast_full_{direc}{undoc_flag}.pdf"),
         plot = g, path = output_path, width = 10.8, height = 6.62, units = "in")
}
