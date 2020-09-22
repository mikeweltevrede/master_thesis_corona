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
library(snakecase)
library(tidyverse)

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
# Maximum incubation period
tau = 14

# Do we want to use a rolling window_size, i.e. only use the most recent
# `window_size` observations?
rolling = TRUE
window_size = 100 + tau

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

if (form %in% c("Linear", "Quadratic", "DownwardsVertex", "UpwardsVertex",
                "Cubic")){
  print(glue("#### Running models while modelling undocumented infections with",
             " the {form} functional form! ####"))
  infective_variable = glue("infectives{form}")
  undoc_flag = glue("_Undoc{form}")
  
} else if (form == ""){
  # Then do not use the undocumented infections modelling
  print("#### Running models WITHOUT modelling undocumented infections! ####")
  infective_variable = "infectives"
  undoc_flag = ""
  
} else {
  sprintf(paste("The variable `form` is %s but it should be one of %s.",
                "Choosing 'infectives' as the infectives variable, so *not*",
                "modelling undocumented infections."),
          form,
          paste(c("Linear", "Quadratic", "DownwardsVertex",
                  "UpwardsVertex", "Cubic", ""), collapse=", ")) %>%
    print
  
  infective_variable = "infectives"
  undoc_flag = ""
}

#### Within ####
fm = glue("{infective_variable} ~ -1 + ",
          "dplyr::lag({infective_variable}, {tau}):",
          "dplyr::lag(susceptibleRate, {tau})+",
          paste(M_regressors, collapse="+")) %>%
  as.formula

firstColor = "#0072B2" # Dark blue
secondColor = "#D55E00" # Orange-brown

start = 20
window_size = start + tau

fit_my_model = function(time_moment, fm, data, rolling, window_size, tau,
                        aic=TRUE) {
  
  unique_infectives = unique(data[(time_moment-window_size+1+tau):time_moment, ][[infective_variable]])
  
  if (length(unique_infectives) == 1) {
    # Then there is no variation and the model cannot estimate a parameter
    return(NA)
  }
  
  if (rolling) {
    if (aic) {
      model = step(lm(fm, data=data[(time_moment-window_size+1):time_moment, ]),
                   k=2, trace=0, scope=list(
                     "lower" = glue("{infective_variable} ~ -1 + ",
                                    "dplyr::lag({infective_variable}, {tau}):",
                                    "dplyr::lag(susceptibleRate, {tau})") %>%
                       as.formula,
                     "upper" = fm))
    } else {
      model = lm(fm, data=data[(time_moment-window_size+1):time_moment, ])
    }
  } else {
    if (aic) {
      model = step(lm(fm, data=head(data, time_moment)), k=2, trace=0,
                   scope=list(
                     "lower" = glue("{infective_variable} ~ -1 + ",
                                    "dplyr::lag({infective_variable}, {tau}):",
                                    "dplyr::lag(susceptibleRate, {tau})") %>%
                       as.formula,
                     "upper" = fm))
    } else {
      model = lm(fm, data=head(data, time_moment))
    }
  }
  
  return(predict(model, newdata=data[1:(time_moment+1), ])[[time_moment+1]])
}

for (direc in unique(df_meta$direction)){
  df_temp = df_long %>%
    filter(direction == !!direc)
  
  regions = unique(df_temp$code)
  plots = vector("list")
  
  for (region in regions){
    # Select only the data for the relevant region
    data = df_temp %>%
      filter(code == !!region)
    
    region_full = filter(df_meta, code == !!region)$regionGH[1]
    
    end = nrow(data)-1
    
    pred_infectives = sapply(window_size:end, fit_my_model, fm=fm, data=data,
                             rolling=TRUE, window_size=window_size, tau=tau,
                             aic=TRUE)
    true_infectives = data[[infective_variable]][window_size:end]
    
    tbl_temp = tibble(
      date = data$date[window_size:end],
      true = true_infectives,
      pred = pred_infectives
    )
    
    inds = which(is.na(tbl_temp$pred))
    starts = c(inds[1], inds[which(diff(inds) != 1)+1])
    ends = c(inds[which(diff(inds) != 1)], tail(inds, 1))
    
    tbl_rects = tibble(
      xstart = tbl_temp$date[starts-1],
      xend = tbl_temp$date[ends+1]
    )
    
    g_temp = tbl_temp %>% 
      ggplot() +
      geom_line(aes(x=date, y=true, color=firstColor)) +
      geom_line(aes(x=date, y=pred, color=secondColor, group = 1))
    
    if (nrow(tbl_rects) > 0) {
      g_temp = g_temp +
        geom_rect(data = tbl_rects, aes(xmin = xstart, xmax = xend, ymin = -Inf,
                                        ymax = Inf), alpha = 0.2)
    }
    
    g_temp = g_temp +
      xlab("") +
      ylab("Number of infectives \n") +
      ggtitle(region_full) +
      scale_colour_manual(name = "Infectives", 
                          labels = c("True", "Predicted"),
                          values = c(firstColor, secondColor)) +
      theme(
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
  
  g = do.call("grid.arrange", c(plots, ncol=floor(sqrt(length(plots))))) %>% 
    plot_grid(mylegend, ncol = 1, rel_heights = c(1, .2))
  
  ggsave(glue("model_within_lag{tau}_forecast_start{start}_{direc}{undoc_flag}",
              "{rolling_flag}.pdf"), plot = g, path = output_path)
}

#### Between ####
#### Data preparation ####
df_sumInc = tibble(date = as.Date(NA),
                   code = character(),
                   sumInfectives = numeric())

regions = unique(df_long$code)

for (region in regions){
  df_sumInc = df_sumInc %>%
    bind_rows(
      df_wide %>%
        select(
          map(regions[regions != region], starts_with, vars = colnames(.)) %>%
            unlist()) %>%
        select(ends_with(infective_variable)) %>%
        mutate_all(dplyr::lag, n=tau) %>%
        transmute(
          date = df_wide$date,
          code = region,
          sumInfectives = rowSums(.)))
}

df_long = df_long %>%
  left_join(df_sumInc, by = c("code", "date"))

#### Make plots ####
fm = glue("{infective_variable} ~ -1 +",
          "dplyr::lag({infective_variable}, {tau}):",
          "dplyr::lag(susceptibleRate, {tau})+",
          "dplyr::lag(susceptibleRate, {tau}):sumInfectives +",
          paste(M_regressors, collapse="+")) %>%
  as.formula

firstColor = "#0072B2" # Dark blue
secondColor = "#D55E00" # Orange-brown

start = 20
window_size = start + tau

fit_my_model = function(time_moment, fm, data, rolling, window_size, tau,
                        aic=TRUE) {
  
  unique_infectives = unique(data[(time_moment-window_size+1+tau):time_moment, ][[infective_variable]])
  
  if (length(unique_infectives) == 1) {
    # Then there is no variation and the model cannot estimate a parameter
    return(NA)
  }
  
  if (rolling) {
    if (aic) {
      model = step(lm(fm, data=data[(time_moment-window_size+1):time_moment, ]),
                   k=2, trace=0,
                     scope=list(
                       "lower" = glue(
                         "{infective_variable} ~ -1 +",
                         "dplyr::lag({infective_variable}, {tau}):",
                         "dplyr::lag(susceptibleRate, {tau})+",
                         "dplyr::lag(susceptibleRate, {tau}):sumInfectives") %>%
                         as.formula,
                       "upper" = fm))
    } else {
      model = lm(fm, data=data[(time_moment-window_size+1):time_moment, ])
    }
  } else {
    if (aic) {
      model = step(lm(fm, data=head(data, time_moment)), k=2, trace=0,
                   scope=list(
                     "lower" = glue(
                       "{infective_variable} ~ -1 +",
                       "dplyr::lag({infective_variable}, {tau}):",
                       "dplyr::lag(susceptibleRate, {tau})+",
                       "dplyr::lag(susceptibleRate, {tau}):sumInfectives") %>%
                       as.formula,
                     "upper" = fm))
    } else {
      model = lm(fm, data=head(data, time_moment))
    }
  }
  
  return(predict(model, newdata=data[1:(time_moment+1), ])[[time_moment+1]])
}

for (direc in unique(df_meta$direction)){
  df_temp = df_long %>%
    filter(direction == !!direc)
  
  regions = unique(df_temp$code)
  plots = vector("list")
  
  for (region in regions){
    # Select only the data for the relevant region
    data = df_temp %>%
      filter(code == !!region)
    
    region_full = filter(df_meta, code == !!region)$regionGH[1]
    
    end = nrow(data)-1
    
    pred_infectives = sapply(window_size:end, fit_my_model, fm=fm, data=data,
                             rolling=TRUE, window_size=window_size, tau=tau,
                             aic=TRUE)
    true_infectives = data[[infective_variable]][window_size:end]
    
    tbl_temp = tibble(
      date = data$date[window_size:end],
      true = true_infectives,
      pred = pred_infectives
    )
    
    inds = which(is.na(tbl_temp$pred))
    starts = c(inds[1], inds[which(diff(inds) != 1)+1])
    ends = c(inds[which(diff(inds) != 1)], tail(inds, 1))
    
    tbl_rects = tibble(
      xstart = tbl_temp$date[starts-1],
      xend = tbl_temp$date[ends+1]
    )
    
    g_temp = tbl_temp %>% 
      ggplot() +
      geom_line(aes(x=date, y=true, color=firstColor)) +
      geom_line(aes(x=date, y=pred, color=secondColor, group = 1))
    
    if (nrow(tbl_rects) > 0) {
      g_temp = g_temp +
        geom_rect(data = tbl_rects, aes(xmin = xstart, xmax = xend, ymin = -Inf,
                                        ymax = Inf), alpha = 0.2)
    }
    
    g_temp = g_temp +
      xlab("") +
      ylab("Number of infectives \n") +
      ggtitle(region_full) +
      scale_colour_manual(name = "Infectives", 
                          labels = c("True", "Predicted"),
                          values = c(firstColor, secondColor)) +
      theme(
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
  
  g = do.call("grid.arrange", c(plots, ncol=floor(sqrt(length(plots))))) %>% 
    plot_grid(mylegend, ncol = 1, rel_heights = c(1, .2))
  
  ggsave(glue("model_between_lag{tau}_forecast_start{start}_{direc}{undoc_flag}",
              "{rolling_flag}.pdf"), plot = g, path = output_path)
}
