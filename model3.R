rm(list=ls())

#### Model 3 - Within and between-region spread ####
# Delta Y_rt = beta_within*Delta Y_rt-tau*S_rt-tau +
#              beta_between*S_rt-tau*\sum_{c!=r}Delta Y_ct-tau +
#              delta*M_rt + nu_rt
# Y_rt is the absolute number of new cases!

#### Set-up ####
# Import standard variables and activate Python environment
source("config.R")

# Import packages
library(glue)
library(latex2exp)
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
# tau: int, maximum incubation period
# rolling: boolean, whether to apply a rolling window
# window_size: int, if using a rolling window, how large?
# form: str, the form of undocumented infections to model with (if any)

# Maximum incubation period
tau = 14

# Do we want to use a rolling window_size, i.e. only use the most recent
# `window_size` observations?
rolling = TRUE
window_size = 100 + tau

if (rolling) {
  rolling_flag = "_rolling"
} else {
  rolling_flag = ""
}

# Determine if we want to model undocumented infectives and, if so, by which
# method. Note that infective_variable is the number of new cases, i.e. Delta X.
form = "Quadratic" %>%
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

#### Data preparation ####
df_sumInc = tibble(date = as.Date(NA),
                   code = character(),
                   sumInfectives = numeric())

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

#### Models without model selection ####
# Retrieve parameter estimates
output_for_table = function(model, significance=6){
  
  get_stars = function(pval) {
    if (pval < 0.01) {
      stars = "***"
    } else if (pval < 0.05) {
      stars = "**"
    } else if (pval < 0.1) {
      stars = "*"
    } else {
      stars = ""
    }
  }
  
  stars = coef(summary(model))[, "Pr(>|t|)"] %>%
    sapply(get_stars)
  estimates = coefficients(model) %>%
    signif(significance) %>%
    paste0(stars)
  names(estimates) = names(coefficients(model))
  
  tvals = coef(summary(model))[, "t value"] %>%
    signif(significance) %>%
    sapply(function(x){paste0("(", x, ")")})
  names(tvals) = names(coefficients(model))
  
  return(list("estimates"=estimates, "tvals"=tvals))
}

# Create results table
results_table = tibble(variables = c(M_regressors, "beta_w", "beta_b"))

# Construct formula
fm = glue("{infective_variable} ~ -1 +",
          "lag({infective_variable}, {tau}):lag(susceptibleRate, {tau})+",
          "lag(susceptibleRate, {tau}):sumInfectives +",
          paste(M_regressors, collapse="+")) %>%
  as.formula

for (region in regions){
  # Select only the data for the relevant region
  data = df_long %>%
    filter(code == !!region)
  
  # Estimate the model by OLS
  if (rolling) {
    model = lm(fm, data=tail(data, window_size))
  } else {
    model = lm(fm, data=data)
  }
  
  # Retrieve parameter estimates
  table_output = output_for_table(model)
  estimates = table_output$estimates
  tvals = table_output$tvals
  
  # Insert parameter estimates in the results table
  results_table = results_table %>%
    left_join(tibble("variables" = c(M_regressors, "beta_w", "beta_b"),
                     !!glue("{region}") := unname(
                       c(estimates[M_regressors],
                         estimates[glue("lag({infective_variable}, {tau}):",
                                        "lag(susceptibleRate, {tau})")],
                         estimates[glue("lag(susceptibleRate, {tau}):",
                                        "sumInfectives")])),
                     !!glue("{region}_tvals") := unname(
                       c(tvals[M_regressors],
                         tvals[glue("lag({infective_variable}, {tau}):",
                                    "lag(susceptibleRate, {tau})")],
                         tvals[glue("lag(susceptibleRate, {tau}):",
                                    "sumInfectives")]))),
              by="variables")
}

# Transpose and reorder columns
results_table = results_table %>%
  gather(region, val, 2:ncol(results_table)) %>%
  spread(names(results_table)[1], val) %>%
  select(region, beta_w, beta_b, everything()) %>%
  column_to_rownames("region")

# Return a LaTeX table
table_no_ms = xtable(results_table, math.style.exponents = TRUE)
print(table_no_ms,
      file=glue("{output_path}/model3_table_no_ms{undoc_flag}{rolling_flag}.txt"))

#### Models with model selection (AIC) ####
results_table_aic = tibble(variables = c(M_regressors, "beta_w", "beta_b"))

# Construct formula
fm = glue("{infective_variable} ~ -1 +",
          "lag({infective_variable}, {tau}):lag(susceptibleRate, {tau})+",
          "lag(susceptibleRate, {tau}):sumInfectives +",
          paste(M_regressors, collapse="+")) %>%
  as.formula

for (region in regions){
  # Select only the data for the relevant region
  data = df_long %>%
    filter(code == !!region)
  
  # Use AIC for model selection
  if (rolling) {
    model = step(lm(fm, data=tail(data, window_size)), k=2, trace=0,
                 scope=list(
                   "lower" = glue(
                     "{infective_variable} ~ -1 +",
                     "lag({infective_variable}, {tau}):",
                     "lag(susceptibleRate, {tau})+",
                     "lag(susceptibleRate, {tau}):sumInfectives") %>%
                     as.formula,
                   "upper" = fm))
  } else {
    model = step(lm(fm, data=data), k=2, trace=0,
                 scope=list(
                   "lower" = glue(
                     "{infective_variable} ~ -1 +",
                     "lag({infective_variable}, {tau}):",
                     "lag(susceptibleRate, {tau})+",
                     "lag(susceptibleRate, {tau}):sumInfectives") %>%
                     as.formula,
                   "upper" = fm))
  }
  
  # Retrieve parameter estimates
  table_output = output_for_table(model)
  estimates = table_output$estimates
  tvals = table_output$tvals
  
  # Insert parameter estimates in the results table
  results_table_aic = results_table_aic %>%
    left_join(tibble("variables" = c(M_regressors, "beta_w", "beta_b"),
                     !!glue("{region}") := unname(
                       c(estimates[M_regressors],
                         estimates[glue("lag({infective_variable}, {tau}):",
                                        "lag(susceptibleRate, {tau})")],
                         estimates[glue("lag(susceptibleRate, {tau}):",
                                        "sumInfectives")])),
                     !!glue("{region}_tvals") := unname(
                       c(tvals[M_regressors],
                         tvals[glue("lag({infective_variable}, {tau}):",
                                    "lag(susceptibleRate, {tau})")],
                         tvals[glue("lag(susceptibleRate, {tau}):",
                                    "sumInfectives")]))),
              by="variables")
}

# Transpose and reorder columns
results_table_aic = results_table_aic %>%
  gather(region, val, 2:ncol(results_table_aic)) %>%
  spread(names(results_table_aic)[1], val) %>%
  select(region, beta_w, beta_b, everything()) %>%
  column_to_rownames("region")

# Return a LaTeX table
table_aic = xtable(results_table_aic, math.style.exponents = TRUE)
print(table_aic,
      file=glue("{output_path}/model3_table_aic{undoc_flag}{rolling_flag}.txt"))

#### Models with model selection (BIC) ####
results_table_bic = tibble(variables = c(M_regressors, "beta_w", "beta_b"))

# Construct formula
fm = glue("{infective_variable} ~ -1 +",
          "lag({infective_variable}, {tau}):lag(susceptibleRate, {tau})+",
          "lag(susceptibleRate, {tau}):sumInfectives +",
          paste(M_regressors, collapse="+")) %>%
  as.formula

for (region in regions){
  # Select only the data for the relevant region
  data = df_long %>%
    filter(code == !!region)
  
  # Use BIC for model selection
  if (rolling) {
    model = step(lm(fm, data=tail(data, window_size)), k=log(window_size), trace=0,
                 scope=list(
                   "lower" = glue(
                     "{infective_variable} ~ -1 +",
                     "lag({infective_variable}, {tau}):",
                     "lag(susceptibleRate, {tau})+",
                     "lag(susceptibleRate, {tau}):sumInfectives") %>%
                     as.formula,
                   "upper" = fm))
  } else {
    model = step(lm(fm, data=data), k=log(nrow(data)), trace=0,
                 scope=list(
                   "lower" = glue(
                     "{infective_variable} ~ -1 +",
                     "lag({infective_variable}, {tau}):",
                     "lag(susceptibleRate, {tau})+",
                     "lag(susceptibleRate, {tau}):sumInfectives") %>%
                     as.formula,
                   "upper" = fm))
  }
  
  # Retrieve parameter estimates
  table_output = output_for_table(model)
  estimates = table_output$estimates
  tvals = table_output$tvals
  
  # Insert parameter estimates in the results table
  results_table_bic = results_table_bic %>%
    left_join(tibble("variables" = c(M_regressors, "beta_w", "beta_b"),
                     !!glue("{region}") := unname(
                       c(estimates[M_regressors],
                         estimates[glue("lag({infective_variable}, {tau}):",
                                        "lag(susceptibleRate, {tau})")],
                         estimates[glue("lag(susceptibleRate, {tau}):",
                                        "sumInfectives")])),
                     !!glue("{region}_tvals") := unname(
                       c(tvals[M_regressors],
                         tvals[glue("lag({infective_variable}, {tau}):",
                                    "lag(susceptibleRate, {tau})")],
                         tvals[glue("lag(susceptibleRate, {tau}):",
                                    "sumInfectives")]))),
              by="variables")
}

# Transpose and reorder columns
results_table_bic = results_table_bic %>%
  gather(region, val, 2:ncol(results_table_bic)) %>%
  spread(names(results_table_bic)[1], val) %>%
  select(region, beta_w, beta_b, everything()) %>%
  column_to_rownames("region")

# Return a LaTeX table
table_bic = xtable(results_table_bic, math.style.exponents = TRUE)
print(table_bic,
      file=glue("{output_path}/model3_table_bic{undoc_flag}{rolling_flag}.txt"))

#### Plot beta over time ####
fm = glue("{infective_variable} ~ -1 +",
          "lag({infective_variable}, {tau}):lag(susceptibleRate, {tau})+",
          "lag(susceptibleRate, {tau}):sumInfectives +",
          paste(M_regressors, collapse="+")) %>%
  as.formula

withinColor = "#0072B2" # Dark blue
betweenColor = "#D55E00" # Orange-brown

#### Without model selection ####
tbl_beta = tibble(date = as.Date(NA), Within=numeric(0), Between=numeric(0),
                  code=character(0))

# Find the estimates of beta per region over time
for (region in regions){
  betas_w = vector("double")
  betas_b = vector("double")
  dates = vector("character")
  
  # Select only the data for the relevant region
  data = df_long %>%
    filter(code == !!region)
  
  for (t in window_size:nrow(data)){
    # Estimate the model by OLS
    if (rolling) {
      model = lm(fm, data=data[(t-window_size+1):t, ])
    } else {
      model = lm(fm, data=head(data, t))
    }
    
    # Retrieve the beta estimate and append this to the list of betas
    beta_w = coef(model)[[glue("lag({infective_variable}, {tau}):",
                               "lag(susceptibleRate, {tau})")]]
    beta_b = coef(model)[[glue("lag(susceptibleRate, {tau}):sumInfectives")]]
    betas_w = c(betas_w, beta_w)
    betas_b = c(betas_b, beta_b)
  }
  
  # Append the results to the table
  tbl_beta = tbl_beta %>%
    bind_rows(tibble(date = data$date[window_size:nrow(data)],
                     Within = betas_w,
                     Between = betas_b,
                     code = region) %>%
                drop_na())
}

# Add region and direction to the table
tbl_beta = tbl_beta %>%
  left_join(df_meta %>% select(c(regionGH, code, direction)), by="code")

# Make a plot per direction; facet_wrap by regions
for (sub_tbl in split(tbl_beta, tbl_beta$direction)){
  direc = sub_tbl$direction[1]
  coeff = mean(sub_tbl$Within / sub_tbl$Between)
  
  g = sub_tbl %>%
    ggplot(aes(x = date)) + 
    geom_point(aes(y = Within, color = withinColor)) +
    geom_smooth(aes(y = Within, color = withinColor), method="loess", span=0.3,
                se=FALSE)  +
    facet_wrap("regionGH") +
    geom_point(aes(y = Between*coeff, color = betweenColor)) +
    geom_smooth(aes(y = Between*coeff, color = betweenColor), method="loess",
                span=0.3, se=FALSE)  +
    xlab("") +
    scale_y_continuous(name = TeX("$\\beta_{within}$\n"),
                       sec.axis = sec_axis(
                         ~./coeff, TeX("$\\beta_{between}$\n"))) +
    scale_colour_manual(name = "Beta", 
                        labels = c("Within", "Between"),
                        values=c(withinColor, betweenColor)) + 
    theme(
      axis.text.y = element_text(color = withinColor),
      axis.text.y.right = element_text(color = betweenColor),
      axis.title.y = element_text(color = withinColor, size=12),
      axis.title.y.right = element_text(color = betweenColor, size=12),
      panel.spacing = unit(0.8, "lines"))
  print(g)
  
  ggsave(
    glue("model3_lag{tau}_betas_{direc}{undoc_flag}{rolling_flag}.pdf"),
    path=output_path, width = 10.8, height = 6.62, units = "in")
}

#### With model selection (AIC) ####
tbl_beta = tibble(date = as.Date(NA), Within=numeric(0), Between=numeric(0),
                  code=character(0))

# Find the estimates of beta per region over time
for (region in regions){
  betas_w = vector("double")
  betas_b = vector("double")
  dates = vector("character")
  
  # Select only the data for the relevant region
  data = df_long %>%
    filter(code == !!region)
  
  for (t in window_size:nrow(data)){
    # Use AIC for model selection - scope says we want to always keep
    # beta_within in
    if (rolling) {
      model = step(lm(fm, data=data[(t-window_size+1):t, ]), k=2, trace=0,
                   scope=list(
                     "lower" = glue(
                       "{infective_variable} ~ -1 +",
                       "lag({infective_variable}, {tau}):",
                       "lag(susceptibleRate, {tau})+",
                       "lag(susceptibleRate, {tau}):sumInfectives") %>%
                       as.formula,
                     "upper" = fm))
    } else {
      model = step(lm(fm, data=head(data, t)), k=2, trace=0,
                   scope=list(
                     "lower" = glue(
                       "{infective_variable} ~ -1 +",
                       "lag({infective_variable}, {tau}):",
                       "lag(susceptibleRate, {tau})+",
                       "lag(susceptibleRate, {tau}):sumInfectives") %>%
                       as.formula,
                     "upper" = fm))
    }
    
    # Retrieve the beta estimate and append this to the list of betas
    beta_w = coef(model)[[glue("lag({infective_variable}, {tau}):",
                               "lag(susceptibleRate, {tau})")]]
    beta_b = coef(model)[[glue("lag(susceptibleRate, {tau}):sumInfectives")]]
    betas_w = c(betas_w, beta_w)
    betas_b = c(betas_b, beta_b)
  }
  
  # Append the results to the table
  tbl_beta = tbl_beta %>%
    bind_rows(tibble(date = data$date[window_size:nrow(data)],
                     Within = betas_w,
                     Between = betas_b,
                     code = region) %>%
                drop_na())
}

# Add region and direction to the table
tbl_beta = tbl_beta %>%
  left_join(df_meta %>% select(c(regionGH, code, direction)), by="code")

# Make a plot per direction; facet_wrap by regions
for (sub_tbl in split(tbl_beta, tbl_beta$direction)){
  direc = sub_tbl$direction[1]
  coeff = mean(sub_tbl$Within / sub_tbl$Between)
  
  g = sub_tbl %>%
    ggplot(aes(x = date)) + 
    geom_point(aes(y = Within, color = withinColor)) +
    geom_smooth(aes(y = Within, color = withinColor), method="loess", span=0.3,
                se=FALSE)  +
    facet_wrap("regionGH") +
    geom_point(aes(y = Between*coeff, color = betweenColor)) +
    geom_smooth(aes(y = Between*coeff, color = betweenColor), method="loess",
                span=0.3, se=FALSE)  +
    xlab("") +
    scale_y_continuous(name = TeX("$\\beta_{within}$\n"),
                       sec.axis = sec_axis(
                         ~./coeff, TeX("$\\beta_{between}$\n"))) +
    scale_colour_manual(name = "Beta", 
                        labels = c("Within", "Between"),
                        values=c(withinColor, betweenColor)) + 
    theme(
      axis.text.y = element_text(color = withinColor),
      axis.text.y.right = element_text(color = betweenColor),
      axis.title.y = element_text(color = withinColor, size=12),
      axis.title.y.right = element_text(color = betweenColor, size=12),
      panel.spacing = unit(0.8, "lines"))
  print(g)
  
  ggsave(
    glue("model3_lag{tau}_betas_{direc}_aic{undoc_flag}{rolling_flag}.pdf"),
    path=output_path, width = 10.8, height = 6.62, units = "in")
}

#### With model selection (BIC) ####
tbl_beta = tibble(date = as.Date(NA), Within=numeric(0), Between=numeric(0),
                  code=character(0))

# Find the estimates of beta per region over time
for (region in regions){
  betas_w = vector("double")
  betas_b = vector("double")
  dates = vector("character")
  
  # Select only the data for the relevant region
  data = df_long %>%
    filter(code == !!region)
  
  for (t in window_size:nrow(data)){
    # Use BIC for model selection - scope says we want to always keep
    # beta_within in
    if (rolling) {
      model = step(lm(fm, data=data[(t-window_size+1):t, ]), k=log(window_size), trace=0,
                   scope=list(
                     "lower" = glue(
                       "{infective_variable} ~ -1 +",
                       "lag({infective_variable}, {tau}):",
                       "lag(susceptibleRate, {tau})+",
                       "lag(susceptibleRate, {tau}):sumInfectives") %>%
                       as.formula,
                     "upper" = fm))
    } else {
      model = step(lm(fm, data=head(data, t)), k=log(t), trace=0,
                   scope=list(
                     "lower" = glue(
                       "{infective_variable} ~ -1 +",
                       "lag({infective_variable}, {tau}):",
                       "lag(susceptibleRate, {tau})+",
                       "lag(susceptibleRate, {tau}):sumInfectives") %>%
                       as.formula,
                     "upper" = fm))
    }
    
    # Retrieve the beta estimate and append this to the list of betas
    beta_w = coef(model)[[glue("lag({infective_variable}, {tau}):",
                               "lag(susceptibleRate, {tau})")]]
    beta_b = coef(model)[[glue("lag(susceptibleRate, {tau}):sumInfectives")]]
    betas_w = c(betas_w, beta_w)
    betas_b = c(betas_b, beta_b)
  }
  
  # Append the results to the table
  tbl_beta = tbl_beta %>%
    bind_rows(tibble(date = data$date[window_size:nrow(data)],
                     Within = betas_w,
                     Between = betas_b,
                     code = region) %>%
                drop_na())
}

# Add region and direction to the table
tbl_beta = tbl_beta %>%
  left_join(df_meta %>% select(c(regionGH, code, direction)), by="code")

# Make a plot per direction; facet_wrap by regions
for (sub_tbl in split(tbl_beta, tbl_beta$direction)){
  direc = sub_tbl$direction[1]
  coeff = mean(sub_tbl$Within / sub_tbl$Between)
  
  g = sub_tbl %>%
    ggplot(aes(x = date)) + 
    geom_point(aes(y = Within, color = withinColor)) +
    geom_smooth(aes(y = Within, color = withinColor), method="loess", span=0.3,
                se=FALSE)  +
    facet_wrap("regionGH") +
    geom_point(aes(y = Between*coeff, color = betweenColor)) +
    geom_smooth(aes(y = Between*coeff, color = betweenColor), method="loess",
                span=0.3, se=FALSE)  +
    xlab("") +
    scale_y_continuous(name = TeX("$\\beta_{within}$\n"),
                       sec.axis = sec_axis(
                         ~./coeff, TeX("$\\beta_{between}$\n"))) +
    scale_colour_manual(name = "Beta", 
                        labels = c("Within", "Between"),
                        values=c(withinColor, betweenColor)) + 
    theme(
      axis.text.y = element_text(color = withinColor),
      axis.text.y.right = element_text(color = betweenColor),
      axis.title.y = element_text(color = withinColor, size=12),
      axis.title.y.right = element_text(color = betweenColor, size=12),
      panel.spacing = unit(0.8, "lines"))
  print(g)
  
  ggsave(
    glue("model3_lag{tau}_betas_{direc}_bic{undoc_flag}{rolling_flag}.pdf"),
    path=output_path, width = 10.8, height = 6.62, units = "in")
}
