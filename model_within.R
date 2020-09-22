rm(list=ls())

#### Within-region spread model ####
# Delta i_rt = beta_within*S_rt-tau*Delta i_rt-tau + delta*X_rt + nu_rt
# i_rt is the absolute number of cases!

#### Set-up ####
# Import standard variables and activate Python environment
source("config.R")

# Import packages
library(glue)
library(gtools)
library(latex2exp)
library(plm)
library(snakecase)
library(tidyverse)
library(xtable)

# Import data
df_long = readr::read_csv(path_full_long, col_types = do.call(
  cols, list(date = col_date(format = "%Y-%m-%d"))))
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

# Function to add stars to the estimates according to the p-values, as well as
# putting the t-statistics between parentheses
output_for_table = function(model, method="ols", significance=6){
  
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
  
  if (method == "ols") {
    tvar = "t"
    tvardash = paste0(tvar, " ")
    
  } else {
    
    if (method == "random") {
      tvar = "z"
    } else {
      tvar = "t"
    }
    
    tvardash = paste0(tvar, "-")
  }
  
  stars = coef(summary(model))[, glue("Pr(>|{tvar}|)")] %>%
    sapply(get_stars)
  estimates = coefficients(model) %>%
    signif(significance) %>%
    paste0(stars)
  names(estimates) = names(coefficients(model))
  
  tvals = coef(summary(model))[, glue("{tvardash}value")] %>%
    signif(significance) %>%
    sapply(function(x){paste0("(", x, ")")})
  names(tvals) = names(coefficients(model))
  
  return(list("estimates"=estimates, "tvals"=tvals))
}

#### Models without model selection ####
# Create results table
results_table = tibble(variables = c(M_regressors, "beta"))

#### National model ####
# Construct formula
fm = glue("infectivesNational{form} ~ -1 + ",
          "dplyr::lag(infectivesNational{form}, {tau}):",
          "dplyr::lag(susceptibleRateNational, {tau}) +",
  paste(M_regressors, collapse="+")) %>%
  as.formula

if (rolling) {
  model = lm(fm, data=tail(df_wide, window_size))
} else {
  model = lm(fm, data=df_wide)
}

summary(model)

#### Make parameter table ####
# Retrieve parameter estimates
table_output = output_for_table(model)
estimates = table_output$estimates
tvals = table_output$tvals

# Insert parameter estimates in the results table
results_table = results_table %>%
  left_join(tibble("variables" = c(M_regressors, "beta"),
                   "National" = unname(
                     c(estimates[M_regressors],
                       estimates[
                         glue("dplyr::lag(infectivesNational{form}, {tau}):",
                              "dplyr::lag(susceptibleRateNational, {tau})")])),
                   "National_tvals" = unname(
                     c(tvals[M_regressors],
                       tvals[
                         glue("dplyr::lag(infectivesNational{form}, {tau}):",
                              "dplyr::lag(susceptibleRateNational, {tau})")]))),
            by="variables")

#### Pooled OLS ####
# Construct formula
fm = glue("{infective_variable} ~ -1 + ",
          "dplyr::lag({infective_variable}, {tau}):",
          "dplyr::lag(susceptibleRate, {tau})+",
          paste(M_regressors, collapse="+")) %>%
  as.formula

if (rolling) {
  model = plm(fm, data=df_long_trimmed, model="pooling",
              index = c("code", "date"))
} else {
  model = plm(fm, data=df_long, model="pooling",
              index = c("code", "date"))
}

summary(model)

table_output = output_for_table(model, method="pooling")
estimates = table_output$estimates
tvals = table_output$tvals

results_table = results_table %>%
  left_join(tibble("variables" = c(M_regressors, "beta"),
                   "National_POLS" := unname(
                     c(estimates[M_regressors],
                       estimates[
                         glue("dplyr::lag({infective_variable}, {tau}):",
                              "dplyr::lag(susceptibleRate, {tau})")])),
                   "National_tvals_POLS" = unname(
                     c(tvals[M_regressors],
                       tvals[
                         glue("dplyr::lag({infective_variable}, {tau}):",
                              "dplyr::lag(susceptibleRate, {tau})")]))),
            by="variables")

#### Regional models ####
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
    left_join(tibble("variables" = c(M_regressors, "beta"),
                     !!glue("{region}") := unname(
                       c(estimates[M_regressors],
                         estimates[
                           glue("dplyr::lag({infective_variable}, {tau}):",
                                "dplyr::lag(susceptibleRate, {tau})")])),
                     !!glue("{region}_tvals") := unname(
                       c(tvals[M_regressors],
                         tvals[
                           glue("dplyr::lag({infective_variable}, {tau}):",
                                "dplyr::lag(susceptibleRate, {tau})")]))),
              by="variables")
}

# Transpose and reorder columns
results_table = results_table %>%
  gather(region, val, 2:ncol(results_table)) %>%
  spread(names(results_table)[1], val) %>%
  select(region, beta, everything())

# Put the national results at the top
results_table = rbind(
  results_table %>% filter(str_detect(region, "National")),
  results_table %>% filter(!str_detect(region, "National"))) %>%
  column_to_rownames("region")

# Return a LaTeX table
table_no_ms = xtable(results_table, math.style.exponents = TRUE)
print(table_no_ms,
      file=glue("{output_path}/model_within_table_no_ms{undoc_flag}",
                "{rolling_flag}.txt"))

#### Models with model selection (AIC) ####
results_table_aic = tibble(variables = c(M_regressors, "beta"))

#### National model ####
# Construct formula
fm = glue("infectivesNational{form} ~ -1 + ",
          "dplyr::lag(infectivesNational{form}, {tau}):",
          "dplyr::lag(susceptibleRateNational, {tau})+",
          paste(M_regressors, collapse="+")) %>%
  as.formula

# Use AIC for model selection - scope says we want to always keep beta_within in
if (rolling) {
  model = step(lm(fm, data=tail(df_wide, window_size)), k=2, trace=0,
               scope=list(
                 "lower" = glue("infectivesNational{form} ~ -1 + ",
                                "dplyr::lag(infectivesNational{form}, {tau}):",
                                "dplyr::lag(susceptibleRateNational, {tau})") %>%
                   as.formula,
                 "upper" = fm))
} else {
  model = step(lm(fm, data=df_wide), k=2, trace=0,
               scope=list(
                 "lower" = glue("infectivesNational{form} ~ -1 + ",
                                "dplyr::lag(infectivesNational{form}, {tau}):",
                                "dplyr::lag(susceptibleRateNational, {tau})") %>%
                   as.formula,
                 "upper" = fm))
}

# Retrieve parameter estimates
table_output = output_for_table(model)
estimates = table_output$estimates
tvals = table_output$tvals

# Insert parameter estimates in the results table
results_table_aic = results_table_aic %>%
  left_join(tibble("variables" = c(M_regressors, "beta"),
                   "National" = unname(
                     c(estimates[M_regressors],
                       estimates[
                         glue("dplyr::lag(infectivesNational{form}, {tau}):",
                              "dplyr::lag(susceptibleRateNational, {tau})")])),
                   "National_tvals" = unname(
                     c(tvals[M_regressors],
                       tvals[
                         glue("dplyr::lag(infectivesNational{form}, {tau}):",
                              "dplyr::lag(susceptibleRateNational, {tau})")]))),
            by="variables")

#### Pooled OLS ####
aicbic_plm <- function(object, criterion) {
  # Credit: Rookie @ StackOverflow: https://stackoverflow.com/users/2525157/rookie
  # https://stackoverflow.com/questions/46186527/how-to-calculate-bic-and-aic-for-a-gmm-model-in-r-using-plm
  
  model_summary = summary(object)
  
  u.hat = residuals(model_summary) # extract residuals
  
  np = length(model_summary$coefficients[, 1]) # number of parameters
  n = nrow(model_summary$model) # number of data
  SSR  = log( (sum(u.hat^2)/(n))) # log sum of squares
  
  if (criterion == "AIC") {
    return(round(2*np + n*(log(2*pi) + SSR  + 1), 1))
  } else if (criterion == "BIC") {
    return(round(log(n)*np + n*(log(2*pi) + SSR + 1), 1))
  }
}

# Construct formula
fm_base = glue("{infective_variable} ~ -1 + ",
               "dplyr::lag({infective_variable}, {tau}):",
               "dplyr::lag(susceptibleRate, {tau})") %>%
  as.formula

if (rolling) {
  model = plm(fm_base, data=df_long_trimmed, model="pooling",
              index = c("code", "date"))
} else {
  model = plm(fm_base, data=df_long, model="pooling",
              index = c("code", "date"))
}

aic_best = aicbic_plm(model, "AIC")

for (i in 1:length(M_regressors)) {
  combos = combinations(length(M_regressors), i, v=M_regressors)
  
  for (r in 1:nrow(combos)) {
    fm_temp = glue("{infective_variable} ~ -1 + ",
                      "dplyr::lag({infective_variable}, {tau}):",
                      "dplyr::lag(susceptibleRate, {tau})+",
                      paste(combos[r, ], collapse="+")) %>%
      as.formula
    
    if (rolling) {
      model_temp = plm(fm_temp, data=df_long_trimmed, model="pooling",
                  index = c("code", "date"))
    } else {
      model_temp = plm(fm_temp, data=df_long, model="pooling",
                  index = c("code", "date"))
    }
    
    aic = aicbic_plm(model_temp, "AIC")
    
    if (aic < aic_best) {
      aic_best = aic
      model = model_temp
    }
  }
}

table_output = output_for_table(model, method="pooling")
estimates = table_output$estimates
tvals = table_output$tvals

results_table_aic = results_table_aic %>%
  left_join(tibble("variables" = c(M_regressors, "beta"),
                   "National_POLS" := unname(
                     c(estimates[M_regressors],
                       estimates[
                         glue("dplyr::lag({infective_variable}, {tau}):",
                              "dplyr::lag(susceptibleRate, {tau})")])),
                   "National_tvals_POLS" = unname(
                     c(tvals[M_regressors],
                       tvals[
                         glue("dplyr::lag({infective_variable}, {tau}):",
                              "dplyr::lag(susceptibleRate, {tau})")]))),
            by="variables")

#### Regional models ####
# Construct formula
fm = glue("{infective_variable} ~ -1 + ",
          "dplyr::lag({infective_variable}, {tau}):",
          "dplyr::lag(susceptibleRate, {tau})+",
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
                   "lower" = glue("{infective_variable} ~ -1 + ",
                                  "dplyr::lag({infective_variable}, {tau}):",
                                  "dplyr::lag(susceptibleRate, {tau})") %>%
                     as.formula,
                   "upper" = fm))
  } else {
    model = step(lm(fm, data=data), k=2, trace=0,
                 scope=list(
                   "lower" = glue("{infective_variable} ~ -1 + ",
                                  "dplyr::lag({infective_variable}, {tau}):",
                                  "dplyr::lag(susceptibleRate, {tau})") %>%
                     as.formula,
                   "upper" = fm))
  }
  
  # Retrieve parameter estimates
  table_output = output_for_table(model)
  estimates = table_output$estimates
  tvals = table_output$tvals
  
  # Insert parameter estimates in the results table
  results_table_aic = results_table_aic %>%
    left_join(tibble("variables" = c(M_regressors, "beta"),
                     !!glue("{region}") := unname(
                       c(estimates[M_regressors],
                         estimates[
                           glue("dplyr::lag({infective_variable}, {tau}):",
                                "dplyr::lag(susceptibleRate, {tau})")])),
                     !!glue("{region}_tvals") := unname(
                       c(tvals[M_regressors],
                         tvals[glue("dplyr::lag({infective_variable}, {tau}):",
                                    "dplyr::lag(susceptibleRate, {tau})")]))),
              by="variables")
}

# Transpose and reorder columns
results_table_aic = results_table_aic %>%
  gather(region, val, 2:ncol(results_table_aic)) %>%
  spread(names(results_table_aic)[1], val) %>%
  select(region, beta, everything())

# Put the national results at the top
results_table_aic = rbind(
  results_table_aic %>% filter(str_detect(region, "National")),
  results_table_aic %>% filter(!str_detect(region, "National"))) %>%
  column_to_rownames("region")

# Return a LaTeX table
table_ms_aic = xtable(results_table_aic, math.style.exponents = TRUE)
print(table_ms_aic,
      file=glue("{output_path}/model_within_table_aic{undoc_flag}",
                "{rolling_flag}.txt"))

#### Models with model selection (BIC) ####
results_table_bic = tibble(variables = c(M_regressors, "beta"))

#### National model ####
# Construct formula
fm = glue("infectivesNational{form} ~ -1 + ",
          "dplyr::lag(infectivesNational{form}, {tau}):",
          "dplyr::lag(susceptibleRateNational, {tau})+",
          paste(M_regressors, collapse="+")) %>%
  as.formula

# Use BIC for model selection - scope says we want to always keep beta_within in
if (rolling) {
  model = step(lm(fm, data=tail(df_wide, window_size)),
               k=log(window_size), trace=0, scope=list(
                 "lower" = glue("infectivesNational{form} ~ -1 + ",
                                "dplyr::lag(infectivesNational{form}, {tau}):",
                                "dplyr::lag(susceptibleRateNational, {tau})") %>%
                   as.formula,
                 "upper" = fm))
} else {
  model = step(lm(fm, data=df_wide), k=log(nrow(df_wide)), trace=0,
               scope=list(
                 "lower" = glue("infectivesNational{form} ~ -1 + ",
                                "dplyr::lag(infectivesNational{form}, {tau}):",
                                "dplyr::lag(susceptibleRateNational, {tau})") %>%
                   as.formula,
                 "upper" = fm))
}

# Retrieve parameter estimates
table_output = output_for_table(model)
estimates = table_output$estimates
tvals = table_output$tvals

# Insert parameter estimates in the results table
results_table_bic = results_table_bic %>%
  left_join(tibble("variables" = c(M_regressors, "beta"),
                   "National" = unname(
                     c(estimates[M_regressors],
                       estimates[
                         glue("dplyr::lag(infectivesNational{form}, {tau}):",
                              "dplyr::lag(susceptibleRateNational, {tau})")])),
                   "National_tvals" = unname(
                     c(tvals[M_regressors],
                       tvals[
                         glue("dplyr::lag(infectivesNational{form}, {tau}):",
                              "dplyr::lag(susceptibleRateNational, {tau})")]))),
            by="variables")

#### Pooled OLS ####
# Construct formula
fm_base = glue("{infective_variable} ~ -1 + ",
               "dplyr::lag({infective_variable}, {tau}):",
               "dplyr::lag(susceptibleRate, {tau})") %>%
  as.formula

if (rolling) {
  model = plm(fm_base, data=df_long_trimmed, model="pooling",
              index = c("code", "date"))
} else {
  model = plm(fm_base, data=df_long, model="pooling",
              index = c("code", "date"))
}

bic_best = aicbic_plm(model, "BIC")

for (i in 1:length(M_regressors)) {
  combos = combinations(length(M_regressors), i, v=M_regressors)
  
  for (r in 1:nrow(combos)) {
    model_temp = glue("{infective_variable} ~ -1 + ",
                      "dplyr::lag({infective_variable}, {tau}):",
                      "dplyr::lag(susceptibleRate, {tau})+",
                      paste(combos[r, ], collapse="+")) %>%
      as.formula
    
    if (rolling) {
      model_temp = plm(fm_temp, data=df_long_trimmed, model="pooling",
                       index = c("code", "date"))
    } else {
      model_temp = plm(fm_temp, data=df_long, model="pooling",
                       index = c("code", "date"))
    }
    
    bic = aicbic_plm(model_temp, "BIC")
    
    if (bic < bic_best) {
      bic_best = bic
      model = model_temp
    }
  }
}

table_output = output_for_table(model, method="pooling")
estimates = table_output$estimates
tvals = table_output$tvals

results_table_bic = results_table_bic %>%
  left_join(tibble("variables" = c(M_regressors, "beta"),
                   "National_POLS" := unname(
                     c(estimates[M_regressors],
                       estimates[
                         glue("dplyr::lag({infective_variable}, {tau}):",
                              "dplyr::lag(susceptibleRate, {tau})")])),
                   "National_tvals_POLS" = unname(
                     c(tvals[M_regressors],
                       tvals[
                         glue("dplyr::lag({infective_variable}, {tau}):",
                              "dplyr::lag(susceptibleRate, {tau})")]))),
            by="variables")

#### Regional models ####
# Construct formula
fm = glue("{infective_variable} ~ -1 + ",
          "dplyr::lag({infective_variable}, {tau}):",
          "dplyr::lag(susceptibleRate, {tau})+",
          paste(M_regressors, collapse="+")) %>%
  as.formula

for (region in regions){
  # Select only the data for the relevant region
  data = df_long %>%
    filter(code == !!region)
  
  # Use BIC for model selection
  if (rolling) {
    model = step(lm(fm, data=tail(data, window_size)), k=log(window_size),
                 trace=0, scope=list(
                   "lower" = glue("{infective_variable} ~ -1 + ",
                                  "dplyr::lag({infective_variable}, {tau}):",
                                  "dplyr::lag(susceptibleRate, {tau})") %>%
                     as.formula,
                   "upper" = fm))
  } else {
    model = step(lm(fm, data=data), k=log(nrow(data)), trace=0,
                 scope=list(
                   "lower" = glue("{infective_variable} ~ -1 + ",
                                  "dplyr::lag({infective_variable}, {tau}):",
                                  "dplyr::lag(susceptibleRate, {tau})") %>%
                     as.formula,
                   "upper" = fm))
  }
  
  # Retrieve parameter estimates
  table_output = output_for_table(model)
  estimates = table_output$estimates
  tvals = table_output$tvals
  
  # Insert parameter estimates in the results table
  results_table_bic = results_table_bic %>%
    left_join(tibble("variables" = c(M_regressors, "beta"),
                     !!glue("{region}") := unname(
                       c(estimates[M_regressors],
                         estimates[glue("dplyr::lag({infective_variable}, {tau}):",
                                        "dplyr::lag(susceptibleRate, {tau})")])),
                     !!glue("{region}_tvals") := unname(
                       c(tvals[M_regressors],
                         tvals[glue("dplyr::lag({infective_variable}, {tau}):",
                                    "dplyr::lag(susceptibleRate, {tau})")]))),
              by="variables")
}

# Transpose and reorder columns
results_table_bic = results_table_bic %>%
  gather(region, val, 2:ncol(results_table_bic)) %>%
  spread(names(results_table_bic)[1], val) %>%
  select(region, beta, everything())

# Put the national results at the top
results_table_bic = rbind(
  results_table_bic %>% filter(str_detect(region, "National")),
  results_table_bic %>% filter(!str_detect(region, "National"))) %>%
  column_to_rownames("region")

# Return a LaTeX table
table_ms_bic = xtable(results_table_bic, math.style.exponents = TRUE)
print(table_ms_bic,
      file=glue("{output_path}/model_within_table_bic{undoc_flag}",
                "{rolling_flag}.txt"))

#### Plot beta over time ####
df_meta = readxl::read_xlsx(path_metadata, sheet = "Metadata")

fm = glue("{infective_variable} ~ -1 + ",
          "dplyr::lag({infective_variable}, {tau}):",
          "dplyr::lag(susceptibleRate, {tau})+",
          paste(M_regressors, collapse="+")) %>%
  as.formula

#### Without model selection ####
tbl_beta = tibble(date = as.Date(NA), betas=numeric(0), code=character(0))

# Find the estimates of beta per region over time
for (region in regions){
  betas = vector("double")
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
    beta = coef(model)[[glue("dplyr::lag({infective_variable}, {tau}):",
                              "dplyr::lag(susceptibleRate, {tau})")]]
    betas = c(betas, beta)
  }
  
  # Append the results to the table
  tbl_beta = tbl_beta %>%
    bind_rows(tibble(date = data$date[window_size:nrow(data)],
                     betas = betas,
                     code = region) %>%
                drop_na())
}

# Add region and direction to the table
tbl_beta = tbl_beta %>%
  left_join(df_meta %>% select(c(region, regionGH, code, direction)), by="code")

# Make a plot per direction
for (sub_tbl in split(tbl_beta, tbl_beta$direction)){
  direc = sub_tbl$direction[1]
  g = ggplot(sub_tbl, aes(date, betas, colour = regionGH)) + 
    geom_point() +
    geom_smooth(method="loess", span=0.3, se=FALSE) +
    xlab("") +
    ylab(TeX("$\\beta_{within}$")) +
    scale_colour_manual(name = "Region",
                        values=c("#0072B2", # Dark blue
                                 "#D55E00", # Orange-brown
                                 "#CC79A7", # Pink
                                 "#009E73", # Green
                                 "#56B4E9", # Light blue
                                 "#E69F00")) + # Yellow
    theme(
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12)
    )
  
  ggsave(glue("model_within_lag{tau}_betawithin_{direc}{undoc_flag}",
              "{rolling_flag}.pdf"), plot = g,
         path=output_path, width = 10.8, height = 6.62, units = "in")
}

#### With model selection (AIC) ####
tbl_beta = tibble(date = as.Date(NA), betas=numeric(0), code=character(0))

# Find the estimates of beta per region over time
for (region in regions){
  betas = vector("double")
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
                     "lower" = glue("{infective_variable} ~ -1 + ",
                                    "dplyr::lag({infective_variable}, {tau}):",
                                    "dplyr::lag(susceptibleRate, {tau})") %>%
                       as.formula,
                     "upper" = fm))
    } else {
      model = step(lm(fm, data=head(data, t)), k=2, trace=0,
                   scope=list(
                     "lower" = glue("{infective_variable} ~ -1 + ",
                                    "dplyr::lag({infective_variable}, {tau}):",
                                    "dplyr::lag(susceptibleRate, {tau})") %>%
                       as.formula,
                     "upper" = fm))
    }
    
    # Retrieve the beta estimate and append this to the list of betas
    beta = coef(model)[[glue("dplyr::lag({infective_variable}, {tau}):",
                              "dplyr::lag(susceptibleRate, {tau})")]]
    betas = c(betas, beta)
  }
  
  # Append the results to the table
  tbl_beta = tbl_beta %>%
    bind_rows(tibble(date = data$date[window_size:nrow(data)],
                     betas = betas,
                     code = region) %>%
                drop_na())
}

# Add region and direction to the table
tbl_beta = tbl_beta %>%
  left_join(df_meta %>% select(c(region, regionGH, code, direction)), by="code")

# Make a plot per direction
for (sub_tbl in split(tbl_beta, tbl_beta$direction)){
  direc = sub_tbl$direction[1]
  g = ggplot(sub_tbl, aes(date, betas, colour = regionGH)) + 
    geom_point() +
    geom_smooth(method="loess", span=0.3, se=FALSE) +
    xlab("") +
    ylab(TeX("$\\beta_{within}$")) +
    scale_colour_manual(name = "Region",
                        values=c("#0072B2", # Dark blue
                                 "#D55E00", # Orange-brown
                                 "#CC79A7", # Pink
                                 "#009E73", # Green
                                 "#56B4E9", # Light blue
                                 "#E69F00")) + # Yellow
    theme(
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12)
    )
  
  ggsave(glue("model_within_lag{tau}_betawithin_{direc}_aic{undoc_flag}",
              "{rolling_flag}.pdf"), plot = g,
         path=output_path, width = 10.8, height = 6.62, units = "in")
}

#### With model selection (BIC) ####
tbl_beta = tibble(date = as.Date(NA), betas=numeric(0), code=character(0))

# Find the estimates of beta per region over time
for (region in regions){
  betas = vector("double")
  dates = vector("character")
  
  # Select only the data for the relevant region
  data = df_long %>% filter(code == !!region)
  
  for (t in window_size:nrow(data)){
    # Use BIC for model selection - scope says we want to always keep
    # beta_within in
    if (rolling) {
      model = step(lm(fm, data=data[(t-window_size+1):t, ]), k=log(window_size),
                   trace=0, scope=list(
                     "lower" = glue("{infective_variable} ~ -1 + ",
                                    "dplyr::lag({infective_variable}, {tau}):",
                                    "dplyr::lag(susceptibleRate, {tau})") %>%
                       as.formula,
                     "upper" = fm))
    } else {
      model = step(lm(fm, data=head(data, t)), k=log(t), trace=0,
                   scope=list(
                     "lower" = glue("{infective_variable} ~ -1 + ",
                                    "dplyr::lag({infective_variable}, {tau}):",
                                    "dplyr::lag(susceptibleRate, {tau})") %>%
                       as.formula,
                     "upper" = fm))
      
    }
    # Retrieve the beta estimate and append this to the list of betas
    beta = coef(model)[[glue("dplyr::lag({infective_variable}, {tau}):",
                              "dplyr::lag(susceptibleRate, {tau})")]]
    betas = c(betas, beta)
  }
  
  # Append the results to the table
  tbl_beta = tbl_beta %>%
    bind_rows(tibble(date = data$date[window_size:nrow(data)],
                     betas = betas,
                     code = region) %>%
                drop_na())
}

# Add region and direction to the table
tbl_beta = tbl_beta %>%
  left_join(df_meta %>% select(c(region, regionGH, code, direction)), by="code")

# Make a plot per direction
for (sub_tbl in split(tbl_beta, tbl_beta$direction)){
  direc = sub_tbl$direction[1]
  g = ggplot(sub_tbl, aes(date, betas, colour = regionGH)) + 
    geom_point() +
    geom_smooth(method="loess", span=0.3, se=FALSE) +
    xlab("") +
    ylab(TeX("$\\beta_{within}$")) +
    scale_colour_manual(values=c("#0072B2", # Dark blue
                                 "#D55E00", # Orange-brown
                                 "#CC79A7", # Pink
                                 "#009E73", # Green
                                 "#56B4E9", # Light blue
                                 "#E69F00")) + # Yellow
    theme(
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12)
    )
  
  ggsave(glue("model_within_lag{tau}_betawithin_{direc}_bic{undoc_flag}",
              "{rolling_flag}.pdf"), plot = g,
         path=output_path, width = 10.8, height = 6.62, units = "in")
}
#### Forecasts ####
firstColor = "#0072B2" # Dark blue
secondColor = "#D55E00" # Orange-brown

start = 20
window_size = start + tau

fit_my_model = function(time_moment, fm, data, rolling, window_size, tau,
                        aic=FALSE) {
  
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
                                    "dplyr::lag(susceptibleRate{undoc_type}, {tau})") %>%
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
                                    "dplyr::lag(susceptibleRate{undoc_type}, {tau})") %>%
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
    
    pred_infectives = numeric(0)
    
    for (time_moment in window_size:end) {
      pred_infectives = c(pred_infectives,
                          fit_my_model(time_moment, fm=fm, data=data,
                                       rolling=rolling, window_size=window_size,
                                       tau=tau, aic=FALSE))
    }
    
    # pred_infectives = sapply(window_size:end, fit_my_model, fm=fm, data=data,
    #                          rolling=TRUE, window_size=window_size, tau=tau,
    #                          aic=TRUE)
    true_infectives = data[[infective_variable]][window_size:end]
    
    tbl_temp = tibble(
      date = data$date[window_size:end],
      true = true_infectives,
      pred = pred_infectives
    )
    
    inds = which(is.na(tbl_temp$pred))
    starts = c(inds[1], inds[which(diff(inds) != 1)+1])
    ends = c(inds[which(diff(inds) != 1)], tail(inds, 1))
    
    tbl_rects_NA = tibble(
      xstart = tbl_temp$date[starts-1],
      xend = tbl_temp$date[ends+1]
    )
    
    g_temp = tbl_temp %>% 
      ggplot() +
      geom_line(aes(x=date, y=true, color=firstColor)) +
      geom_line(aes(x=date, y=pred, color=secondColor, group = 1))
    
    if (nrow(tbl_rects_NA) > 0) {
      g_temp = g_temp +
        geom_rect(data = tbl_rects_NA, aes(xmin = xstart, xmax = xend,
                                           ymin = -Inf, ymax = Inf), alpha = 0.2)
    }
    
    if (lockdown_end >= min(tbl_temp$date)) {
      g_temp = g_temp +
        geom_rect(data = tbl_rects, aes(xmin = xstart, xmax = xend, ymin = -Inf, # Shadow over plot
                                        ymax = Inf), alpha = 0.2)
    }
    
    g_temp = g_temp +
      xlab("") +
      ylab("Infectives \n") +
      ggtitle(region_full) +
      scale_colour_manual(name = "Infectives", 
                          labels = c("True", "Predicted"),
                          values = c(firstColor, secondColor)) +
      theme(
        axis.title = element_text(size=16),
        axis.text = element_text(size=14),
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
      theme(
        legend.position = "none",
        axis.title = element_text(size=16),
        axis.text = element_text(size=14)
        )
    
    plots[[region]] = g_temp
  }
  
  g = do.call("grid.arrange", c(plots, ncol=floor(sqrt(length(plots))))) %>% 
    plot_grid(mylegend, ncol = 1, rel_heights = c(1, .2))
    
  ggsave(glue("model_within_lag{tau}_forecast_start{start}_{direc}{undoc_flag}",
              "{rolling_flag}{window_flag}.pdf"), plot = g, path = output_path, width = 10.8,
         height = 6.62, units = "in")
}
