rm(list=ls())

#### Model 3 - Within and between-region spread ####
# I_rt = alpha_within*I_rt-tau*S_rt-tau +
#        alpha_within*S_rt-tau*\sum_{c!=r}I_ct-tau +
#        X_rt*delta + nu_rt

#### Setup ####
# Import standard variables
source("config.R")

# Import packages
library(glue)
library(latex2exp)
library(restriktor)
library(snakecase)
library(tidyverse)
library(xtable)

# Import data
df_long = readr::read_csv(path_full_long, col_types = do.call(
  cols, list(date = col_date(format = "%Y-%m-%d"))))
df_wide = readr::read_csv(path_full_wide, col_types = do.call(
  cols, list(date=col_date(format="%Y-%m-%d"))))

pivot_to_df_wide = function(df_long) {
  # Turn long data into wide data
  df_wide = df_long %>%
    pivot_wider(names_from = code,
                values_from = all_of(colnames(df_long)[-c(1,2)]))
  
  # Column names are now of the form "variable_regionCode". We want to have them
  # in the form "regionCode_variable".
  colnames(df_wide) = colnames(df_wide) %>%
    sapply(function(s) {
      s %>%
        str_split("_", simplify = TRUE) %>%
        rev %>%
        paste(collapse="_")
    }, USE.NAMES=FALSE)
  
  # Order the columns alphabetically, keeping date at the start
  df_wide = df_wide[, c("date", sort(colnames(df_wide[-1])))]
  
  return(df_wide)
}

# Incubation period
lag = 5

# Get all region abbreviations
regions = df_long$code %>% unique

# Select regressors
X_regressors = c("weekend")
all_variables = c("(Intercept)", paste0(X_regressors, 1))

# Should we apply the restriktor package to ensure that the alphas are positive?
restrict = TRUE

#### Data preprocessing ####
# Transform the variable to include undocumented infections, if applicable. Note
# that this does not depend on the region specifically but only the values at
# that moment of time. As such, we do not need to loop and can simply apply the
# function to each row.
form = "" %>%
  to_upper_camel_case

if (form %in% c("Linear", "Quadratic", "DownwardsVertex", "UpwardsVertex",
                "Cubic")){
  infective_variable = glue("infectives{form}") %>%
    snakecase::to_lower_camel_case()
  
  undoc_flag = glue("_Undoc{form}")
  
  print(glue("####Running models while modelling undocumented infections with ",
             "the {form} functional form!####"))
  
} else if (form == ""){
  # Then do not use the undocumented infections modelling
  infective_variable = "infectives"
  undoc_flag = ""
  
  print("####Running models WITHOUT modelling undocumented infections!####")
  
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
        mutate_all(dplyr::lag, n=lag) %>%
        transmute(
          date = df_wide$date,
          code = region,
          sumInfectives = rowSums(.)))
}

df_long = df_long %>%
  left_join(df_sumInc, by = c("code", "date"))

# Add weekend effect
df_long = df_long %>%
  mutate(populationDensity = totalPopulation/area)

df_wide = pivot_to_df_wide(df_long) %>%
  mutate(weekend =
           lubridate::wday(date, label = TRUE) %in% c("Sat", "Sun") %>%
           as.integer %>% as.factor,
         lockdown =
           ifelse(date > as.Date("2020-03-10", format = "%Y-%m-%d") &
                    date < as.Date("2020-06-03", format = "%Y-%m-%d"),
                  1, 0) %>%
           as.factor)

df_long = df_long %>%
  mutate(weekend =
           lubridate::wday(date, label = TRUE) %in% c("Sat", "Sun") %>%
           as.integer %>% as.factor,
         lockdown =
           ifelse(date > as.Date("2020-03-10", format = "%Y-%m-%d") &
                    date < as.Date("2020-06-03", format = "%Y-%m-%d"),
                  1, 0) %>%
           as.factor)

# Add nationwide variables by summing the individual regions' variables
susceptibleTotal = df_wide %>%
  select(ends_with("susceptiblePopulation")) %>%
  rowSums
total = df_wide %>%
  select(ends_with("totalPopulation")) %>%
  rowSums
infectivesTotal = df_wide %>% 
  select(ends_with(glue("_{infective_variable}"))) %>% 
  rowSums
area = df_wide %>%
  select(ends_with("area")) %>%
  rowSums
df_wide = df_wide %>%
  mutate(susceptibleRateTotal = susceptibleTotal/total,
         infectivesTotal = infectivesTotal,
         populationDensityTotal = total/area)

#### Models without model selection ####
# Note that a national model does not make sense in model 3
results_table = tibble(variables = c(all_variables, "alpha_w", "alpha_b"))

# Construct formula
fm = glue("{infective_variable} ~ ",
          "lag({infective_variable}, {lag}):lag(susceptibleRate, {lag})+",
          "lag(susceptibleRate, {lag}):sumInfectives +",
          paste(X_regressors, collapse="+")) %>%
  as.formula

for (region in regions){
  # Select only the data for the relevant region
  data = df_long %>% filter(code == !!region)
  
  # Estimate the model by OLS
  model = lm(fm, data=data)
  
  # pdf(glue("{output_path}/model3_lag{lag}_lmplot_{region}{undoc_flag}.pdf"))
  # par(mfrow=c(2,2))
  # plot(model)
  # par(mfrow=c(1,1))
  # dev.off()
  
  # TODO: Check the maths behind this
  if (restrict) {
    model = restriktor(model,
                       constraints = rbind(c(0, 1, 0, 0), # alpha_within > 0
                                           c(0, 0, 1, 0)), # alpha_between > 0
                       rhs = c(0,0))
  }
  
  # Retrieve parameter estimates
  estimates = coef(summary(model))[, "Estimate"]
  pvals = coef(summary(model))[, "Pr(>|t|)"] # TODO: SE with stars
  
  # Insert parameter estimates in the results table
  results_table = results_table %>%
    left_join(tibble("variables" = c(all_variables, "alpha_w", "alpha_b"),
                     !!glue("{region}") := unname(
                       c(estimates[all_variables],
                         estimates[glue("lag({infective_variable}, {lag}):",
                                        "lag(susceptibleRate, {lag})")],
                         estimates[glue("lag(susceptibleRate, {lag}):",
                                        "sumInfectives")])),
                     !!glue("{region}_pvals") := unname(
                       c(pvals[all_variables],
                         pvals[glue("lag({infective_variable}, {lag}):",
                                    "lag(susceptibleRate, {lag})")],
                         pvals[glue("lag(susceptibleRate, {lag}):",
                                    "sumInfectives")]))),
              by="variables")
}

# Transpose and reorder columns
results_table = results_table %>%
  gather(region, val, 2:ncol(results_table)) %>%
  spread(names(results_table)[1], val) %>%
  select(region, alpha_w, alpha_b, everything()) %>%
  column_to_rownames("region")

# Return a LaTeX table
table_no_ms = xtable(results_table, math.style.exponents = TRUE)
table_no_ms

#### Models with model selection (AIC) ####
results_table_aic = tibble(variables = c(all_variables, "alpha_w", "alpha_b"))

# Construct formula
fm = glue("{infective_variable} ~ ",
          "lag({infective_variable}, {lag}):lag(susceptibleRate, {lag})+",
          "lag(susceptibleRate, {lag}):sumInfectives +",
          paste(X_regressors, collapse="+")) %>%
  as.formula

for (region in regions){
  # Select only the data for the relevant region
  data = df_long %>% filter(code == !!region)
  
  # Use AIC for model selection
  model = step(lm(fm, data=data), k=2, trace=0,
               scope=list("lower" = glue(
                 "{infective_variable} ~ ",
                 "lag({infective_variable}, {lag}):lag(susceptibleRate, {lag})+",
                 "lag(susceptibleRate, {lag}):sumInfectives") %>%
                   as.formula,
                 "upper" = fm))
  
  # pdf(glue("{output_path}/model3_lag{lag}_lmplot_{region}_aic{undoc_flag}.pdf"))
  # par(mfrow=c(2,2))
  # plot(model)
  # par(mfrow=c(1,1))
  # dev.off()
  
  # TODO: Check the maths behind this
  if (restrict) {
    model = restriktor(model,
                       constraints = rbind(c(0, 1, 0, 0), # alpha_within > 0
                                           c(0, 0, 1, 0)), # alpha_between > 0
                       rhs = c(0,0))
  }
  
  # Retrieve parameter estimates
  estimates = coef(summary(model))[, "Estimate"]
  pvals = coef(summary(model))[, "Pr(>|t|)"] # TODO: SE with stars
  
  # Insert parameter estimates in the results table
  results_table_aic = results_table_aic %>%
    left_join(tibble("variables" = c(all_variables, "alpha_w", "alpha_b"),
                     !!glue("{region}") := unname(
                       c(estimates[all_variables],
                         estimates[glue("lag({infective_variable}, {lag}):",
                                        "lag(susceptibleRate, {lag})")],
                         estimates[glue("lag(susceptibleRate, {lag}):",
                                        "sumInfectives")])),
                     !!glue("{region}_pvals") := unname(
                       c(pvals[all_variables],
                         pvals[glue("lag({infective_variable}, {lag}):",
                                    "lag(susceptibleRate, {lag})")],
                         pvals[glue("lag(susceptibleRate, {lag}):",
                                    "sumInfectives")]))),
              by="variables")
}

# Transpose and reorder columns
results_table_aic = results_table_aic %>%
  gather(region, val, 2:ncol(results_table_aic)) %>%
  spread(names(results_table_aic)[1], val) %>%
  select(region, alpha_w, alpha_b, everything()) %>%
  column_to_rownames("region")

# Return a LaTeX table
table_aic = xtable(results_table_aic, math.style.exponents = TRUE)
table_aic

#### Models with model selection (BIC) ####
results_table_bic = tibble(variables = c(all_variables, "alpha_w", "alpha_b"))

# Construct formula
fm = glue("{infective_variable} ~ ",
          "lag({infective_variable}, {lag}):lag(susceptibleRate, {lag})+",
          "lag(susceptibleRate, {lag}):sumInfectives +",
          paste(X_regressors, collapse="+")) %>%
  as.formula

for (region in regions){
  # Select only the data for the relevant region
  data = df_long %>% filter(code == !!region)
  
  # Use BIC for model selection
  model = step(lm(fm, data=data), k=log(nrow(data)), trace=0,
               scope=list("lower" = glue(
                 "{infective_variable} ~ ",
                 "lag({infective_variable}, {lag}):lag(susceptibleRate, {lag})+",
                 "lag(susceptibleRate, {lag}):sumInfectives") %>%
                   as.formula,
                 "upper" = fm))
  
  # pdf(glue("{output_path}/model3_lag{lag}_lmplot_{region}_bic{undoc_flag}.pdf"))
  # par(mfrow=c(2,2))
  # plot(model)
  # par(mfrow=c(1,1))
  # dev.off()
  
  # TODO: Check the maths behind this
  if (restrict) {
    model = restriktor(model,
                       constraints = rbind(c(0, 1, 0, 0), # alpha_within > 0
                                           c(0, 0, 1, 0)), # alpha_between > 0
                       rhs = c(0,0))
  }
  
  # Retrieve parameter estimates
  estimates = coef(summary(model))[, "Estimate"]
  pvals = coef(summary(model))[, "Pr(>|t|)"] # TODO: SE with stars
  
  # Insert parameter estimates in the results table
  results_table_bic = results_table_bic %>%
    left_join(tibble("variables" = c(all_variables, "alpha_w", "alpha_b"),
                     !!glue("{region}") := unname(
                       c(estimates[all_variables],
                         estimates[glue("lag({infective_variable}, {lag}):",
                                        "lag(susceptibleRate, {lag})")],
                         estimates[glue("lag(susceptibleRate, {lag}):",
                                        "sumInfectives")])),
                     !!glue("{region}_pvals") := unname(
                       c(pvals[all_variables],
                         pvals[glue("lag({infective_variable}, {lag}):",
                                    "lag(susceptibleRate, {lag})")],
                         pvals[glue("lag(susceptibleRate, {lag}):",
                                    "sumInfectives")]))),
              by="variables")
}

# Transpose and reorder columns
results_table_bic = results_table_bic %>%
  gather(region, val, 2:ncol(results_table_bic)) %>%
  spread(names(results_table_bic)[1], val) %>%
  select(region, alpha_w, alpha_b, everything()) %>%
  column_to_rownames("region")

# Return a LaTeX table
table_bic = xtable(results_table_bic, math.style.exponents = TRUE)
table_bic
