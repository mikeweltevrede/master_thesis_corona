#### Model 1 - Within-region spread ####
# We start with a simple model ignoring effects across regions:
# I_rt = alpha_within*I_rt-tau*S_rt-tau + X_rt*delta + nu_rt

#### Setup ####
# source("clean_full.R") # May error; if so, run it by hand

# Import standard variables
source("config.R")

# Import packages
library(tidyverse)
library(glue)
library(xtable)
library(latex2exp)

# Import data
df_long = readr::read_csv(path_full_long, col_types = do.call(
  cols, list(date = col_date(format = "%Y-%m-%d"))))
df_wide_full = readr::read_csv(path_full_wide, col_types = do.call(
  cols, list(date=col_date(format="%Y-%m-%d"))))

# Incubation period
lag = 5

X_regressors = c("weekend", "weekNumber")
all_variables = c("(Intercept)", X_regressors) %>%
  str_replace("weekend", "weekend1")

#### Data preprocessing ####
# Add nationwide variables by summing the individual regions' variables
susceptible = df_wide_full %>%
  select(ends_with("susceptiblePopulation")) %>%
  rowSums
total = df_wide_full %>%
  select(ends_with("totalPopulation")) %>%
  rowSums

# We use the *_testedPositive columns because these have been cleaned
confirmed = df_wide_full %>% 
  select(ends_with("_testedPositive")) %>% 
  rowSums
df_wide_full = df_wide_full %>%
  mutate(susceptibleRateTotal = susceptible/total,
         confirmedTotal = confirmed)

# Add weekend and weekday effect
df_wide_full = df_wide_full %>%
  mutate(weekNumber = lubridate::week(df_wide_full$date)) %>%
  mutate(weekend = lubridate::wday(df_wide_full$date, label = TRUE)
         %in% c("Sat", "Sun") %>% as.integer %>% as.factor)

df_long = df_long %>%
  mutate(weekNumber = lubridate::week(df_long$date)) %>%
  mutate(weekend = lubridate::wday(df_long$date, label = TRUE)
         %in% c("Sat", "Sun") %>% as.integer %>% as.factor)

#### Least Squares Dummy Variables (LSDV) regression ####
fm = paste("confirmedTotal ~ ",
           glue("lag(confirmedTotal, {lag}):lag(susceptibleRateTotal, ",
                "{lag})+"),
           paste(X_regressors,
                 collapse="+")) %>%
  paste("+factor(code)") %>%
  as.formula

model = lm(fm, data=df_long)
summary(model)

png(glue("{output_path}/model1_lag{lag}_lmplot_lsdv.png"))
par(mfrow=c(2,2))
plot(model)
par(mfrow=c(1,1))
dev.off()

#### Models without model selection ####
results_table = tibble(variables = c(all_variables, "alpha"))

#### National model ####
# Construct formula
fm = paste("confirmedTotal ~ ",
           glue("lag(confirmedTotal, {lag}):lag(susceptibleRateTotal, ",
                "{lag})+"),
           paste(X_regressors,
                 collapse="+")) %>%
  as.formula

# Estimate the model by OLS
model = lm(fm, data=df_wide_full)

png(glue("{output_path}/model1_lag{lag}_lmplot_national.png"))
par(mfrow=c(2,2))
plot(model)
par(mfrow=c(1,1))
dev.off()

# Retrieve parameter estimates
estimates = coef(summary(model))[, "Estimate"]
pvals = coef(summary(model))[, "Pr(>|t|)"] # TODO: SE with stars

# Insert parameter estimates in the results table
results_table = results_table %>%
  left_join(tibble("variables" = c(all_variables, "alpha"),
                   "National" = unname(
                     c(estimates[all_variables], estimates[
                       glue("lag(confirmedTotal, {lag}):",
                            "lag(susceptibleRateTotal, {lag})")])),
                   "National_pvals" = unname(
                     c(pvals[all_variables], pvals[
                       glue("lag(confirmedTotal, {lag}):",
                            "lag(susceptibleRateTotal, {lag})")]))),
            by="variables")

#### Regional models ####
regions = df_long$code %>% unique

# Construct formula
fm = paste("testedPositive ~ ",
           glue("lag(testedPositive, {lag}):lag(susceptibleRate, {lag})+"),
           paste(X_regressors, collapse="+")) %>%
  as.formula

for (region in regions){
  # Select only the data for the relevant region
  data = df_long %>% filter(code == !!region)
  
  # Estimate the model by OLS
  model = lm(fm, data=data)
  
  png(glue("{output_path}/model1_lag{lag}_lmplot_{region}.png"))
  par(mfrow=c(2,2))
  plot(model)
  par(mfrow=c(1,1))
  dev.off()
  
  # Retrieve parameter estimates
  estimates = coef(summary(model))[, "Estimate"]
  pvals = coef(summary(model))[, "Pr(>|t|)"] # TODO: SE with stars
  
  # Insert parameter estimates in the results table
  results_table = results_table %>%
    left_join(tibble("variables" = c(all_variables, "alpha"),
                     !!glue("{region}") := unname(
                       c(estimates[all_variables], estimates[
                         glue("lag(testedPositive, {lag}):",
                              "lag(susceptibleRate, {lag})")])),
                     !!glue("{region}_pvals") := unname(
                       c(pvals[all_variables], pvals[
                         glue("lag(testedPositive, {lag}):",
                              "lag(susceptibleRate, {lag})")]))),
              by="variables")
}

# Transpose and reorder columns
results_table = results_table %>%
  gather(region, val, 2:ncol(results_table)) %>%
  spread(names(results_table)[1], val) %>%
  select(region, alpha, everything())

# Put the national results at the top
results_table = rbind(
  results_table %>% filter(str_detect(region, "National")),
  results_table %>% filter(!str_detect(region, "National"))) %>%
  column_to_rownames("region")

# Return a LaTeX table
table_no_ms = xtable(results_table, math.style.exponents = TRUE)
table_no_ms

#### Models with model selection (BIC) ####
results_table_ms = tibble(variables = c(all_variables, "alpha"))

#### National model ####
# Construct formula
fm = paste("confirmedTotal ~ ",
           glue("lag(confirmedTotal, {lag}):lag(susceptibleRateTotal, ",
                "{lag})+"),
           paste(X_regressors,
                 collapse="+")) %>%
  as.formula

# Use BIC for model selection - scope says we want to always keep
# alpha_within in
model = step(lm(fm, data=df_wide_full), k=log(nrow(df_wide_full)), trace=0,
             scope=list("lower" = paste("testedPositive ~ ",
                                        glue("lag(testedPositive, {lag}):",
                                             "lag(susceptibleRate, {lag})")) %>%
                          as.formula,
                        "upper" = fm))

png(glue("{output_path}/model1_lag{lag}_lmplot_national_bic.png"))
par(mfrow=c(2,2))
plot(model)
par(mfrow=c(1,1))
dev.off()

# Retrieve parameter estimates
estimates = coef(summary(model))[, "Estimate"]
pvals = coef(summary(model))[, "Pr(>|t|)"] # TODO: SE with stars

# Insert parameter estimates in the results table
results_table_ms = results_table_ms %>%
  left_join(tibble("variables" = c(all_variables, "alpha"),
                   "National" = unname(
                     c(estimates[all_variables], estimates[
                       glue("lag(confirmedTotal, {lag}):",
                            "lag(susceptibleRateTotal, {lag})")])),
                   "National_pvals" = unname(
                     c(pvals[all_variables], pvals[
                       glue("lag(confirmedTotal, {lag}):",
                            "lag(susceptibleRateTotal, {lag})")]))),
            by="variables")

#### Regional models ####
# Construct formula
fm = paste("testedPositive ~ ",
           glue("lag(testedPositive, {lag}):lag(susceptibleRate, {lag})+"),
           paste(X_regressors, collapse="+")) %>%
  as.formula

for (region in regions){
  # Select only the data for the relevant region
  data = df_long %>% filter(code == !!region)
  
  # Use BIC for model selection
  model = step(lm(fm, data=data), k=log(nrow(data)), trace=0,
               scope=list("lower" = paste("testedPositive ~ ",
                                          glue("lag(testedPositive, {lag}):",
                                               "lag(susceptibleRate, {lag})")) %>%
                            as.formula,
                          "upper" = fm))
  
  png(glue("{output_path}/model1_lag{lag}_lmplot_{region}_bic.png"))
  par(mfrow=c(2,2))
  plot(model)
  par(mfrow=c(1,1))
  dev.off()
  
  # Retrieve parameter estimates
  estimates = coef(summary(model))[, "Estimate"]
  pvals = coef(summary(model))[, "Pr(>|t|)"] # TODO: SE with stars
  
  # Insert parameter estimates in the results table
  results_table_ms = results_table_ms %>%
    left_join(tibble("variables" = c(all_variables, "alpha"),
                     !!glue("{region}") := unname(
                       c(estimates[all_variables], estimates[
                         glue("lag(testedPositive, {lag}):",
                              "lag(susceptibleRate, {lag})")])),
                     !!glue("{region}_pvals") := unname(
                       c(pvals[all_variables], pvals[
                         glue("lag(testedPositive, {lag}):",
                              "lag(susceptibleRate, {lag})")]))),
              by="variables")
}

# Transpose and reorder columns
results_table_ms = results_table_ms %>%
  gather(region, val, 2:ncol(results_table_ms)) %>%
  spread(names(results_table_ms)[1], val) %>%
  select(region, alpha, everything())

# Put the national results at the top
results_table_ms = rbind(
  results_table_ms %>% filter(str_detect(region, "National")),
  results_table_ms %>% filter(!str_detect(region, "National"))) %>%
  column_to_rownames("region")

# Return a LaTeX table
table_ms = xtable(results_table_ms, math.style.exponents = TRUE)
table_ms

#### Plot alpha over time ####
df_meta = readxl::read_xlsx(path_wiki, sheet = "Metadata")

# Starting index - we want at least this number of observations
start = 50

fm = paste("testedPositive ~ ",
           glue("lag(testedPositive, {lag}):lag(susceptibleRate, {lag})+"),
           paste(X_regressors, collapse="+")) %>%
  as.formula

#### Without model selection ####
tbl_alpha = tibble(date = as.Date(NA), alphas=numeric(0), code=character(0))

# Find the estimates of alpha per region over time
for (region in regions){
  alphas = vector("double")
  dates = vector("character")
  
  # Select only the data for the relevant region
  data = df_long %>% filter(code == !!region)
  
  for (t in start:nrow(data)){
    # Estimate the model by OLS
    model = lm(fm, data=data[1:t, ])
    
    # Retrieve the alpha estimate and append this to the list of alphas
    # TODO: Unhardcode the lag
    alpha = model$coefficients[[glue("lag(testedPositive, {lag}):",
                                     "lag(susceptibleRate, {lag})")]]
    alphas = c(alphas, alpha)
  }
  
  # Append the results to the table
  tbl_alpha = tbl_alpha %>%
    bind_rows(tibble(date = data$date[start:nrow(data)],
                     alphas = alphas,
                     code = region) %>%
                drop_na())
}

# Add region and direction to the table
tbl_alpha = tbl_alpha %>%
  left_join(df_meta %>% select(c(region, code, direction)), by="code")

# Make a plot per direction
for (sub_tbl in split(tbl_alpha, tbl_alpha$direction)){
  direc = sub_tbl$direction[1]
  g = ggplot(sub_tbl, aes(date, alphas, colour = region)) + 
    geom_point() +
    geom_smooth(method="loess", span=0.3, se=FALSE) +
    xlab("") +
    ylab(TeX("$\\alpha_{within}$")) +
    scale_colour_manual(values=c("#E69F00", "#56B4E9", "#009E73", "#0072B2",
                                 "#D55E00", "#CC79A7"))
  print(g)
  ggsave(glue("model1_lag{lag}_alphawithin_{direc}.png"), path=output_path)
}

#### With model selection BIC ####
tbl_alpha = tibble(date = as.Date(NA), alphas=numeric(0), code=character(0))

# Find the estimates of alpha per region over time
for (region in regions){
  alphas = vector("double")
  dates = vector("character")
  
  # Select only the data for the relevant region
  data = df_long %>% filter(code == !!region)
  
  for (t in start:nrow(data)){
    # Use BIC for model selection - scope says we want to always keep
    # alpha_within in
    model = step(lm(fm, data=data[1:t, ]), k=log(t), trace=0,
                 scope=list("lower" = paste("testedPositive ~ ",
                                            glue("lag(testedPositive, {lag}):",
                                                 "lag(susceptibleRate, {lag})")) %>%
                              as.formula,
                            "upper" = fm))
    model
    
    # Retrieve the alpha estimate and append this to the list of alphas
    # TODO: Unhardcode the lag
    alpha = model$coefficients[[glue("lag(testedPositive, {lag}):",
                                     "lag(susceptibleRate, {lag})")]]
    alphas = c(alphas, alpha)
  }
  
  # Append the results to the table
  tbl_alpha = tbl_alpha %>%
    bind_rows(tibble(date = data$date[start:nrow(data)],
                     alphas = alphas,
                     code = region) %>%
                drop_na())
}

# Add region and direction to the table
tbl_alpha = tbl_alpha %>%
  left_join(df_meta %>% select(c(region, code, direction)), by="code")

# Make a plot per direction
for (sub_tbl in split(tbl_alpha, tbl_alpha$direction)){
  direc = sub_tbl$direction[1]
  g = ggplot(sub_tbl, aes(date, alphas, colour = region)) + 
    geom_point() +
    geom_smooth(method="loess", span=0.3, se=FALSE) +
    xlab("") +
    ylab(TeX("$\\alpha_{within}$")) +
    scale_colour_manual(values=c("#E69F00", "#56B4E9", "#009E73", "#0072B2",
                                 "#D55E00", "#CC79A7"))
  print(g)
  ggsave(glue("model1_lag{lag}_alphawithin_{direc}_bic.png"), path=output_path)
}
