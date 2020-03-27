rm(list=ls())

library(rvest)
library(tidyverse)

url_nl = "https://en.wikipedia.org/wiki/2020_coronavirus_pandemic_in_the_Netherlands"
xpath_nl = '//*[@id="mw-content-text"]/div/table[4]' # New COVID-19 cases in the Netherlands by province

#### Read the data ####
scrape_table = function(url, xpath) {
  tbl = url %>%
    read_html() %>%
    html_nodes(xpath=xpath) %>%
    html_table(fill=T) %>%
    `[[`(1)
  
  return(tbl)
}

table_nl = scrape_table(url_nl, xpath_nl)

# Set column names to Date, province names, Unknown, Total
colnames(table_nl) = table_nl[1, ]
colnames(table_nl) = c("Date", colnames(table_nl)[
  sapply(colnames(table_nl), function(x){nchar(x) == 2})], "Unknown", "Total")

# Subset the data to keep only the province-level data and total
table_nl = table_nl[2:(nrow(table_nl)-2), 1:15] # 1 date column, 12 provinces, 1 unknown, 1 total

#### Clean the data ####
## Step 1: Remove the source indicators from Wikipedia: [e], for instance
table_nl = data.frame(lapply(table_nl,
                             function(x) { str_replace(x, "\\[\\w+\\]", "") }))

## Step 2: Change the data type from character to numeric
table_nl[2:ncol(table_nl)] = lapply(table_nl[2:ncol(table_nl)],
                                    function(x) as.numeric(as.character(x)))

## Step 3: Change the data type for the Date column from character to Date
table_nl$Date = as.Date.factor(table_nl$Date, "%Y-%m-%d")
