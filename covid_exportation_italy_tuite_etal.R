######################################################################################
### backcalculating number of infections in source country based on exported cases ###
###                 formula based on Fraser et al., Science 2009                   ###
######################################################################################

### Use data on exported COVID cases to estimate incidence in source country

source.country <- "Italy"

### read in travel data for 2015

trav.dat <- read_delim("data/PHAC1502.csv", 
                       ";", escape_double = FALSE, trim_ws = TRUE) %>%
            select(-`1`)

names(trav.dat)<-c("Period","Origin","Leg1","Leg2","Leg3","Leg4","Leg5","Destination","MC1","MC2","MC3","MC4","MC5","MC6","Classe","Pax","Estimated_Pax","Average_Fare")

locations.dat <- read_delim("data/locationsregions_modif.csv", 
                            ";", escape_double = FALSE, trim_ws = TRUE) %>%
                select(loc, name, ctry, ctryname) %>%
                distinct(loc, .keep_all = TRUE)


trav.dat %>%
  select(Origin, Destination, Estimated_Pax) %>%
  left_join(locations.dat, by=c("Origin"="loc")) %>%
  rename(origin_ctry = ctryname,
         origin_name = name) %>%
  left_join(locations.dat, by=c("Destination"="loc")) %>%
  rename(dest_ctry = ctryname,
         dest_name = name) %>%
  filter(origin_ctry=="Italy") %>%
  group_by(dest_ctry) %>%
  summarize(paxTot = sum(Estimated_Pax)) %>%
  arrange(-paxTot) %>%
  filter(dest_ctry!="Italy") -> it.trav.dat
  
paxTot2015 <- sum(it.trav.dat$paxTot)


### 2019 IATA travel data
trav.dat <- read_excel("data/final-destination-city-volume - Italy to the World - February 2019.xlsx") %>%
            filter(!is.na(Country))

paxTot2019 <- sum(trav.dat$`Total volume`)


export.data <- read_excel("data/exported_cases_from_italy.xlsx") 

export.dat <-  export.data %>% 
               arrange(country) %>%
               pull(country)
excl.overland.dat <- export.data %>%
                filter(shared_border_or_documented_overland==0) %>%
                pull(country)

export.data %>%
  filter(shared_border_or_documented_overland==1)
export.data %>%
  group_by(shared_border_or_documented_overland) %>%
  summarize(expCases = sum(exported_cases)) %>%
  pull(expCases) -> exported.cases


#####

###DEFINE PARAMETERS

inbound_vol <- 48576 
outbound_vol <- 28460 
ratio <- inbound_vol/ (inbound_vol + outbound_vol)

T_v_in <- 3.4/30 
N_r = 60.48e6  #population of Italy
T_r = 1# outbreak duration (in months) -- assume that cases occurred in past month
T_v <- (ratio * T_v_in) + (1-ratio)*T_r #average time at risk for exported cases (in months)

estimate.size <- function(include_overland, adjust_to_2019, include_all_countries){
  
  inf_v <- if_else(include_overland==TRUE, sum(exported.cases), exported.cases[1])  # 40 excluding countries with shared borders/documented overland travel; 46 if include all
  
  if(include_all_countries==TRUE){
    it.trav.dat %>% 
      summarize(tot = sum(paxTot)) %>%
      pull() -> N_v
  } else {
  it.trav.dat %>%  #total outbound passengers to countries with exported cases (excluding overland/shared borders)
    filter(if (include_overland==TRUE) dest_ctry %in% export.dat else dest_ctry %in% excl.overland.dat) %>%
    summarize(tot = sum(paxTot)) %>%
    pull() -> N_v
  }
  
  N_v <- if_else(adjust_to_2019==TRUE, N_v*paxTot2019/paxTot2015, N_v)
  
  lambda = inf_v / (N_v * T_v) # average monthly rate of infection
  
  lambda.params<-poisson.test(x=inf_v, T=(N_v * T_v))
  lambdas <- c(lambda_lower = lambda.params$conf.int[1], 
               lambda_mean = lambda.params$estimate,
               lambda_upper = lambda.params$conf.int[2])
  
  ### calculate range of expected number of (cumulative) infections in source population ###
  ###               using lower, mean, and upper bounds for lambda in travelers          ###
  
  X_r = NULL
  
  for(i in 1:length(lambdas)){
    X_r[i] = lambdas[i] * N_r * T_r
  }
  
  return(X_r)
}



#######################

x1<-estimate.size(include_overland=TRUE, adjust_to_2019=FALSE, include_all_countries = FALSE)
x2 <- estimate.size(include_overland=TRUE, adjust_to_2019=TRUE, include_all_countries = FALSE)
x3<- estimate.size(include_overland=FALSE, adjust_to_2019=FALSE, include_all_countries = FALSE)
x4 <- estimate.size(include_overland=FALSE, adjust_to_2019=TRUE, include_all_countries = FALSE)
x5 <- estimate.size(include_overland=FALSE, adjust_to_2019=FALSE, include_all_countries = TRUE)
x6 <- estimate.size(include_overland=FALSE, adjust_to_2019=TRUE, include_all_countries = TRUE)
data.frame(rbind(x1, x2, x3, x4, x5, x6)) -> est.dat


est.dat %>%
  rename(LCI = X1, 
         mean = X2, 
         UCI = X3) %>%
  rownames_to_column(var="scenario") %>%
  mutate(scenario = factor(scenario, levels = c("x1", "x2", "x3", "x4", "x5", "x6"), 
                           labels = c("All cases", 
                                      "All cases (adjusted to 2019)",
                                      "Exclude bordering countries and overland travel", 
                                      "Exclude bordering countries and overland travel (adjusted to 2019)", 
                                      "Include travel to all countries",
                                      "Include travel to all countries (adjusted to 2019)"
                                      )), 
         basecase = if_else(scenario=="All cases", "TRUE", "FALSE")) -> est.data


ggplot(data=est.data, aes(x = Scenario, y= mean, ymin=LCI, ymax=UCI, color= basecase)) +
           geom_pointrange(size = 2) +
           coord_flip() +
           theme_half_open() +
           background_grid() +
           theme(text = element_text(size=24), 
                 legend.position = "none") +
           scale_color_manual(values = c("grey", "blue"))+
           labs(x = "Exported case locations\n", y="\nEstimated outbreak size in Italy") ->p
            




