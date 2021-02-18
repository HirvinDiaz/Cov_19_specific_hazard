#### Final project script. Author: Hirvin Azael Diaz Zepeda ####
#### Relative survival over time ####

####  Clean Work - space  ####
rm(list = ls()) # to clean the work - space

####  Load Packages  ####
library(dplyr)
library(ggplot2)
library(utils)
library(scales)
library(dampack)
library(chron)
library(tibble)
library(knitr)
library(pander)
library(kableExtra)
library(CEAutil)
library(tidyr)
library(dampack)
library(data.table)
library(reshape2)
library(survival)
library(survPen)
library(relsurv)
library(readr)
library(readxl)
library(lubridate)
library(muhaz)
library(ggpubr)

######  Load data bases  ######

# Data base of suspected people with COVID-19 published by Mexico´s Ministry of 
# Health 

load("data/cov_24_01.Rdata") 

#### Load ratetable ####
load("data/rate_mx.Rdata")

#### Covid_rr_hos_res for intubated individuals higher than 44 years in national DB ####
x <- as.Date(max(Cov$date_admission), format = "%Y-%m-%d")

Covid_rr <- Cov %>% 
  filter(intubated == 1 | type == 2) %>% 
  filter(age >= 45) %>% 
  select(c("sex", "age", "date_admission", "date_death", "date_symptoms", "state"))%>% 
  mutate(sex = ifelse(sex == 1, 2, 1))%>% # recode 1 is male and 2 female
  mutate(time = ifelse(is.na(date_death),x - date_symptoms, 
                       date_death - date_symptoms)) %>% # time of observation is 
  mutate(diag = date_symptoms) %>%                 # date of death minus date of 
  mutate(diag = as.numeric(diag, origin = "1960-01-01")) %>%
  filter(time >= 1) %>% # eliminate negative times
  mutate(original_age = age) %>% 
  mutate(age = original_age*365.25) %>%  # age needs to be in days for relsurv
  mutate(age_range = ifelse(original_age >= 45 & original_age <= 54, "45 - 54", 
                            ifelse(original_age >= 55 & original_age <= 64, "55 - 64",
                                   ifelse(original_age >= 65 & original_age <= 69, 
                                          "65 - 69", "70 +")))) %>% 
  mutate(month = ifelse(date_symptoms > "2020-01-01" & date_symptoms < "2020-04-01", "March", 
                        ifelse(date_symptoms >= "2020-04-01" & date_symptoms <= "2020-04-30", "April", 
                               ifelse(date_symptoms >= "2020-05-01" & date_symptoms <= "2020-05-31", "May",
                                      ifelse(date_symptoms >= "2020-06-01" & date_symptoms <= "2020-06-30", "June",
                                             ifelse(date_symptoms >= "2020-07-01" & date_symptoms <= "2020-07-31", "July",
                                                    ifelse(date_symptoms >= "2020-08-01" & date_symptoms <= "2020-08-31", "August",
                                                           ifelse(date_symptoms >= "2020-09-01" & date_symptoms <= "2020-09-30", "September",
                                                                  ifelse(date_symptoms >= "2020-10-01" & date_symptoms <= "2020-10-31", "October",
                                                                         ifelse(date_symptoms >= "2020-11-01" & date_symptoms <= "2020-11-30", "November",
                                                                                ifelse(date_symptoms >= "2020-12-01" & date_symptoms <= "2020-12-31", "December", "January 21")))))))))))


#### Apply rs.surv function
KM_surv_sex_age <- rs.surv.exp(Surv(time, state) ~ sex,
                               data = Covid_rr,
                               ratetable = rate_exp_mx, 
                               type="kaplan-meier",
                               conf.type="log",
                               conf.int=0.95,
                               rmap = list(age = age, 
                                           sex = sex, 
                                           year = diag))

plot(KM_surv_sex_age)
