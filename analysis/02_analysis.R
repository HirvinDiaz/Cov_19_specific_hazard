#### Final project script. Author: Hirvin Azael Diaz Zepeda ####

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

load("data/cov_14_01.Rdata") 

##### Disease-specific hazard #####

# Data base with national mortality data from repo of PADECI demog-mx created 
# from CONAPO data bases
df_mortrate_state_age_sex <- read_csv("data-raw/df_mortrate_state_age_sex.csv")

# Create Mortality rates for Mexico
# Data frame with population mortality rates in 2020 and 2021 at a national level 
# by sex and age

Mort_mx_2020 <- df_mortrate_state_age_sex %>% 
  filter(year == 2020) %>% 
  group_by( age, sex) %>% 
  summarise(mort_rate = mean(mort_rate)) %>% 
  filter(age <= 100) 

Mort_mx_2020 <- as.data.frame(Mort_mx_2020)

#
Mort_mx_2021 <- df_mortrate_state_age_sex %>% 
  filter(year == 2021) %>% 
  group_by(age, sex) %>% 
  summarise(mort_rate = mean(mort_rate)) %>% 
  filter(age <= 100) 

Mort_mx_2021 <- as.data.frame(Mort_mx_2021)

# Reshape data to have every row by sex and age
Mort_mx_wide_2020 <- reshape(data = Mort_mx, 
                        idvar = "age", 
                        v.names = "mort_rate",
                        timevar = "sex",
                        direction = "wide")

Mort_mx_wide_2020 <- Mort_mx_wide_2020 %>% 
  rename(female = mort_rate.Female,
         male = mort_rate.Male)

#
Mort_mx_wide_2021 <- reshape(data = Mort_mx, 
                             idvar = "age", 
                             v.names = "mort_rate",
                             timevar = "sex",
                             direction = "wide")

Mort_mx_wide_2021 <- Mort_mx_wide_2021 %>% 
  rename(female = mort_rate.Female,
         male = mort_rate.Male)


#### Population survival rates
male_mex <- cbind(exp(- Mort_mx_wide_2020$male), exp(- Mort_mx_wide_2021$male))
female_mex <- cbind(exp(- Mort_mx_wide_2020$female), exp(- Mort_mx_wide_2021$female))

# Create rate data frame in the format relsurv needs
rate_exp_mx <- transrate(men =  male_mex, 
                         women = female_mex,
                         yearlim = c(2020,2021), 
                         int.length = 1)
rate_exp_mx

save(rate_exp_mx, file = "data/rate_mx.Rdata")

# Create database of the cohort in the format of relsurv
x <- as.Date(max(Cov$date_admission), format = "%Y-%m-%d")

#Covid_rr <- Cov %>% 
#  filter(age >= 45) %>% 
#  select(c("sex", "age", "date_admission", "date_death", "date_symptoms", "state"))%>% 
#  mutate(sex = ifelse(sex == 1, 2, 1))%>% # recode 1 is male and 2 female
#  mutate(time = ifelse(is.na(date_death),x - date_admission, 
#                       date_death - date_admission)) %>% # time of observation is 
#  rename(diag = date_symptoms) %>%                       # date of death minus date of 
#  # admission, if censored
#  # last date minus doa.
#  mutate(diag = 
#           as.numeric(diag, origin = "1960-01-01")) %>%  # diagnostic date = date_symptoms
#  filter(time >= 1) %>% # eliminate negative times
#  mutate(original_age = age) %>% 
#  mutate(age = original_age*365.25) %>%  # age needs to be in days for relsurv
#  mutate(age_range = ifelse(original_age >= 45 & original_age <= 54, "45 - 54", 
#                            ifelse(original_age >= 55 & original_age <= 64, "55 - 64",
#                                   ifelse(original_age >= 65 & original_age <= 69, 
#                                          "65 - 69", "70 +"))))

#### Kaplan-Meier with date_admission as initial date ####

# KM_surv <- rs.surv(Surv(time, state) ~ 1,
#         data = Covid_rr,
#         ratetable = rate_exp_mx_2020, 
#         type="kaplan-meier",
#         conf.type="log",
#         conf.int=0.95,
#         rmap = list(age = age, sex = sex, year = diag))
# 
# plot(KM_surv, 
#      conf.int = FALSE, 
#      xlim = c(0,75), 
#      ylim = c(0.80,1))


#### Covid_rr for all individuals higher than 44 years in national DB ####

Covid_rr <- Cov %>% 
#  filter(intubated == 1) %>% 
  filter(age >= 45) %>% 
  select(c("sex", "age", "date_admission", "date_death", "date_symptoms", "state"))%>% 
  mutate(sex = ifelse(sex == 1, 2, 1))%>% # recode 1 is male and 2 female
  mutate(time = ifelse(is.na(date_death),x - date_symptoms, 
                       date_death - date_symptoms)) %>% # time of observation is 
  rename(diag = date_symptoms) %>%                 # date of death minus date of 
  # admission, if censored
  # last date minus doa.
  mutate(diag = 
           as.numeric(diag, origin = "1960-01-01")) %>%  # diagnostic date = date_symptoms
  filter(time >= 1) %>% # eliminate negative times
  mutate(original_age = age) %>% 
  mutate(age = original_age*365.25) %>%  # age needs to be in days for relsurv
  mutate(age_range = ifelse(original_age >= 45 & original_age <= 54, "45 - 54", 
                            ifelse(original_age >= 55 & original_age <= 64, "55 - 64",
                                   ifelse(original_age >= 65 & original_age <= 69, 
                                          "65 - 69", "70 +"))))

KM_surv_sym <- rs.surv(Surv(time, state) ~ sex,
                   data = Covid_rr,
                   ratetable = rate_exp_mx, 
                   type="kaplan-meier",
                   conf.type="log",
                   conf.int=0.95,
                   rmap = list(age = age, sex = sex, year = diag))

KM_surv_n <- head(KM_surv$surv, 61)
KM_surv_sym_n <- head(KM_surv_sym$surv, 61)
time <- seq(0,60, by = 1)

plot(KM_surv_sym)

df_survival <- data.frame(time, KM_surv_n, KM_surv_sym_n)

df_survival_long <- gather(data = df_survival, 
                           key = "Survival", 
                           value = "Prob",
                           -time)
ggplot(df_survival_long, 
       aes(x = time, 
           y = Prob)) +
  geom_line() + 
  facet_grid(~Survival)

# a matrix containing the yearly (conditional) probabilities of one year survival for men.
# Rows represent age (increasing 1 year per line,starting with 0), 
# the columns represent cohort years (the limits are in yearlim, the increase is in int.length.

#### Hospitalized ####

Cov_hosp <- Cov %>% 
  filter(type == 2)

Covid_rr <- Cov_hosp %>% 
  filter(age >= 45) %>% 
  select(c("sex", "age", "date_admission", "date_death", "date_symptoms", "state"))%>% 
  mutate(sex = ifelse(sex == 1, 2, 1))%>% # recode 1 is male and 2 female
  mutate(time = ifelse(is.na(date_death),x - date_symptoms, 
                       date_death - date_symptoms)) %>% # time of observation is 
  rename(diag = date_symptoms) %>%                       # date of death minus date of 
  # admission, if censored
  # last date minus doa.
  mutate(diag = 
           as.numeric(diag, origin = "1960-01-01")) %>%  # diagnostic date = date_symptoms
  filter(time >= 1) %>% # eliminate negative times
  mutate(original_age = age) %>% 
  mutate(age = original_age*365.25) %>%  # age needs to be in days for relsurv
  mutate(age_range = ifelse(original_age >= 45 & original_age <= 54, "45 - 54", 
                            ifelse(original_age >= 55 & original_age <= 64, "55 - 64",
                                   ifelse(original_age >= 65 & original_age <= 69, 
                                          "65 - 69", "70 +"))))

KM_surv_hosp <- rs.surv(Surv(time, state) ~ 1,
                       data = Covid_rr,
                       ratetable = rate_exp_mx_2020, 
                       type="kaplan-meier",
                       conf.type="log",
                       conf.int=0.95,
                       rmap = list(age = age, sex = sex, year = diag))

KM_surv_hosp_n <- head(KM_surv_hosp$surv, 61)

df_survival <- data.frame(time, KM_surv_n, KM_surv_sym_n, KM_surv_hosp_n)

df_survival_long <- gather(data = df_survival, 
                           key = "Survival", 
                           value = "Prob",
                           -time)

ggplot(df_survival_long, 
       aes(x = time, 
           y = Prob)) +
  geom_line() + 
  facet_grid(~Survival)

#### ICU ####
Cov_icu <- Cov %>% 
  filter(icu == 1)

Covid_rr <- Cov_icu %>% 
  filter(age >= 45) %>% 
  select(c("sex", "age", "date_admission", "date_death", "date_symptoms", "state"))%>% 
  mutate(sex = ifelse(sex == 1, 2, 1))%>% # recode 1 is male and 2 female
  mutate(time = ifelse(is.na(date_death),x - date_symptoms, 
                       date_death - date_symptoms)) %>% # time of observation is 
  rename(diag = date_symptoms) %>%                       # date of death minus date of 
  # admission, if censored
  # last date minus doa.
  mutate(diag = 
           as.numeric(diag, origin = "1960-01-01")) %>%  # diagnostic date = date_symptoms
  filter(time >= 1) %>% # eliminate negative times
  mutate(original_age = age) %>% 
  mutate(age = original_age*365.25) %>%  # age needs to be in days for relsurv
  mutate(age_range = ifelse(original_age >= 45 & original_age <= 54, "45 - 54", 
                            ifelse(original_age >= 55 & original_age <= 64, "55 - 64",
                                   ifelse(original_age >= 65 & original_age <= 69, 
                                          "65 - 69", "70 +"))))

KM_surv_icu <- rs.surv(Surv(time, state) ~ 1,
                        data = Covid_rr,
                        ratetable = rate_exp_mx_2020, 
                        type="kaplan-meier",
                        conf.type="log",
                        conf.int=0.95,
                        rmap = list(age = age, sex = sex, year = diag))

KM_surv_icu_n <- head(KM_surv_icu$surv, 61)

df_survival <- data.frame(time, KM_surv_n, KM_surv_sym_n, KM_surv_hosp_n, KM_surv_icu_n)

df_survival_long <- gather(data = df_survival, 
                           key = "Survival", 
                           value = "Prob",
                           -time)


data(rdata)

View(rdata)
ggplot(df_survival_long, 
       aes(x = time, 
           y = Prob)) +
  geom_line() + 
  facet_grid(~Survival)

