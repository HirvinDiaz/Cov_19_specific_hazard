#### Creation of data tables ####

####  Clean Work - space  ####
rm(list = ls()) # to clean the work - space

####  Load Packages  ####
library(dplyr)
library(ggplot2)
library(utils)
library(scales)
library(dampack)
library(CEAutil)
library(tidyr)
library(data.table)
library(reshape2)
library(readr)
library(readxl)
library(lubridate)

#### Load data ####
load("data/df_hazards_ICU_hosp.Rdata")

#### Load ratetable ####
load("data/rate_mx.Rdata")

# Public Data base of suspected people with COVID-19 
load("data/cov_09_02.Rdata") 

#### Create df with probabilities ####
d_p_HD <- df_hazards_ICU_hosp %>% 
  filter(Type != "Observed")

d_p_HD <- reshape(d_p_HD,
                  idvar = c("time","group", "sex", "state"),
                  timevar = "Type",
                  direction = "wide")

d_p_HD <- d_p_HD %>% 
  mutate(p_dCoV = 1 - exp(-Hazard.Covid_19)) %>% 
  mutate(p_dPop = 1 - exp(- Hazard.Population)) %>% 
  arrange(group, time)%>% 
  select(c("time", "group", "sex","state", "p_dCoV", "p_dPop"))

d_p_HD_hosp <- d_p_HD %>% 
  filter(state == "Hosp")

d_p_HD_icu <- d_p_HD %>% 
  filter(state == "ICU")

#### Create Cohorts ####
x <- as.Date(max(Cov$date_admission), format = "%Y-%m-%d")

Cohort <- Cov %>% 
  filter(type == 2 | intubated == 1) %>% 
  filter(age >= 45) %>% 
  select(c("sex", "age", "state", "type","date_death","date_symptoms", "intubated"))%>% 
  mutate(sex = ifelse(sex == 1, "female", "male"))%>% 
  mutate(time = ifelse(is.na(date_death),x - date_symptoms, 
                       date_death - date_symptoms)) %>% # time of observation is 
  # date of death minus date of 
  # admission, if censored
  # last date minus doa.
  rename(diag = date_symptoms) %>%  
  mutate(diag = as.numeric(diag, origin = "1960-01-01")) %>%  
  filter(time >= 1) %>% # eliminate negative times
  mutate(age_range = ifelse(age >= 45 & age <= 54, "45 - 54", 
                            ifelse(age >= 55 & age <= 64, "55 - 64",
                                   ifelse(age >= 65 & age <= 69, "65 - 69", 
                                          "70 +"))))

Cohort <- Cohort %>% 
  mutate(state = ifelse(intubated == 1, "ICU", "Hosp")) %>%
  select(c("sex","age_range","state")) %>%
  rename(group = age_range) %>% 
  mutate(time = 0)

Cohort_hosp <- Cohort %>% 
  filter(state == "Hosp")

Cohort_ICU <- Cohort %>% 
  filter(state == "ICU")


#### CEA Hospitalized ####
#### Parameters set ####
n_i              <- length(Cohort_hosp$sex)# number of simulated individuals 
n_t              <- 50 # time horizon, 50 days
v_names_states   <- c("Cov19+", "CoV19_Dead", "O_Causes_Dead") # state names
n_states         <- length(v_names_states)   # number of states
d_c              <- d_e <- 0.0082 # daily discount rates 
# calculate discount weights for costs for each cycle based on discount rate d_c
v_dwc <- 1 / (1 + d_c) ^ (0:n_t) 
# calculate discount weights for effectiveness for each cycle based on discount? 
v_dwe <- 1 / (1 + d_e) ^ (0:n_t)

## Costs and utilities inputs (in MX pesos)
# Average cost for patients that require hospitalized care
c_hosp     <- ((35000 + 50000 + 70000 + 80000)/4)
c_sCov     <- c_hosp + c_amb 
c_dCov     <- 0       # cost of remaining one cycle Dead
c_dPop     <- 0       # cost of remaining one cycle Dead
c_Trt      <- 0
# Mean QALD (Quality Adjusted Life Days) loss.
m_QALD     <-  2.5
u_sCov     <- (100 - m_QALD)/100    # utility when Sick 
u_dCov     <- 0                     # utility when Dead
u_dPop     <- 0                     # utility when Dead

# Convert data frame to a data table for efficiency
dt_p_CoV <- data.table(d_p_HD)

# set the data table to be indexed by age, day, sex and ICU
setkey(dt_p_CoV, group, time, sex, state)

# Create data frame of population from cohort. All begin in day 0
df_X <- Cohort
