#### 10 Deterministic Sensitibity Analysis  ####
# Load packages
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

# Load data 
load("data/df_hazards_ICU_hosp.Rdata")

# Load ratetable 
load("data/rate_mx.Rdata")

# Public Data base of suspected people with COVID-19 
load("data/cov_09_02.Rdata") 

# Load functions
source("R/Functions_Microsim.R")
source("R/Functions.R")

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
  mutate(state = "Hosp")

#### Create df with probabilities ####
d_h_HD_hosp <- df_hazards_ICU_hosp %>% 
  filter(Type != "Observed") %>% 
  filter(state == "Hosp")

d_h_HD_hosp <- reshape(d_h_HD_hosp,
                       idvar = c("time","group", "sex", "state"),
                       timevar = "Type",
                       direction = "wide")
# Create data frame of population from cohort. All begin in day 0
df_X <- Cohort_hosp

# Erase data frames not useful in this point
rm(Cohort,Cohort_hosp,Cov)

# List of parameters
l_params <- list( #El elemento que estamos creando es una lista
  Remd_effect    = 0.73,  # probability to die when healthy
  Remd_Ba_effect = 0.65,   # probability to become sick when healthy
  c_hosp         = 8141.63,    # probability to become healthy when sick
  c_remd         = 8331.85,  # probability to become sicker when sick
  c_bari         = 3672.35,      # hazard ratio of death in sick vs healthy
  c_dead         = 0,     # hazard ratio of death in sicker vs healthy
  l_sick         = 1,   # cost of remaining one cycle in the healthy state
  l_dead         = 0   # cost of remaining one cycle in the sick state
)

# Strategy names
v_names_str <- c("No Treatment", "Remdesivir", "Remdesivir & Bariticinib") 

# Number of strategies
n_str <- length(v_names_str)

CEA_function(l_params = l_params, df_X = df_X, df_h = d_h_HD_hosp)

options(scipen = 999) # disabling scientific notation in R

c_hospit <- 8141.63

c_hospit_max <- 8141.63*1.25
c_hospit_min <- 8141.63*0.75
c_remd_max <- 8331.85*1.25
c_remd_min <- 8331.85*0.75

# dataframe containing all parameters, their basecase values, and the min and 
# max values of the parameters of interest 
df_params_owsa <- data.frame(pars = c("Remd_effect", "Remd_Ba_effect", "c_hosp", "c_remd"),
                             min  = c(0.52 ,  0.39 , 6106.2225, 6248.8875), # min parameter values
                             max  = c(1.03, 1.09 , 10177.0375, 10414.8125)  # max parameter values
)

owsa_nmb  <- run_owsa_det(params_range     = df_params_owsa,   # dataframe with parameters for owsa
                          params_basecase  = l_params,     # list with all parameters
                          nsamp            = 100,              # number of parameter values
                          # Entre por ejemplo .05 y 0.1555
                          FUN              = CEA_function(l_params, df_X, d_h_HD_hosp), # function to compute outputs
                          outcomes         = c("NMB"),         # output to do the OWSA on
                          strategies       = v_names_str,      # names of the strategies
                          n_wtp            = 120000)           # extra argument to pass to FUN



