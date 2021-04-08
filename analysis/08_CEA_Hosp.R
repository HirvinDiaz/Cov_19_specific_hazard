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

# Load data 
load("data/df_hazards_ICU_hosp.Rdata")

# Load ratetable 
load("data/rate_mx.Rdata")

# Public Data base of suspected people with COVID-19 
load("data/cov_04_04.Rdata") 

# Data with life expectancy
load("data/Life_expectancy.Rdata")

# Load functions for microsimulation
source("R/Functions.R")

#### Filter to get only 2020 ####

Cov <- Cov %>% 
  filter(date_admission >= "2020-03-01" & date_admission <= "2020-12-31")

#### Create df with probabilities ####
d_h_HD_hosp <- df_hazards_ICU_hosp %>% 
  filter(Type != "Observed") %>% 
  filter(state == "Hosp") %>% 
  filter(time <= 50)

d_h_HD_hosp <- reshape(d_h_HD_hosp,
                  idvar = c("time","group", "sex", "state"),
                  timevar = "Type",
                  direction = "wide")

hr_Remd <- 0.73
hr_RemdBa_vs_Remd <- 0.65

d_h_HD_hosp <- d_h_HD_hosp %>% 
  mutate(haz_Remd = exp(log(Hazard.Covid_19) + log(hr_Remd))) %>% 
  mutate(haz_Remd_Ba = exp(log(haz_Remd) + log(hr_RemdBa_vs_Remd)))
  
d_p_HD_hosp <- d_h_HD_hosp %>% 
  mutate(p_dCoV = 1 - exp(-Hazard.Covid_19)) %>% 
  mutate(p_dPop = 1 - exp(- Hazard.Population)) %>% 
  mutate(p_dCoV_Remd = 1 - exp(- haz_Remd)) %>% 
  mutate(p_dCoV_Remd_Ba = 1 - exp(- haz_Remd_Ba)) %>% 
  arrange(group, time)%>% 
  select(c("time", "group", "sex", "p_dCoV", "p_dPop", "p_dCoV_Remd", 
           "p_dCoV_Remd_Ba"))

#### Create Cohorts ####
x <- as.Date("2021-03-29", format = "%Y-%m-%d")

Cohort <- Cov %>% 
  filter(type == 2 & intubated != 1) %>% 
  filter(age >= 45 & age <= 100) %>% 
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
df_pop_ch <- Cohort %>% 
  select(c("sex", "age"))

Cohort <- Cohort %>% 
  mutate(state = ifelse(intubated == 1, "ICU", "Hosp")) %>%
  select(c("sex","age_range","state")) %>%
  rename(group = age_range) %>% 
  mutate(time = 0)

Cohort_hosp <- Cohort %>% 
  mutate(state = "Hosp")

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
v_dwe <- 1 / (1) ^ (0:n_t)

## Costs and utilities inputs (in MX pesos)
# Average cost for patients that require hospitalized care
c_hosp     <- 9272     # cost of remaining one cycle hospitalized 
c_resp     <- 18954.18 # cost of remaining one cycle intubated
c_remd     <- 6188     # cost of Remdesivir Treatment for one Cycle
c_Bari     <- 3672     # cost of Baricitinib Treatment for one Cycle
c_dead     <- 0        # cost of remaining one cycle Dead

# Life Days outcome.
l_sick     <- 1 # utility when Sick 
l_dead     <- 0 # utility when Dead

# Convert data frame to a data table for efficiency
dt_p_CoV <- data.table(d_p_HD_hosp)

# set the data table to be indexed by age, day, sex and ICU
setkey(dt_p_CoV, group, time, sex)

# Create data frame of population from cohort. All begin in day 0
df_X <- Cohort_hosp

Probs <- function(v_M_t, df_X, t, Trt = FALSE) { # t <- 1
  # Arguments:
  # v_M_t: health state occupied at cycle t (character variable)
  # df_X: data frame with individual characteristics data 
  # v_Ts: vector with the duration of being sick
  # t: cycle
  # transition probabilities for that cycle
  # create matrix of state transition probabilities
  m_p_t           <- matrix(0, nrow = n_states, ncol = n_i)  
  # give the state names to the rows
  rownames(m_p_t) <-  v_names_states
  # Lookup baseline probability of dying from Covid-19 or other causes  
  if (Trt == "Remd"){
    p_die_CoV_all <- dt_p_CoV[.(df_X$group, df_X$time + t, df_X$sex), p_dCoV_Remd] 
  } else if (Trt == "Remd_Ba"){
    p_die_CoV_all <- dt_p_CoV[.(df_X$group, df_X$time + t, df_X$sex), p_dCoV_Remd_Ba] 
  } else if (Trt == "None" | Trt == FALSE){
    p_die_CoV_all <- dt_p_CoV[.(df_X$group, df_X$time + t, df_X$sex), p_dCoV] 
  } else {
    message("Choose a treatment within the options")
  }
  p_die_Pop_all <- dt_p_CoV[.(df_X$group, df_X$time + t, df_X$sex), p_dPop] 
  
  p_die_CoV     <- p_die_CoV_all[v_M_t == "Cov19+"]  
  p_die_Pop     <- p_die_Pop_all[v_M_t == "Cov19+"] 
  # update m_p_t with the appropriate probabilities   
  # transition probabilities when healthy
  m_p_t[, v_M_t == "Cov19+"] <- rbind(1 - (p_die_CoV + p_die_Pop), p_die_CoV, p_die_Pop)    
  # transition probabilities when sick 
  m_p_t[, v_M_t == "CoV19_Dead"] <- rbind(0, 1 ,0)
  # transition probabilities when sicker 
  m_p_t[, v_M_t == "O_Causes_Dead"] <- rbind(0, 0, 1)
  return(t(m_p_t))
}     

Costs <- function (v_M_t, Trt = FALSE) {
  # v_M_t: current health state
  c_t <- c()
  if (Trt == "Remd"){
    c_t[v_M_t == "Cov19+"]  <- c_hosp + c_remd
  } else if (Trt == "Remd_Ba"){
    c_t[v_M_t == "Cov19+"]  <- c_hosp + c_remd + c_Bari
  } else if (Trt == "None" | Trt == FALSE){
    c_t[v_M_t == "Cov19+"]  <- c_hosp 
  } else {
    message("Choose a treatment within the options")
  }
  # costs accrued by being Sick this cycle
  c_t[v_M_t == "CoV19_Dead"] <- c_dead  # costs accrued by being sick this cycle
  c_t[v_M_t == "O_Causes_Dead"] <- c_dead  # costs accrued by being sicker this cycle
  return(c_t)  # return costs accrued this cycle
}

# Switch, leer cual es el argumento de treatment.
# Hacerlo con vectores, para utilizarlo con multiplicadores.
# Generar las distribuciones, con mis intervalos de confianza, en el espacio logartmico 
# Derivamos los errores estandar. Saco el lower y el upper bound, 

Effs <- function (v_M_t) {
  # v_M_t: current health state
  q_t <- c() 
  q_t[v_M_t == "Cov19+"] <- l_sick   
  q_t[v_M_t == "CoV19_Dead"] <- l_dead     
  q_t[v_M_t == "O_Causes_Dead"]  <- l_dead     
  return(q_t)  # return the LDs accrued this cycle
}

#### 04.2 Dynamic characteristics 
# These are just starting conditions - they will change with the simulation
v_M_init  <- rep("Cov19+", n_i)       # everyone begins in the sick state

MicroSim <- function(n_i, df_X, df_pop_ch, life_expectancy, Trt = FALSE) { #t <- 1
  # set.seed(02021989)
  m_M <- m_C <- m_E <-  matrix(NA, nrow = n_i, ncol = n_t + 1, 
                               dimnames = list(paste("ind"  , 1:n_i, sep = " "), 
                                               paste("cycle", 0:n_t, sep = " "))) 
  
  m_M[, 1] <- v_M_init          # initial health state
  m_C[, 1] <- Costs(m_M[, 1])   # costs accrued during cycle 0
  m_E[, 1] <- Effs(m_M[, 1])    # QALYs accrued during cycle 0
  
  # open a loop for time running cycles 1 to n_t 
  for (t in 1:n_t) { # t <- 1
    # calculate the transition probabilities for the cycle based on health state t
    m_P <- Probs(v_M_t = m_M[, t], df_X, t = t, Trt = Trt) 
    # sample the current health state and store that state in matrix m_M
    m_M[, t + 1]  <- samplev(m_P, 1)    
    # calculate costs per individual during cycle t + 1
    m_C[, t + 1]  <- Costs(m_M[, t + 1], Trt = Trt)  
    # calculate QALYs per individual during cycle t + 1
    m_E[, t + 1]  <- Effs(m_M[, t + 1])  
    
    # Display simulation progress
    if(t/(n_t/10) == round(t/(n_t/10), 0)) { # display progress every 10%
      cat('\r', paste(t/n_t * 100, "% done", sep = " "))
    }
    
  } # close the loop for the time points 
  
  if (Trt == "Remd"){
    m_C[, 11:n_t + 1]  <- m_C[, 11:n_t + 1] - c_remd
  } else if (Trt == "Remd_Ba"){
    m_C[, 11:n_t + 1]  <- m_C[, 11:n_t + 1] - c_remd
    m_C[, 15:n_t + 1]  <- m_C[, 15:n_t + 1] - c_Bari
  } else {
    m_C <- m_C
  }
  m_C[m_C < 0] <- 0
  
  Final_outcome <- as.data.frame(m_E) 
  Final_outcome <- Final_outcome %>% 
    select(`cycle 50`)
  
  Effects <- cbind(Final_outcome,df_pop_ch)
  
  Effects <- Effects %>% 
    merge(Life_expectancy, 
          by = c("sex" = "sex",
                 "age" = "age")) %>% 
    mutate(years = ifelse(`cycle 50` == 0, 0, Years_dis))
  
  # calculate  
  tc <- m_C %*% v_dwc           # total (discounted) cost per individual
  te <- m_E %*% v_dwe           # total LDs per individual 
  tc_hat <- mean(tc)            # average (discounted) cost 
  te_hat <- mean(te)            # average LDs per individual
  lf_hat <- mean(Effects$years) # average Life years Gained
  tc_sum <- sum(tc)             # sum (discounted) cost
  te_sum <- sum(te)             # sum LDs per individual
  lf_y_s <- sum(Effects$years)  # total Life years Gained
  
  # store the results from the simulation in a list
  results <- list(m_M = m_M, 
                  m_C = m_C, 
                  m_E = m_E, 
                  tc = tc, 
                  te = te, 
                  tc_hat = tc_hat, 
                  te_hat = te_hat,
                  tc_sum = tc_sum, 
                  te_sum = te_sum,
                  Ly_s   = lf_hat)   
  
  return(results)  # return the results
  
} # end of the MicroSim function  

outcomes_remd_ba <- MicroSim(n_i, df_X, df_pop_ch, Life_expectancy, Trt = "Remd_Ba")

results_remd_ba  <- data.frame("Total Cost" = outcomes_remd_ba$tc_hat, 
                         "Total LYs" = outcomes_remd_ba$Ly_s)
results_remd_ba <- results_remd_ba %>% 
  mutate(Cost_eff = Total.Cost/Total.LYs)

#### Test model ####

# Load data 
load("data/df_hazards_ICU_hosp.Rdata")

# Load ratetable 
load("data/rate_mx.Rdata")

d_h_HD_hosp <- d_h_HD_hosp %>% 
  mutate(haz_Remd = exp(log(Hazard.Covid_19) + log(hr_Remd))) %>% 
  mutate(haz_Remd_Ba = exp(log(haz_Remd) + log(hr_RemdBa_vs_Remd)))

d_p_HD_hosp <- d_h_HD_hosp %>% 
  mutate(p_dCoV = 0) %>% 
  mutate(p_dPop = 1 - exp(- Hazard.Population)) %>% 
  mutate(p_dCoV_Remd = 1 - exp(- haz_Remd)) %>% 
  mutate(p_dCoV_Remd_Ba = 1 - exp(- haz_Remd_Ba)) %>% 
  arrange(group, time)%>% 
  select(c("time", "group", "sex", "p_dCoV", "p_dPop", "p_dCoV_Remd", 
           "p_dCoV_Remd_Ba"))

#### Create Cohorts ####
x <- as.Date("2021-03-29", format = "%Y-%m-%d")

Cohort <- Cov %>% 
  filter(type == 2 & intubated != 1) %>% 
  filter(age >= 45 & age <= 100) %>% 
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
df_pop_ch <- Cohort %>% 
  select(c("sex", "age"))

Cohort <- Cohort %>% 
  mutate(state = ifelse(intubated == 1, "ICU", "Hosp")) %>%
  select(c("sex","age_range","state")) %>%
  rename(group = age_range) %>% 
  mutate(time = 0)

Cohort_hosp <- Cohort %>% 
  mutate(state = "Hosp")

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
v_dwe <- 1 / (1) ^ (0:n_t)

## Costs and utilities inputs (in MX pesos)
# Average cost for patients that require hospitalized care
c_hosp     <- 9272     # cost of remaining one cycle hospitalized 
c_resp     <- 18954.18 # cost of remaining one cycle intubated
c_remd     <- 6188     # cost of Remdesivir Treatment for one Cycle
c_Bari     <- 3672     # cost of Baricitinib Treatment for one Cycle
c_dead     <- 0        # cost of remaining one cycle Dead

# Life Days outcome.
l_sick     <- 1 # utility when Sick 
l_dead     <- 0 # utility when Dead

# Convert data frame to a data table for efficiency
dt_p_CoV <- data.table(d_p_HD_hosp)

# set the data table to be indexed by age, day, sex and ICU
setkey(dt_p_CoV, group, time, sex)

# Create data frame of population from cohort. All begin in day 0
df_X <- Cohort_hosp

Probs <- function(v_M_t, df_X, t, Trt = FALSE) { # t <- 1
  # Arguments:
  # v_M_t: health state occupied at cycle t (character variable)
  # df_X: data frame with individual characteristics data 
  # v_Ts: vector with the duration of being sick
  # t: cycle
  # transition probabilities for that cycle
  # create matrix of state transition probabilities
  m_p_t           <- matrix(0, nrow = n_states, ncol = n_i)  
  # give the state names to the rows
  rownames(m_p_t) <-  v_names_states
  # Lookup baseline probability of dying from Covid-19 or other causes  
  if (Trt == "Remd"){
    p_die_CoV_all <- dt_p_CoV[.(df_X$group, df_X$time + t, df_X$sex), p_dCoV_Remd] 
  } else if (Trt == "Remd_Ba"){
    p_die_CoV_all <- dt_p_CoV[.(df_X$group, df_X$time + t, df_X$sex), p_dCoV_Remd_Ba] 
  } else if (Trt == "None" | Trt == FALSE){
    p_die_CoV_all <- dt_p_CoV[.(df_X$group, df_X$time + t, df_X$sex), p_dCoV] 
  } else {
    message("Choose a treatment within the options")
  }
  p_die_Pop_all <- dt_p_CoV[.(df_X$group, df_X$time + t, df_X$sex), p_dPop] 
  
  p_die_CoV     <- p_die_CoV_all[v_M_t == "Cov19+"]  
  p_die_Pop     <- p_die_Pop_all[v_M_t == "Cov19+"] 
  # update m_p_t with the appropriate probabilities   
  # transition probabilities when healthy
  m_p_t[, v_M_t == "Cov19+"] <- rbind(1 - (p_die_CoV + p_die_Pop), p_die_CoV, p_die_Pop)    
  # transition probabilities when sick 
  m_p_t[, v_M_t == "CoV19_Dead"] <- rbind(0, 1 ,0)
  # transition probabilities when sicker 
  m_p_t[, v_M_t == "O_Causes_Dead"] <- rbind(0, 0, 1)
  return(t(m_p_t))
}     

Costs <- function (v_M_t, Trt = FALSE) {
  # v_M_t: current health state
  c_t <- c()
  if (Trt == "Remd"){
    c_t[v_M_t == "Cov19+"]  <- c_hosp + c_remd
  } else if (Trt == "Remd_Ba"){
    c_t[v_M_t == "Cov19+"]  <- c_hosp + c_remd + c_Bari
  } else if (Trt == "None" | Trt == FALSE){
    c_t[v_M_t == "Cov19+"]  <- c_hosp 
  } else {
    message("Choose a treatment within the options")
  }
  # costs accrued by being Sick this cycle
  c_t[v_M_t == "CoV19_Dead"] <- c_dead  # costs accrued by being sick this cycle
  c_t[v_M_t == "O_Causes_Dead"] <- c_dead  # costs accrued by being sicker this cycle
  return(c_t)  # return costs accrued this cycle
}

# Switch, leer cual es el argumento de treatment.
# Hacerlo con vectores, para utilizarlo con multiplicadores.
# Generar las distribuciones, con mis intervalos de confianza, en el espacio logartmico 
# Derivamos los errores estandar. Saco el lower y el upper bound, 

Effs <- function (v_M_t) {
  # v_M_t: current health state
  q_t <- c() 
  q_t[v_M_t == "Cov19+"] <- l_sick   
  q_t[v_M_t == "CoV19_Dead"] <- l_dead     
  q_t[v_M_t == "O_Causes_Dead"]  <- l_dead     
  return(q_t)  # return the LDs accrued this cycle
}

#### 04.2 Dynamic characteristics 
# These are just starting conditions - they will change with the simulation
v_M_init  <- rep("Cov19+", n_i)       # everyone begins in the sick state

MicroSim <- function(n_i, df_X, df_pop_ch, life_expectancy, Trt = FALSE) { #t <- 1
  # set.seed(02021989)
  m_M <- m_C <- m_E <-  matrix(NA, nrow = n_i, ncol = n_t + 1, 
                               dimnames = list(paste("ind"  , 1:n_i, sep = " "), 
                                               paste("cycle", 0:n_t, sep = " "))) 
  
  m_M[, 1] <- v_M_init          # initial health state
  m_C[, 1] <- Costs(m_M[, 1])   # costs accrued during cycle 0
  m_E[, 1] <- Effs(m_M[, 1])    # QALYs accrued during cycle 0
  
  # open a loop for time running cycles 1 to n_t 
  for (t in 1:n_t) { # t <- 1
    # calculate the transition probabilities for the cycle based on health state t
    m_P <- Probs(v_M_t = m_M[, t], df_X, t = t, Trt = Trt) 
    # sample the current health state and store that state in matrix m_M
    m_M[, t + 1]  <- samplev(m_P, 1)    
    # calculate costs per individual during cycle t + 1
    m_C[, t + 1]  <- Costs(m_M[, t + 1], Trt = Trt)  
    # calculate QALYs per individual during cycle t + 1
    m_E[, t + 1]  <- Effs(m_M[, t + 1])  
    
    # Display simulation progress
    if(t/(n_t/10) == round(t/(n_t/10), 0)) { # display progress every 10%
      cat('\r', paste(t/n_t * 100, "% done", sep = " "))
    }
    
  } # close the loop for the time points 
  
  if (Trt == "Remd"){
    m_C[, 11:n_t + 1]  <- m_C[, 11:n_t + 1] - c_remd
  } else if (Trt == "Remd_Ba"){
    m_C[, 11:n_t + 1]  <- m_C[, 11:n_t + 1] - c_remd
    m_C[, 15:n_t + 1]  <- m_C[, 15:n_t + 1] - c_Bari
  } else {
    m_C <- m_C
  }
  m_C[m_C < 0] <- 0
  
  Final_outcome <- as.data.frame(m_E) 
  Final_outcome <- Final_outcome %>% 
    select(`cycle 50`)
  
  Effects <- cbind(Final_outcome,df_pop_ch)
  
  Effects <- Effects %>% 
    merge(Life_expectancy, 
          by = c("sex" = "sex",
                 "age" = "age")) %>% 
    mutate(years = ifelse(`cycle 50` == 0, 0, Years_dis))
  
  # calculate  
  tc <- m_C %*% v_dwc           # total (discounted) cost per individual
  te <- m_E %*% v_dwe           # total LDs per individual 
  tc_hat <- mean(tc)            # average (discounted) cost 
  te_hat <- mean(te)            # average LDs per individual
  lf_hat <- mean(Effects$years) # average Life years Gained
  tc_sum <- sum(tc)             # sum (discounted) cost
  te_sum <- sum(te)             # sum LDs per individual
  lf_y_s <- sum(Effects$years)  # total Life years Gained
  
  # store the results from the simulation in a list
  results <- list(m_M = m_M, 
                  m_C = m_C, 
                  m_E = m_E, 
                  tc = tc, 
                  te = te, 
                  tc_hat = tc_hat, 
                  te_hat = te_hat,
                  tc_sum = tc_sum, 
                  te_sum = te_sum,
                  Ly_s   = lf_hat)   
  
  return(results)  # return the results
  
} # end of the MicroSim function  

outcomes_0_prob <- MicroSim(n_i, df_X, df_pop_ch, Life_expectancy, Trt = F)

results_0_prob  <- data.frame("Total Cost" = outcomes_0_prob$tc_hat, 
                               "Total LYs" = outcomes_0_prob$Ly_s)
results_0_prob <- results_0_prob %>% 
  mutate(Cost_eff = Total.Cost/Total.LYs)

mean(Life_expectancy$Years_dis)

Cov <- Cov %>% 
  left_join(Life_expectancy, by = "age")
df_pop_ch <- df_pop_ch %>% 
  left_join(Life_expectancy, by = "age")

mean(df_pop_ch$Years_dis)
