#### CEA intubated ####

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
load("data/cov_09_02.Rdata")

# Data with life expectancy
load("data/Life_expectancy.Rdata")

# Load functions for microsimulation
source("R/Functions.R")

#### Create df with probabilities ####
d_h_HD_resp <- df_hazards_ICU_hosp %>% 
  filter(Type != "Observed") %>% 
  filter(state == "ICU") %>% 
  filter(time <= 50)

d_h_HD_resp <- reshape(d_h_HD_resp,
                       idvar = c("time","group", "sex", "state"),
                       timevar = "Type",
                       direction = "wide")

hr_Dexa <- 0.64

d_h_HD_resp <- d_h_HD_resp %>% 
  mutate(haz_Dexa = exp(log(Hazard.Covid_19) + log(hr_Dexa)))

d_p_HD_resp <- d_h_HD_resp %>% 
  mutate(p_dCoV = 1 - exp(-Hazard.Covid_19)) %>% 
  mutate(p_dPop = 1 - exp(- Hazard.Population)) %>% 
  mutate(p_dCoV_Dexa = 1 - exp(- haz_Dexa)) %>%  
  arrange(group, time)%>% 
  select(c("time", "group", "sex", "p_dCoV", "p_dPop", "p_dCoV_Dexa"))

#### Create Cohorts ####
x <- as.Date(max(Cov$date_admission), format = "%Y-%m-%d")

Cohort <- Cov %>% 
  filter(intubated == 1) %>% 
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
df_pop_ch <- Cohort %>% 
  select(c("sex", "age"))

Cohort <- Cohort %>% 
  mutate(state = ifelse(intubated == 1, "ICU", "Hosp")) %>%
  select(c("sex","age_range","state")) %>%
  rename(group = age_range) %>% 
  mutate(time = 0)

Cohort_ICU <- Cohort %>% 
  filter(state == "ICU")

#### CEA ICU ####
#### Parameters set ####
n_i              <- length(Cohort_ICU$sex)# number of simulated individuals 
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
c_hosp     <- 16550    # cost of remaining one cycle hospitalized 
c_resp     <- 73542.50 # cost of remaining one cycle intubated
c_dexa     <- 30       # cost of Dexamethasone Treatment for one Cycle
c_dead     <- 0        # cost of remaining one cycle Dead

# Life Days outcome.
l_sick     <- 1 # utility when Sick 
l_dead     <- 0 # utility when Dead

# Convert data frame to a data table for efficiency
dt_p_CoV <- data.table(d_p_HD_resp)

# set the data table to be indexed by age, day, sex and ICU
setkey(dt_p_CoV, group, time, sex)

# Create data frame of population from cohort. All begin in day 0
df_X <- Cohort_ICU

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
  if (Trt == "Dexa"){
    p_die_CoV_all <- dt_p_CoV[.(df_X$group, df_X$time + t, df_X$sex), p_dCoV_Dexa] 
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
  if (Trt == "Dexa"){
    c_t[v_M_t == "Cov19+"]  <- c_resp + c_dexa
  } else if (Trt == "None" | Trt == FALSE){
    c_t[v_M_t == "Cov19+"]  <- c_resp 
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
  set.seed(02021989)
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
  
  if (Trt == "Dexa"){
    m_C[, 11:n_t + 1]  <- m_C[, 11:n_t + 1] - c_dexa
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
    mutate(years = ifelse(`cycle 50` == 0, 0, Years))
  
  # calculate  
  tc <- m_C %*% v_dwc    # total (discounted) cost per individual
  te <- m_E %*% v_dwe    # total LDs per individual 
  tc_hat <- mean(tc)     # average (discounted) cost 
  te_hat <- mean(te)     # average LDs per individual
  tc_sum <- sum(tc)      # sum (discounted) cost
  te_sum <- sum(te)      # sum LDs per individual
  lf_hat <- mean(Effects$years) # average Life years Gained
  
  # store the results from the simulation in a list
  results <- list(m_M = m_M, 
                  m_C = m_C, 
                  m_E = m_E, 
                  tc = tc , 
                  te = te, 
                  tc_hat = tc_hat, 
                  te_hat = te_hat,
                  tc_sum = tc_sum, 
                  te_sum = te_sum,
                  Ly_s   = lf_hat)   
  
  return(results)  # return the results
  
} # end of the MicroSim function  

outcomes <- MicroSim(n_i, df_X, df_pop_ch, Life_expectancy, Trt = FALSE)

results  <- data.frame("Total Cost" = outcomes$tc_hat, 
                       "Total Lys" = outcomes$Ly_s)



