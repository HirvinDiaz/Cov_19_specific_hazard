### Markov Model alive - Dead to calculate QALE (Quality Adjusted Life Expectancy)
### For each patient in the cohort. Model based in the work 
# - Jalal H, Pechlivanoglou P, Krijkamp E, Alarid-Escudero F, Enns E, Hunink MG. 
# An Overview of R in Health Decision Sciences. Med Decis Making. 2017; 37(3): 735-746. 
# https://journals.sagepub.com/doi/abs/10.1177/0272989X16686559
# - Krijkamp EM, Alarid-Escudero F, Enns EA, Jalal HJ, Hunink MGM, Pechlivanoglou P. 
# Microsimulation modeling for health decision sciences using R: A tutorial. 
# Med Decis Making. 2018;38(3):400–22. 
# https://journals.sagepub.com/doi/abs/10.1177/0272989X18754513
### Adapted by Hirvin Diaz 
### Mail: hirvinazael@gmail.com

#### 01 Load packages ####
library(data.table)
library(reshape2)
library(readr)
library(readxl)
library(dplyr)

# Load functions for microsimulation
source("R/Functions.R")

#### 02 Load databases ####

# load("data/cov_04_04_comor.Rdata") 

# <- Cov %>% 
#   filter(type == 2)

# save(Cov, file = "data/Cov_Alive_Dead_Markov")

load("data/Cov_Alive_Dead_Markov")

#### 03 Input Model Parameters ####
set.seed(1)   # set the seed 

Cov$sex[Cov$sex == 1] <- "Female"
Cov$sex[Cov$sex == 2] <- "Male"

#### 04 Data wrangling ####

Cov_H_45 <- Cov %>% 
  filter(age == 45 & sex == "Male") %>% 
  select(c(2,5))

#### Age 45 - 50 Healthy ####
# Model structure
n_i       <- length(Cov_H_45$sex)      # number of individuals
n_t       <- 55                           # number of cycles
v_names_states  <- c("Alive", "Dead")     # state names
n_s  <- length(v_names_states)            # number of states
d_e <- 0.05                               # discount rate for QALYs at 5%
d_c <- 0.05                               # discount rate for costs at 5%

# calculate discount weights for costs for each cycle based on discount rate d_c
v_dwc <- 1 / (1 + d_e) ^ (0:n_t) 
# calculate discount weights for effectiveness for each cycle based on discount
# rate d_e
v_dwe <- 1 / (1 + d_c) ^ (0:n_t)

# Probabilities of death from healthy state is age-dependent
# load probability of death by age from Human Mortality Database

df_hmd_MX_2020 <- read.csv("data-raw/df_mortrate_state_age_sex.csv") %>% 
  filter(year == 2020 & state == "National") 

# create a "long" data frame that we can index by age and sex
df_p_AD <- df_hmd_MX_2020 %>% 
  select("age", "sex", "mort_rate") %>% 
  rename(r_AD = mort_rate) %>% 
  filter(age <= 101)

# calculate probability from death rate
df_p_AD$p_AD <- 1-exp(-df_p_AD$r_AD)

dtb_p_AD <- data.table(df_p_AD) # convert data frame to a data table for efficiency
setkey(dtb_p_AD, age, sex)      # set the data table to be indexed by Age and Sex

# Costs and utilities inputs 
c_A  <- 11483.3                 # cost of remaining one cycle Alive
c_D  <- 0                       # cost of remaining one cycle Dead
u_A  <- 0.827                   # utility when Alive 
u_D  <- 0                       # utility when Sick

#### 04 Sample individual level characteristics ####
#### 04.1 Static characteristics
# Store these characteristics in a single data frame for convenient referencing

# Initial Age

# Store these static individual characteristics in one data frame
df_X <- data.frame(Ind = 1:n_i, age = Cov_H$age, sex = Cov_H$sex)
head(df_X)

#### 04.2 Dynamic characteristics 
#### 05 Define Simulation Functions ####

# ProbsDeath function:
# looks up probability of death for given cycle
Probs <- function(v_M_t, df_X, t) { 
  # M_t  :  current health states
  # v_ai :  initial ages
  # v_sex:  sexes
  # t    :  cycle
  # create matrix of state transition probabilities
  m_p_t           <- matrix(0, nrow = n_s, ncol = n_i) 
  # give the state names to the rows
  rownames(m_p_t) <-  v_names_states
  
  # lookup probability baseline probability of dying based on individual characteristics
  p_die_base <- dtb_p_AD[.(df_X$age + t, df_X$sex), p_AD]  # <- how you index a data table
  p_die      <- p_die_base[v_M_t == "Alive"]
  
  # transition probabilities when Alive
  m_p_t[, v_M_t == "Alive"] <- rbind(1 - (p_die), p_die) 
  # transition probabilities when Dead
  m_p_t[, v_M_t == "Dead"] <- rbind(0, 1)
  
  return(t(m_p_t))
}

# Costs function: 
# calculates the cost accrued by an individual this cycle
Costs <- function (M_t) {
  # M_t: current health state
  c_t <- c()
  c_t[M_t == "Alive"] <- c_A   # costs accrued by being Healthy this cycle
  c_t[M_t == "Dead"] <- c_D   # costs at Dead state
  
  return(c_t)  # return costs accrued this cycle
}

# Effs function:
# calculates QALYs accrued by an individual for this cycle
Effs <- function (M_t) {
  # M_t: current health state
  q_t <- c() 
  q_t[M_t == "Alive"] <- u_A  # QALYs accrued by being Healthy this cycle
  q_t[M_t == "Dead"] <- u_D  # QALYs at Dead state
  return(q_t)  # return the QALYs accrued this cycle
}


#### 06 Run Microsimulation ####

v_M_init  <- rep("Alive", n_i)       # everyone begins in the healthy state

MicroSim <- function(n_i, df_X, seed = 1) { #t <- 1
  
  set.seed(seed) # set the seed
  
  m_M <- m_C <- m_E <-  matrix(NA, nrow = n_i, ncol = n_t + 1, 
                               dimnames = list(paste("ind"  , 1:n_i, sep = " "), 
                                               paste("cycle", 0:n_t, sep = " ")))  
  
  m_M[, 1] <- v_M_init          # initial health state
  m_C[, 1] <- Costs(m_M[, 1])   # costs accrued during cycle 0
  m_E[, 1] <- Effs(m_M[, 1])    # QALYs accrued during cycle 0
  
  # open a loop for time running cycles 1 to n_t 
  for (t in 1:n_t) { # t <- 1
    # calculate the transition probabilities for the cycle based on health state t
    m_P <- Probs(v_M_t = m_M[, t], df_X, t = t) 
    # sample the current health state and store that state in matrix m_M
    m_M[, t + 1]  <- samplev(m_P, 1)    
    # calculate costs per individual during cycle t + 1
    m_C[, t + 1]  <- Costs(m_M[, t + 1])  
    # calculate QALYs per individual during cycle t + 1
    m_E[, t + 1]  <- Effs(m_M[, t + 1])  
    
    # Display simulation progress
    if(t/(n_t/10) == round(t/(n_t/10), 0)) { # display progress every 10%
      cat('\r', paste(t/n_t * 100, "% done", sep = " "))
    }
    
  } # close the loop for the time points 
  
  # calculate  
  tc <- m_C %*% v_dwc    # total (discounted) cost per individual
  te <- m_E %*% v_dwe    # total (discounted) QALYs per individual 
  tc_hat <- mean(tc)     # average (discounted) Cost 
  te_hat <- mean(te)     # average (discounted) QALYs
  
  # store the results from the simulation in a list
  results <- list(m_M = m_M, 
                  m_C = m_C, 
                  m_E = m_E, 
                  tc = tc, 
                  te = te, 
                  tc_hat = tc_hat, 
                  te_hat = te_hat)   
  
  return(results)  # return the results
  
} # end of the MicroSim function  


################################################################################

outcomes <- MicroSim(n_i, df_X, seed = 1)

results  <- data.frame("Total Cost" = outcomes$tc_hat, 
                       "Total QALYs" = outcomes$te_hat)

###############################################################################

df_cost   <- as.data.frame(matrix(0, 
                                nrow = 5,
                                ncol = 1))

df_effect <- as.data.frame(matrix(0, 
                                nrow = 5,
                                ncol = 1))


## Model ####
for (i in 45:49) {
  #### Data wrangling ####
  Cov_H <- Cov %>% 
    filter(age == i & sex == "Male") %>% 
    select(c(2,5))
  
  # Model structure
  n_i       <- length(Cov_H$sex)      # number of individuals
  n_t       <- 101 - i                 # number of cycles
  
  # calculate discount weights for costs for each cycle based on discount rate d_c
  v_dwc <- 1 / (1 + d_e) ^ (0:n_t) 
  # calculate discount weights for effectiveness for each cycle based on discount
  v_dwe <- 1 / (1 + d_c) ^ (0:n_t)
  
  # Costs and utilities inputs 
  c_A  <- 11483.3                 # cost of remaining one cycle Alive
  c_D  <- 0                       # cost of remaining one cycle Dead
  u_A  <- 0.871	                  # utility when Alive 
  u_D  <- 0                       # utility when Sick
  
  Probs <- function(v_M_t, df_X, t) { 
    # M_t  :  current health states
    # v_ai :  initial ages
    # v_sex:  sexes
    # t    :  cycle
    # create matrix of state transition probabilities
    m_p_t           <- matrix(0, nrow = n_s, ncol = n_i) 
    # give the state names to the rows
    rownames(m_p_t) <-  v_names_states
    
    # lookup probability baseline probability of dying based on individual characteristics
    p_die_base <- dtb_p_AD[.(df_X$age + t, df_X$sex), p_AD]  # <- how you index a data table
    p_die      <- p_die_base[v_M_t == "Alive"]
    
    # transition probabilities when Alive
    m_p_t[, v_M_t == "Alive"] <- rbind(1 - (p_die), p_die) 
    # transition probabilities when Dead
    m_p_t[, v_M_t == "Dead"] <- rbind(0, 1)
    
    return(t(m_p_t))
  }
  
  
  # Costs function: 
  # calculates the cost accrued by an individual this cycle
  Costs <- function (M_t) {
    # M_t: current health state
    c_t <- c()
    c_t[M_t == "Alive"] <- c_A   # costs accrued by being Healthy this cycle
    c_t[M_t == "Dead"] <- c_D   # costs at Dead state
    
    return(c_t)  # return costs accrued this cycle
  }
  
  # Effs function:
  # calculates QALYs accrued by an individual for this cycle
  Effs <- function (M_t) {
    # M_t: current health state
    q_t <- c() 
    q_t[M_t == "Alive"] <- u_A  # QALYs accrued by being Healthy this cycle
    q_t[M_t == "Dead"] <- u_D  # QALYs at Dead state
    return(q_t)  # return the QALYs accrued this cycle
  }
  
  
  df_X <- data.frame(Ind = 1:n_i, age = Cov_H$age, sex = Cov_H$sex)
  
  v_M_init  <- rep("Alive", n_i) 
  
  MicroSim <- function(n_i, df_X, seed = 1) { #t <- 1
    
    set.seed(seed) # set the seed
    
    m_M <- m_C <- m_E <-  matrix(NA, nrow = n_i, ncol = n_t + 1, 
                                 dimnames = list(paste("ind"  , 1:n_i, sep = " "), 
                                                 paste("cycle", 0:n_t, sep = " ")))  
    
    m_M[, 1] <- v_M_init          # initial health state
    m_C[, 1] <- Costs(m_M[, 1])   # costs accrued during cycle 0
    m_E[, 1] <- Effs(m_M[, 1])    # QALYs accrued during cycle 0
    
    # open a loop for time running cycles 1 to n_t 
    for (t in 1:n_t) { # t <- 1
      # calculate the transition probabilities for the cycle based on health state t
      m_P <- Probs(v_M_t = m_M[, t], df_X, t = t) 
      # sample the current health state and store that state in matrix m_M
      m_M[, t + 1]  <- samplev(m_P, 1)    
      # calculate costs per individual during cycle t + 1
      m_C[, t + 1]  <- Costs(m_M[, t + 1])  
      # calculate QALYs per individual during cycle t + 1
      m_E[, t + 1]  <- Effs(m_M[, t + 1])  
      
      # Display simulation progress
      if(t/(n_t/10) == round(t/(n_t/10), 0)) { # display progress every 10%
        cat('\r', paste(t/n_t * 100, "% done", sep = " "))
      }
      
    } # close the loop for the time points 
    
    # calculate  
    tc <- m_C %*% v_dwc    # total (discounted) cost per individual
    te <- m_E %*% v_dwe    # total (discounted) QALYs per individual 
    tc_hat <- mean(tc)     # average (discounted) Cost 
    te_hat <- mean(te)     # average (discounted) QALYs
    
    # store the results from the simulation in a list
    results <- list(m_M = m_M, 
                    m_C = m_C, 
                    m_E = m_E, 
                    tc = tc, 
                    te = te, 
                    tc_hat = tc_hat, 
                    te_hat = te_hat)   
    
    return(results)  # return the results
    
  } # end of the MicroSim function  
  outcomes <- MicroSim(n_i, df_X, seed = 1)
  
  df_cost[i-44, 1] <- outcomes$tc_hat
  df_effect[i-44, 1] <- outcomes$te_hat
  
}

df_cost
df_effect



# Modelo Markov, asignando la edad inicial, 
# no hay error en la estimacion de esperanza de vida promedio
# Tabla de age-sepcific hazard y age specific-utilities.