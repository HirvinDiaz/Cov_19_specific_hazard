#### 01 Load packages ####
library(data.table)
library(reshape2)
library(readr)
library(readxl)
library(dplyr)

# Poner un mean population de intubates y not intubated
# The not intubated cohort with an averga age lives bla

# New Markov QALEs# 

source("R/Markov_3state_time_functions.R")

u_45_49_m <- 0.887 
u_50_59_m <- 0.861
u_60_69_m <- 0.840
u_70_79_m <- 0.802
u_80_100_m <- 0.782


u_45_49_f <- 0.863
u_50_59_f <- 0.837
u_60_69_f <- 0.811
u_70_79_f <- 0.771
u_80_100_f <- st0.724


df_hmd_MX_2020 <- read.csv("data-raw/df_mortrate_state_age_sex.csv") %>% 
  filter(year == 2020 & state == "National") 

c_A_m_45_59     <- 10912
c_A_f_45_59     <- 13519
c_A_m_60_mas     <- 24383

# Mort rate, costos para cada año de edad, utilidades y los gastos.

# Weights by age

cero <- sum(df_hmd_MX_2020$population[df_hmd_MX_2020$age == 0])
uno_y_cinco <- sum(df_hmd_MX_2020$population[df_hmd_MX_2020$age >= 1 & df_hmd_MX_2020$age <= 5])
cinco_y_nueve <- sum(df_hmd_MX_2020$population[df_hmd_MX_2020$age >= 5 & df_hmd_MX_2020$age <= 9])
diez_y_diecinueve <- sum(df_hmd_MX_2020$population[df_hmd_MX_2020$age >= 10 & df_hmd_MX_2020$age <= 19])
veinte_y_cinc_nueve_male <- sum(df_hmd_MX_2020$population[df_hmd_MX_2020$age >= 20 & df_hmd_MX_2020$age <= 59 & df_hmd_MX_2020$sex == "Male"])
veinte_y_cinc_nueve_female <- sum(df_hmd_MX_2020$population[df_hmd_MX_2020$age >= 20 & df_hmd_MX_2020$age <= 59 & df_hmd_MX_2020$sex == "Female"])
mas_sesenta <- sum(df_hmd_MX_2020$population[df_hmd_MX_2020$age >= 60])

# Strategy names
v_names_str <- c("Base Case")

# Number of strategies
n_str <- length(v_names_str)

# Markov model parameters
v_n  <- c("Alive","Dead")    # State names
n_states  <- length(v_n)     # Number of states
age     <- 45                # age at baseline
n_t     <- max_age - age     # time horizon, number of cycles
n_s     <- length(v_n)       # number of health states 

# Transition probabilities (per cycle)
haz_mort <- df_hmd_MX_2020 %>% 
  filter(sex == "Male" & age >= 45 & age <= 99) 
haz_mort <- as.vector(haz_mort[["mort_rate"]])
#vector (loop age and sex y que adentro me arroje la probabilidades)        	 
v_prob   <- 1 - exp(-haz_mort) 


# Cost and utility inputs 

u_A     <- u_45_49_m  # Generar un vector de utilidades a la edad 100
u_D     <- 0          # Utility when dead

# Hacer una tabla de utilidades y voy a tener edad y utilidad (por cada edad) (50 - 54) por edad

edad <- rep(seq(45,100, by = 1),2)

sex <- c(rep("Male", 56), rep("Female", 56))

df_utility <- as.data.frame(cbind(edad,sex)) 

df_utility <- df_utility %>% 
  mutate(Utils = ifelse(edad >= 45 & edad <= 49 & sex == "Male", u_45_49_m, 
                        ifelse(edad >= 45 & edad <= 49 & sex == "Female", u_45_49_f, 
                               ifelse(edad >= 50 & edad <= 59 & sex == "Male", u_50_59_m, 
                                      ifelse(edad >= 50 & edad <= 59 & sex == "Female", u_50_59_f,
                                             ifelse(edad >= 60 & edad <= 69 & sex == "Male", u_60_69_m,
                                                    ifelse(edad >= 60 & edad <= 69 & sex == "Female", u_60_69_f,
                                                           ifelse(edad >= 70 & edad <= 79 & sex == "Male", u_70_79_m,
                                                                  ifelse(edad >= 70 & edad <= 79 & sex == "Female", u_70_79_f,
                                                                         ifelse(edad >= 80 & sex == "Male", u_80_100_m,u_80_100_f))))))))))

df_utility <- df_utility %>% 
  mutate(Costs = ifelse(edad >= 45 & edad <= 59 & sex == "Male", c_A_m_45_59,
                        ifelse(edad >= 45 & edad <= 59 & sex == "Female", c_A_f_45_59, c_A_m_60_mas)))

df_utility$edad <- as.numeric(df_utility$edad)

class(df_utility$edad)
# Discounting factor
d_r     <- 0.05 # equal discount of costs and QALYs by 3%
# calculate discount weights for costs for each cycle based on discount rate d_c
v_dwc   <- 1 / (1 + d_r) ^ (0:n_t) 
# Vector de costos por edad
# calculate discount weights for effectiveness for each cycle based on discount rate d_e
v_dwe   <- 1 / (1 + d_r) ^ (0:n_t) 


# Para descontar los costos apkico el vector de costos, 
# La longitud del vector de descuentos cambia dependiendo la edad, el vector se lo aplico a la traza, voy a generar una traza
# 56 filas, la traza (agarro la columna) la multiplico por costos o utilidades, tenemos el vector de costos
# y descuentos tienen el mismo tamaño, los multiplico y lo sumo. 

# 112 modelos markov (56 y 56) Me dan 112 valores de QALEs y de costos.

# SUGERENCIA: Calcular todo con descuento y sin descuento. (Lys y Lys descontados (Lys = columna de la traza de los vivos
# Lys por el vector de descuento.)) 

# Mean age de la muestra, tomar los que sobrevivieron y ver que me da

### Meta_modelaje 

### Linear regression metamodeling for sensitivity analysis. #Dampack

# create the markov trace matrix M capturing the proportion of the cohort 
# in each state at each cycle
#******************************************************************************
#### 04.2 Transition probability array ####
#******************************************************************************

m_P_diag <- matrix(0,
                   nrow     = n_states,
                   ncol     = n_states,
                   dimnames = list(v_n, v_n))
m_P_diag["Alive", "Dead" ]     = ""
m_P_diag["Dead" , "Dead" ]     = ""
m_P_diag["Alive", "Alive" ]     = ""
m_P_diag["Dead" , "Alive" ]     = ""


#******************************************************************************
#### 04.1 Cohort trace ####
#******************************************************************************

# Create the cohort trace
m_M <- matrix(NA,
              # Create Markov trace (n_t + 1 because R doesn't understand Cycle 0)
              nrow = n_t + 1,
              ncol = n_states,
              dimnames = list(0:n_t, v_n))

# Initialize first cycle of Markov trace
m_M[1, ] <- c(1, 0)

#******************************************************************************
#### 04.2 Transition probability array ####
#******************************************************************************

# # Create the transition probability matrix
# m_P  <- matrix(0,
#                nrow     = n_states,
#                ncol     = n_states,
#                dimnames = list(v_n, v_n)) # Name the columns and rows of the transition
#
# # Probability matrix
# m_P

# Create 3-D array
a_P <- array(0,
             dim = c(n_states, n_states, n_t),
             dimnames = list(v_n, v_n, 0:(n_t - 1))) # Name the dimensions of the array

# Fill in the transition probability array:
a_P["Alive","Alive",] <- 1 - v_prob
a_P["Alive", "Dead",]    <- v_prob

# From Dead
a_P["Dead", "Dead", ] <- 1

# Transition probability array
#******************************************************************************
#### 04.3 Check if transition array and probabilities are valid ####
#******************************************************************************

# # Check rows add up to 1
# rowSums(m_P)

# Check if transition probability array is valid (i.e., elements cannot < 0 or > 1)
check_transition_probability(a_P,
                             verbose = TRUE) # Prints a message if there's an error

# Check if transition probability array sum to 1 (i.e., each row should sum to 1)
check_sum_of_transition_array(a_P,
                              n_states = n_states,
                              n_t      = n_t,
                              verbose  = TRUE) # Prints a message if there's an error
for (t in 1:n_t){
  # Estimate the state vector for the next cycle (t + 1)
  m_M[t + 1, ] <- m_M[t, ] %*% a_P[, , t]
}

#### Lys ####
v_os <- 1 - m_M[, "Dead"]
v_le <- sum(v_os)

#### Lys discounted ####
v_ly <- t(v_os) %*% v_dwe

#******************************************************************************
#### 07.1 Mean Costs and QALYs ####
#******************************************************************************

# Per cycle
# Calculate expected costs by multiplying m_M with the cost vector for the different
# Health states

Costs <- df_utility %>% 
  filter(sex == "Male")

v_costs <- as.vector(Costs[["Costs"]])

v_tc <- m_M[,1] * v_costs

# Calculate expected QALYs by multiplying m_M with the utilities for the different
# Health states
Utils <- df_utility %>% 
  filter(sex == "Male")

v_Utils <- cbind(as.vector(Costs[["Utils"]]))

v_tu <- m_M[,1] * v_Utils

#******************************************************************************
#### 07.2 Discounted Mean Costs and QALYs ####
#******************************************************************************

# Discount costs by multiplying the cost vector with discount weights (v_dw)
n_tc_d <-  t(v_tc) %*% v_dwc

# Discount QALYS by multiplying the QALYs vector with discount weights (v_dw)
n_te_d <-  t(v_tu) %*% v_dwe

#******************************************************************************
#### 07.3 Results ####
#******************************************************************************

results <- data.frame( "Total Discounted Cost" = n_tc_d, 
                       "Total Discounted QALYs" = n_te_d, 
                       check.names = F)



df_costs_test <- df_utility %>% 
  filter(sex == "Female" & age >= 52)

costs_test <- as.vector(df_costs_test[["Costs"]])

#### Markov Loop ####
# Max Age
#max_age <- 50               

# Discounting factor
d_r     <- 0.05 

#Sex
Sex <- "Female"

df_dis_cost   <- as.data.frame(matrix(0, 
                                  nrow = 56,
                                  ncol = 1))
colnames(df_dis_cost) <- "dis_costs"

df_cost   <- as.data.frame(matrix(0, 
                                      nrow = 56,
                                      ncol = 1))
colnames(df_cost) <- "costs"

df_qalys <- as.data.frame(matrix(0, 
                                  nrow = 56,
                                  ncol = 1))
colnames(df_qalys) <- "QALYs"

df_lys <- as.data.frame(matrix(0, 
                                  nrow = 56,
                                  ncol = 1))
colnames(df_lys) <- "lys"

df_dis_lys <- as.data.frame(matrix(0, 
                               nrow = 56,
                               ncol = 1))
colnames(df_dis_lys) <- "dis_lys"

for (t in 45:100) { #t <- 45
  age     <- t                # age at baseline
  n_t     <- 101 - age        # time horizon, number of cycles
  n_s     <- length(v_n)       # number of health states 
  # calculate discount weights for costs for each cycle based on discount rate d_c
  v_dwc   <- 1 / (1 + d_r) ^ (0:n_t) 
  v_dwe   <- 1 / (1 + d_r) ^ (0:n_t) 
  
  # Transition probabilities (per cycle)
  haz_mort <- df_hmd_MX_2020 %>% 
    filter(sex == Sex & age >= t & age <= 100) 
  
  haz_mort <- as.vector(haz_mort[["mort_rate"]])
  #vector (loop age and sex y que adentro me arroje la probabilidades)        	 
  v_prob   <- 1 - exp(-haz_mort) 
  
  m_P_diag <- matrix(0,
                     nrow     = n_states,
                     ncol     = n_states,
                     dimnames = list(v_n, v_n))
  m_P_diag["Alive", "Dead" ]     = ""
  m_P_diag["Dead" , "Dead" ]     = ""
  m_P_diag["Alive", "Alive" ]     = ""
  m_P_diag["Dead" , "Alive" ]     = ""
  
  # Create the cohort trace
  m_M <- matrix(NA,
                # Create Markov trace (n_t + 1 because R doesn't understand Cycle 0)
                nrow = n_t + 1,
                ncol = n_states,
                dimnames = list(0:n_t, v_n))
  
  # Initialize first cycle of Markov trace
  m_M[1, ] <- c(1, 0)
  
  # Create 3-D array
  a_P <- array(0,
               dim = c(n_states, n_states, n_t),
               dimnames = list(v_n, v_n, 0:(n_t - 1))) # Name the dimensions of the array
  
  # Fill in the transition probability array:
  a_P["Alive","Alive",] <- 1 - v_prob
  a_P["Alive", "Dead",]    <- v_prob
  
  # From Dead
  a_P["Dead", "Dead", ] <- 1

  for (i in 1:n_t){
    # Estimate the state vector for the next cycle (t + 1)
    m_M[i + 1, ] <- m_M[i, ] %*% a_P[, , i]
  }
  
  #### Lys ####
  v_os <- 1 - m_M[, "Dead"]
  v_le <- sum(v_os)
  
  #### Lys discounted ####
  v_ly <- t(v_os) %*% v_dwe
  
  df_costs_ut <- df_utility  %>% 
    filter(sex == Sex & edad >= t)
  
  v_costs <- as.vector(c(0,df_costs_ut[["Costs"]]))
  v_tc <- m_M[,1] * v_costs
  
  # Calculate expected QALYs by multiplying m_M with the utilities for the different
  # Health states
  
  v_Utils <- as.vector(c(0,df_costs_ut[["Utils"]]))
  
  v_tu <- m_M[,1] * v_Utils
  
  #******************************************************************************
  #### 07.2 Discounted Mean Costs and QALYs ####
  #******************************************************************************
  
  # Discount costs by multiplying the cost vector with discount weights (v_dw)
  n_tc_d <-  t(v_tc) %*% v_dwc
  
  # Discount QALYS by multiplying the QALYs vector with discount weights (v_dw)
  n_te_d <-  t(v_tu) %*% v_dwe
  
  #******************************************************************************
  #### 07.3 Results ####
  #******************************************************************************
  
  results <- data.frame( "Dis_cost" = n_tc_d, 
                         "QALYs" = n_te_d,
                         "cost" = sum(v_tc),
                         "LYs" = v_le,
                         "dis_LYs" = v_ly,
                         check.names = F)
  
  df_dis_cost[t - 44, 1] <- results$Dis_cost
  df_cost[t - 44, 1] <- results$cost
  df_qalys[t- 44, 1] <- results$QALYs
  df_lys[t- 44, 1] <- results$LYs
  df_dis_lys[t- 44, 1] <- results$dis_LYs
}


df_QALEs_females <- bind_cols(df_dis_cost,df_qalys, df_cost, df_dis_lys, df_lys)

df_QALEs_females <- df_QALEs_females %>% 
  mutate(sex = "Female", age = seq(45,100, by = 1)) %>% 
  select(7,6,1:5)

df_QALEs <- bind_rows(df_QALEs_males, df_QALEs_females)

save(df_QALEs, file = "data/df_QALEs_markov.Rdata")
