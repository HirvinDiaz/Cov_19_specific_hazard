#### Creation of data tables ####

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

#### Load data ####
#load("data/hazard_hosp.Rdata")
#load("data/hazard_intubated.Rdata")
load("data/hazard_hos_res.Rdata")

#### Load ratetable ####
load("data/rate_mx.Rdata")

#### Hospitalized population ####
#### transform original database to data table with probabilities ####
#data_table_hosp <- reshape(df_hazards_hosp, 
#        idvar = c("time","group", "sex"), 
#        timevar = "Type", 
#        direction = "wide")
#
#data_table_hosp <- data_table_hosp %>% 
#  mutate(cohort = "hosp") %>% 
#  mutate(p_dCoV = 1 - exp(-Hazard.Covid_19)) %>% 
#  mutate(p_dPop = 1 - exp(- Hazard.Population)) %>% 
#  arrange(group, time) %>% 
#  select(c("time", "group", "sex", "cohort", "p_dCoV", "p_dPop")) %>% 
#  rename(day = time)

#### Respiratory auxiliar population ####
#### transform original database to data table with probabilities ####
#data_table_resp <- reshape(df_hazards_intubated, 
#                           idvar = c("time","group", "sex"), 
#                           timevar = "Type", 
#                           direction = "wide")
#
#data_table_resp <- data_table_resp %>% 
#  mutate(cohort = "resp") %>% 
#  mutate(p_dCoV = 1 - exp(-Hazard.Covid_19)) %>% 
#  mutate(p_dPop = 1 - exp(- Hazard.Population)) %>% 
#  arrange(group, time) %>% 
#  select(c("time", "group", "sex", "cohort", "p_dCoV", "p_dPop")) %>% 
#  rename(day = time)

#d_p_HD <- rbind(data_table_hosp)
                #data_table_resp)

# Hospitalized and 

#### Respiratory auxiliar population ####
#### transform original database to data table with probabilities ####
data_table_hos <- reshape(df_hazards_hos_res, 
                           idvar = c("time","group", "sex"), 
                           timevar = "Type", 
                           direction = "wide")

data_table_res <- reshape(df_hazards_hos_res, 
                          idvar = c("time","group", "sex"), 
                          timevar = "Type", 
                          direction = "wide")


data_table_resp <- data_table_res %>% 
  mutate(pop = "int") %>% 
  mutate(p_dCoV = 1 - exp(-Hazard.Covid_19)) %>% 
  mutate(p_dPop = 1 - exp(- Hazard.Population)) %>% 
  arrange(group, time) %>% 
  select(c("time", "group", "sex", "pop", "p_dCoV", "p_dPop")) %>% 
  rename(day = time)

data_table_hosp <- data_table_hos %>% 
  mutate(pop = "hos") %>% 
  mutate(p_dCoV = 1 - exp(-Hazard.Covid_19)) %>% 
  mutate(p_dPop = 1 - exp(- Hazard.Population)) %>% 
  arrange(group, time) %>% 
  select(c("time", "group", "sex", "pop", "p_dCoV", "p_dPop")) %>% 
  rename(day = time)

d_p_HD <- rbind(data_table_hosp, data_table_resp)

#### Cohort to be simulated ####
load("data/cov_24_01.Rdata") 

##### Covid_rr for intubated and hospitalized individuals higher than 44 years 
# in national DB #####
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

#length(Cohort$intubated[Cohort$intubated == 1])
#length(Cov$intubated[Cov$intubated == 1 & Cov$type != 2])


Cohort <- Cohort %>% 
  mutate(pop = ifelse(intubated == 1, "int", "hos")) %>%
  select(c("sex","age_range","pop")) %>%
  rename(group = age_range) %>% 
  mutate(day = 0)

# length(Cohort$sex[Cohort$pop == "int"])

#### Microsimulation without treatment ####
# Parameters set
# number of simulated individuals, in this case 488,866
n_i              <- length(Cohort$sex) 
n_t              <- 60 # time horizon, 60 days
# model states: Positive Case - Cov19+, Death by Covid-19 - CoV19_Dead, 
# Death by other causes - O_Causes_Dead
v_names_states   <- c("Cov19+", "CoV19_Dead", "O_Causes_Dead") 
n_states         <- length(v_names_states) # the number of states
d_c              <- d_e <- (0.017 + 0.016)/2 # Daily discount rates for costs and utilities
# calculate discount weights for costs for each cycle based on discount rate d_c
v_dwc <- 1 / (1 + d_c) ^ (0:n_t) 
# calculate discount weights for effectiveness for each cycle based on discount 
# rate d_e
v_dwe <- 1 / (1 + d_e) ^ (0:n_t)

## Costs and utilities inputs (in MX pesos)
c_amb <- 14500 # Average cost by ambulatory patient
# Average cost for patients that require hospitalized care
c_hosp <- ((35000 + 50000 + 70000+ 80000)/4)
# compute proportion of ambulatory and hospitalized patients  
# Proportion of hospitalized patients
# cost of remaining one cycle sick with COVID-19
c_sCov     <- c_hosp + c_amb 
c_dCov     <- 0       # cost of remaining one cycle Dead
c_dPop     <- 0       # cost of remaining one cycle Dead
c_Trt      <- 0
# Mean QALD (Quality Adjusted Life Days) loss.
m_QALD     <-  2.5
u_sCov     <- (100 - m_QALD)/100    # utility when Sick 
u_dCov     <- 0       # utility when Dead
u_dPop     <- 0       # utility when Dead

# Convert data frame to a data table for efficiency
dt_p_CoV <- data.table(d_p_HD)

# set the data table to be indexed by age, day and Sex
setkey(dt_p_CoV, group, day, sex, pop)

# Create data frame of population from cohort. All begin in day 0
df_X <- Cohort

source("R/Functions.R")

Probs <- function(v_M_t, df_X, t) { # t <- 1
  # Arguments:
  # v_M_t: health state occupied at cycle t (character variable)
  # df_X: data frame with individual characteristics data 
  # v_Ts: vector with the duration of being sick
  # t: cycle
  # Returns: 
  # transition probabilities for that cycle
  # create matrix of state transition probabilities
  m_p_t           <- matrix(0, nrow = n_states, ncol = n_i)  
  # give the state names to the rows
  rownames(m_p_t) <-  v_names_states
  
  # Lookup baseline probability of dying from Covid-19 or other causes based on  
  # individual characteristics of day, sex and age 
  p_die_CoV_all <- dt_p_CoV[.(df_X$group, df_X$day + t, df_X$sex, df_X$pop), p_dCoV] 
  p_die_Pop_all <- dt_p_CoV[.(df_X$group, df_X$day + t, df_X$sex, df_X$pop), p_dPop] 
  p_die_CoV     <- p_die_CoV_all[v_M_t == "Cov19+"]  
  p_die_Pop     <- p_die_Pop_all[v_M_t == "Cov19+"] 
  
  # update m_p_t with the appropriate probabilities   
  # transition probabilities when healthy
  m_p_t[, v_M_t == "Cov19+"] <- rbind(1 - (p_die_CoV + p_die_Pop), p_die_CoV, p_die_Pop)    
  # transition probabilities when sick 
  m_p_t[, v_M_t == "CoV19_Dead"] <- rbind(0, 1 ,0)
  # transition probabilities when sicker 
  m_p_t[, v_M_t == "O_Causes_Dead"] <- rbind(0, 0, 1)
  
  
  # t(m_p_t[, 1:10]) # Show the probabilities for the first 10 individuals 
  
  return(t(m_p_t))
}     

Costs <- function (v_M_t, Trt = FALSE) {
  # v_M_t: current health state
  c_t <- c()
  c_t[v_M_t == "Cov19+"]  <- c_sCov + (c_Trt * Trt)    # costs accrued by being healthy this cycle
  c_t[v_M_t == "CoV19_Dead"] <- c_dCov  # costs accrued by being sick this cycle
  c_t[v_M_t == "O_Causes_Dead"] <- c_dCov  # costs accrued by being sicker this cycle
  return(c_t)  # return costs accrued this cycle
}

Effs <- function (v_M_t, Trt = FALSE) {
  # v_M_t: current health state
  q_t <- c() 
  q_t[v_M_t == "Cov19+"] <- u_sCov   
  # QALYs accrued by being healthy this cycle
  q_t[v_M_t == "CoV19_Dead"] <- u_dCov     # QALYs accrued by being healthy this cycle
  q_t[v_M_t == "O_Causes_Dead"]  <- u_dCov     # QALYs accrued by being sick this cycle
  return(q_t)  # return the QALYs accrued this cycle
}

#### 04.2 Dynamic characteristics 
# These are just starting conditions - they will change with the simulation
v_M_init  <- rep("Cov19+", n_i)       # everyone begins in the healthy state

MicroSim <- function(n_i, df_X) { #t <- 1
  
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
  tc_hat <- mean(tc)     # average (discounted) cost 
  te_hat <- mean(te)     # average (discounted) QALYs
  tc_sum <- sum(tc)      # sum (discounted) cost
  te_sum <- sum(te)      # sum (discounted) QALYs
  
  # store the results from the simulation in a list
  results <- list(m_M = m_M, 
                  m_C = m_C, 
                  m_E = m_E, 
                  tc = tc , 
                  te = te, 
                  tc_hat = tc_hat, 
                  te_hat = te_hat,
                  tc_sum = tc_sum, 
                  te_sum = te_sum )   
  
  return(results)  # return the results
  
} # end of the MicroSim function  

outcomes <- MicroSim(n_i, df_X)

results  <- data.frame("Total Cost" = outcomes$tc_hat, 
                       "Total QALYs" = outcomes$te_hat)

# Create dataframe to transform: alive =1 and death = 0
state_m_s <- as.data.frame(outcomes$m_M) 
state_m_s[state_m_s == "Cov19+"] <- 1
state_m_s[state_m_s == "CoV19_Dead"] <- 0 
state_m_s[state_m_s == "O_Causes_Dead"] <- 0 
state_m_s <- sapply(state_m_s, as.numeric)

Cycle_sum_NT <- 
  as.matrix(colSums(state_m_s)/n_i) # Compute survival probability for each cycle
print(Cycle_sum_NT)

length(Cycle_sum_NT)

#### Microsimulate 1000 times ####

#d <- as.data.frame(matrix(nrow = 61, 
#                          ncol = 100, 
#                          data = NA))
#G <- seq(1, 61, by = 1)
#
#d[,2]
#
#for (i in 1:100) {
# G <- seq(1, 61, by = 1)
# d[,i] <- G
#}

KM_micro_df <- as.data.frame(matrix(nrow = 61, 
                                    ncol = 1000, 
                                    data = NA))

for (i in 1:1000) {
  outcomes <- MicroSim(n_i, df_X)
  
  results  <- data.frame("Total Cost" = outcomes$tc_hat, 
                         "Total QALYs" = outcomes$te_hat)
  
  # Create dataframe to transform: alive =1 and death = 0
  state_m_s <- as.data.frame(outcomes$m_M) 
  state_m_s[state_m_s == "Cov19+"] <- 1
  state_m_s[state_m_s == "CoV19_Dead"] <- 0 
  state_m_s[state_m_s == "O_Causes_Dead"] <- 0 
  state_m_s <- sapply(state_m_s, as.numeric)
  
  Cycle_sum_NT <- 
    as.matrix(colSums(state_m_s)/n_i) # Compute survival probability for each cycle
  
  KM_micro_df[,i] <- Cycle_sum_NT
}

save(KM_micro_df, file = "data/KaplanM_100micro.Rdata")
load("data/KaplanM_100micro.Rdata")

min(KM_micro_df)
#V83
max(KM_micro_df[61,])
#V45

mean_sym <- rowMeans(KM_micro_df)
KM_micro_df$min <- do.call(pmin,KM_micro_df)
KM_micro_df$max <- do.call(pmax,KM_micro_df)

time <- seq(0,60, by = 1)

Covid_rr <- Cov %>% 
  filter(type == 2 | intubated == 1) %>%
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

fit_net_ages_sex <- rs.surv(Surv(time, state) ~ 1, 
                            data = Covid_rr,
                            ratetable = rate_exp_mx, 
                            method = "pohar-perme", 
                            type="kaplan-meier",
                            conf.type="log",
                            conf.int=0.95,
                            rmap = list(age = age, sex = sex, year = diag))

Emp_KM <- c(1, fit_net_ages_sex$surv)
Emp_KM <- head(Emp_KM, 61)

KM_micro_df <- as.data.frame(cbind(time, 
                                   KM_micro_df$min,
                                   KM_micro_df$max,
                                   mean_sym, 
                                   Emp_KM))
KM_micro_df <- KM_micro_df %>% 
  rename(UB = V2,
         LB = V3,
         Microsimulation = mean_sym,
         `Observed cohort` = Emp_KM)

KM_micro_long <- gather(data = KM_micro_df, 
                        key = "KM", 
                        value = "S_Probability", 
                        - c(time, LB, UB))

ggplot(data = KM_micro_long, 
       aes(x = time, 
           y = S_Probability,
           color = KM))+ 
  geom_line(size = 1.1)+
  geom_ribbon(data = filter(KM_micro_long,
                            KM == "Microsimulation"), 
              aes(ymin = LB, 
                  ymax = UB), 
              alpha = 0.02) +
  theme(plot.title = element_text(face = "bold", 
                                  size = 16,
                                  family =),
        plot.caption = element_text(hjust = 0,
                                    colour = "#777777",
                                    size = 10),
        panel.background = element_rect(fill = "white", 
                                        colour = "gray", 
                                        size = 0.15, 
                                        linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, 
                                        linetype = 'solid',
                                        colour = "gray"), 
        axis.text.x = element_text(angle = 90, hjust = 0))+
  scale_x_continuous(breaks = number_ticks(6))+
  # scale_color_manual(values=c ("#11bda3", "#e8e40c", "#e3094a","#113abf", 
  #                              "#e918f0", "#9a18f0" ,"#6309e0","#113abf",
  #                              "#1380bf", "#11b9bf", "#11bda3", "#0fbd71",
  #                              "#0be357", "#5be809", "#a9e309", "#e8e40c")) +
  labs(title = "Kaplan-Meier curves",
       x = "Days",
       y = "Survival Probability",
       color = "Method")

ggsave(paste0("figs/Kaplan_Meier_simulation",
              format(Sys.Date(), "%F"), ".pdf"), 
       width = 7, height = 5)
