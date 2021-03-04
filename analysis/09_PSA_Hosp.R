#### 09 Probabilistic Sensitibity Analysis  ####

# Based in the code developed by the Decision Analysis in R for Technologies in 
# Health (DARTH) workgroup:

# Load data 
load("data/df_hazards_ICU_hosp.Rdata")

# Load ratetable 
load("data/rate_mx.Rdata")

# Public Data base of suspected people with COVID-19 
load("data/cov_09_02.Rdata") 

# Source functions
source("R/Functions_Microsim.R")
source("R/Functions.R")

### Paralelization packages ###
library("foreach")
library(dlookr)
library(doParallel)
library(doSNOW)

n_sim <- 1000
# Function to generate PSA parameters ()
psa_params <- function(n_sim){
  set.seed(02021989) # set a seed to be able to reproduce the same results
  df_parameters <- data.frame(
    # Hazard ratios 
    Remd_effect      = rlnorm(n_sim,
                          meanlog = log(0.73), 
                          sd = (log(1.03)-log(0.52))/(2*1.96)), 
    # Hazard Ratio of Remdesivir
    Remd_Ba_effect   = rlnorm(n_sim, 
                          meanlog = log(0.65), 
                          sd = (log(1.09)-log(0.39))/(2*1.96)),
    # Hazard Ratio of Remdesivir + Ba 
    # Costs
    c_hosp     = rgamma(n_sim, shape = 81.4163, scale = 100),  # cost of remaining one cycle hospitalized 
    c_remd     = rgamma(n_sim, shape = 106.743, scale = 100),
    # cost of Remdesivir Treatment for one Cycle
    c_bari     = rgamma(n_sim, shape = 36.72, scale = 100),     
    # cost of Baricitinib Treatment for one Cycle
    c_dead     = rep(0, n_sim),        # cost of remaining one cycle Dead
    # Utilities
    l_sick     = rep(1, n_sim), # utility when Sick 
    l_dead     = rep(0, n_sim) # utility when Dead                              
  )
  return(df_parameters)
}

# Try it
df_params <- psa_params(n_sim) 

l_params <- list( #El elemento que estamos creando es una lista
  Remd_effect    = 0.73,  # probability to die when healthy
  Remd_Ba_effect = 0.65,   # probability to become sick when healthy
  c_hosp         = 8100,    # probability to become healthy when sick
  c_remd         = 10674.3,  # probability to become sicker when sick
  c_bari         = 3672.35,      # hazard ratio of death in sick vs healthy
  c_dead         = 0,     # hazard ratio of death in sicker vs healthy
  l_sick         = 1,   # cost of remaining one cycle in the healthy state
  l_dead         = 0   # cost of remaining one cycle in the sick state
)

# Strategy names
v_names_str <- c("No Treatment", "Remdesivir", "Remdesivir & Bariticinib") 

# Number of strategies
n_str <- length(v_names_str)

#### Create df with probabilities ####
d_h_HD_hosp <- df_hazards_ICU_hosp %>% 
  filter(Type != "Observed") %>% 
  filter(state == "Hosp")

d_h_HD_hosp <- reshape(d_h_HD_hosp,
                       idvar = c("time","group", "sex", "state"),
                       timevar = "Type",
                       direction = "wide")

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

# Create data frame of population from cohort. All begin in day 0
df_X <- Cohort_hosp

rm(Cohort,Cohort_hosp,Cov)
# Dataframe of costs
df_cost <- as.data.frame(matrix(0, 
                                nrow = n_sim,
                                ncol = n_str))
colnames(df_cost) <- v_names_str
# Dataframe of effectiveness
df_effect <- as.data.frame(matrix(0, 
                                  nrow = n_sim,
                                  ncol = n_str))
colnames(df_effect) <- v_names_str
i <- 1
# Run Markov model on each parameter set of PSA input dataset
p = Sys.time()
for(i in 1:n_sim){
  l_params_psa <- updt_prm_list(l_params, df_params[i,])
  df_out_psa  <- CEA_function(l_params_psa, df_X = df_X, df_h = d_h_HD_hosp)
  df_cost[i, ] <- df_out_psa$Cost
  df_effect[i, ] <- df_out_psa$Effect
  # Display simulation progress
  if(i/(n_sim/10) == round(i/(n_sim/10), 0)) { # display progress every 10%
    cat('\r', paste(i/n_sim * 100, "% done", sep = " "))
  }
}
comu_time = Sys.time() - p

comu_time

#### Microsimulation of 1,000 different sets of parameters ####
#### Parallelization ####
no_cores <- detectCores()-1

# Make cluster object
cl <- makeCluster(no_cores)

# Register clusters
registerDoSNOW(cl)
options <- list(attachExportEnv = TRUE)

p = Sys.time()
df_ceas <- foreach(i = 1:n_sim, .combine = rbind,
        .export = ls(globalenv()),
        .packages=c("dampack", 
                    "dplyr", 
                    "data.table", 
                    "reshape2"),
        .options.snow = options) %dopar% {
        l_params_psa <- updt_prm_list(l_params,df_params[i,])
        df_out_psa  <- CEA_function(l_params_psa, 
                                    df_X = df_X, 
                                    df_h = d_h_HD_hosp)
        df_ceas <- c(df_out_psa$Cost, df_out_psa$Effect)}
comu_time_par = Sys.time() - p
stopCluster(cl)

ls(Cohort)
