#### 12 Probabilistic Sensitibity Analysis patients in ICU  ####

# Based in the code developed by the Decision Analysis in R for Technologies in 
# Health (DARTH) workgroup:

####  Clean Work - space  ####
rm(list = ls()) # to clean the work - space

# Load data 
load("data/df_hazards_ICU_hosp.Rdata")

# Load ratetable 
load("data/rate_mx.Rdata")

# Public Data base of suspected people with COVID-19 
load("data/cov_04_04.Rdata") 

# Data with life expectancy
load("data/Life_expectancy.Rdata")

# Source functions
source("R/Functions_Microsim.R")
source("R/Functions.R")

### Paralelization packages ###
library(foreach)
library(dlookr)
library(doParallel)
library(doSNOW)

#### Filter to get only 2020 ####

Cov <- Cov %>% 
  filter(date_admission >= "2020-03-01" & date_admission <= "2020-12-31")

# Number of simulations
n_sim <- 1000

# Strategy names
v_names_str <- c("No Treatment", "Dexamethasone") 

# Number of strategies
n_str <- length(v_names_str)

#### Create df with probabilities ####
d_h_HD_resp <- df_hazards_ICU_hosp %>% 
  filter(Type != "Observed") %>% 
  filter(state == "ICU")%>% 
  filter(time <= 50)

d_h_HD_resp <- reshape(d_h_HD_resp,
                       idvar = c("time","group", "sex", "state"),
                       timevar = "Type",
                       direction = "wide")

#### Create Cohorts ####
x <- as.Date("2021-03-29", format = "%Y-%m-%d")

Cohort <- Cov %>% 
  filter(intubated == 1) %>% 
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

# Create data frame of population from cohort. All begin in day 0
df_X <- Cohort

# List with parameters
l_params <- list( # El elemento que estamos creando es una lista
  hr_Dexa        = 0.64,  
  c_resp         = 44151,    
  c_dexa         = 4,  
  c_dead         = 0,     
  l_sick         = 1,  # cost of remaining one cycle in the healthy state
  l_dead         = 0   # cost of remaining one cycle in the sick state
)

# Function to generate PSA parameters ()
psa_params <- function(n_sim){
  # set a seed to be able to reproduce the same results
  df_parameters <- data.frame(
    # Hazard ratios of dexamthasone 
    hr_Dexa    = rlnorm(n_sim,
                        meanlog = log(0.64), 
                        sd = (log(0.81)-log(0.51))/(2*1.96)), 
    # Costs
    c_resp     = rgamma(n_sim, shape = 8, scale = 5518.875),  
    c_dexa     = 4,
    # cost of Baricitinib 
    c_dead     = rep(0, n_sim),        # cost of remaining one cycle Dead
    # Utilities
    l_sick     = rep(1, n_sim), # utility when Sick 
    l_dead     = rep(0, n_sim) # utility when Dead
  )
  return(df_parameters)
}

# Data-frame with parameters
df_params <- psa_params(n_sim) 
rm(Cohort,Cov)

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

# Run Markov model on each parameter set of PSA input dataset
p = Sys.time()
for(i in 1:n_sim){ #i <- 1 
  l_params_psa <- updt_prm_list(l_params, df_params[i,])
  df_out_psa  <- CEA_function_icu(l_params_psa, 
                                  df_X = df_X, 
                                  df_h = d_h_HD_resp,
                                  df_pop_ch =  df_pop_ch, 
                                  life_expectancy =  Life_expectancy)
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
df_ceas_icu <- foreach(i = 1:n_sim, .combine = rbind,
                   .export = ls(globalenv()),
                   .packages=c("dampack", 
                               "dplyr", 
                               "data.table", 
                               "reshape2"),
                   .options.snow = options) %dopar% {
                     l_params_psa <- updt_prm_list(l_params,df_params[i,])
                     df_out_psa  <- CEA_function_icu(l_params_psa, 
                                                     df_X = df_X, 
                                                     df_h = d_h_HD_resp,
                                                     df_pop_ch =  df_pop_ch, 
                                                     life_expectancy =  Life_expectancy)
                     df_ceas <- c(df_out_psa$Cost, df_out_psa$Effect)}
comu_time_par = Sys.time() - p
stopCluster(cl)

# Save and data frames

save(df_ceas_icu, file = "data/df_ceas_icu_3.Rdata")
save(df_params, file = "data/df_params_icu_3.Rdata")

load("data/df_CEA_1000_icu.Rdata")
load("data/df_params_icu.Rdata")

df_costs <- as.data.frame(df_ceas_icu[,1:2])

df_costs <- df_costs %>% 
  rename(`No Treatment` = V1,
         `Dexamethasone` = V2)

df_effects <- as.data.frame(df_ceas_icu[,3:4])

df_effects <- df_effects %>% 
  rename(`No Treatment` = V1,
         `Dexamethasone` = V2)

df_costs <- df_costs %>% 
  mutate(name = row_number())

df_effects <- df_effects %>% 
  mutate(name = row_number())

df_effects <- df_effects %>% 
  select(- "name")

df_costs <- df_costs %>% 
  select(- "name")

v_strategies <- c("No Treatment", "Dexamethasone")

psa_obj <- make_psa_obj(cost = df_costs, 
                        effectiveness = df_effects, 
                        parameters = df_params, 
                        strategies = v_strategies)

n_strategies <- length(v_strategies)

save(df_params, df_costs, df_effects, v_strategies, n_strategies, psa_obj,
     file = "data/PSA_dataset_third_trial_icu.RData")

load(file = "data/PSA_dataset_third_trial_icu.RData")

plot(psa_obj)

PIB_pc <- 9946*22.10
PIB_pc_3 <- PIB_pc*3

v_wtp <- seq(0, PIB_pc_3, by = (PIB_pc_3/20))

# Compute expected costs and effects for each strategy from the PSA
df_ce_psa <- summary(psa_obj)

# Calculate incremental cost-effectiveness ratios (ICERs)
df_cea_psa <- calculate_icers(cost       = df_ce_psa$meanCost, 
                              effect     = df_ce_psa$meanEffect,
                              strategies = df_ce_psa$Strategy)

# Save CEA table with ICERs
# As .RData
save(df_cea_psa, 
     file = "data/ICER_results_third_trial_icu.RData")

plot(df_cea_psa)

#ggplot(df_cea_psa, 
#       aes(x = Effect, y = Cost))+
#  geom_line()+
#  geom_point()

ceac_obj <- ceac(wtp = v_wtp, psa = psa_obj)

# Regions of highest probability of cost-effectiveness for each strategy
summary(ceac_obj)
# CEAC & CEAF plot
plot(ceac_obj, currency = "Mexican pesos", 
     txtsize = 11)

## 09.4.3 Expected Loss Curves (ELCs)
elc_obj <- calc_exp_loss(wtp = v_wtp, psa = psa_obj)
# ELC plot
plot(elc_obj, log_y = FALSE)
