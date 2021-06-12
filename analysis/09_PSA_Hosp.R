#### 09 Probabilistic Sensitibity Analysis  ####

# Based in the code developed by the Decision Analysis in R for Technologies in 
# Health (DARTH) workgroup:

####  Load Packages  ####
library(ggplot2)
library(utils)
library(scales)
library(dampack)
library(CEAutil)
library(tidyverse)
library(dampack)
library(data.table)
library(reshape2)
library(survPen)
library(relsurv)
library(readr)
library(readxl)
library(lubridate)
library(muhaz)
library(ggpubr)

####  Clean Work - space  ####
rm(list = ls()) # to clean the work - space

# Load data 
load("data/df_hazards_ICU_hosp.Rdata")

# Load ratetable 
load("data/rate_mx.Rdata")

# Public Data base of suspected people with COVID-19 
load("data/cov_04_04.Rdata") 

# Data with QALEs expenditure for survivals
load("data/df_QALEs_Markov.Rdata")

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

set.seed(29051993)
# Strategy names
v_names_str <- c("No Treatment", "Remdesivir", "Remdesivir & Bariticinib") 

# Number of strategies
n_str <- length(v_names_str)

#### Create df with probabilities ####
d_h_HD_hosp <- df_hazards_ICU_hosp %>% 
  filter(state == "Hosp") %>% 
  filter(time <= 50)

d_h_HD_hosp <- reshape(d_h_HD_hosp,
                       idvar = c("time","group", "sex", "state"),
                       timevar = "Type",
                       direction = "wide")

#### Create Cohorts ####
x <- as.Date("2021-03-29", format = "%Y-%m-%d")

Cohort <- Cov %>% 
  filter(type == 2 & intubated != 1) %>% 
  filter(age >= 45 & age <= 100) %>% 
  select(c("sex", "age", "state", "type","date_death","date_symptoms", 
           "intubated"))%>% 
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

length(Cohort$sex[Cohort$age_range == "65 - 69" & Cohort$sex == "female"])

df_pop_ch <- Cohort %>% 
  select(c("sex", "age"))

Cohort <- Cohort %>% 
  mutate(state = ifelse(intubated == 1, "ICU", "Hosp")) %>%
  select(c("sex","age_range","state","age")) %>%
  rename(group = age_range) %>% 
  mutate(time = 0)

# Create data frame of population from cohort. All begin in day 0
df_X <- Cohort

# List of parameters
l_params <- list( 
  hr_Remd           = 0.73,    # probability to die when healthy
  hr_RemdBa_vs_Remd = 0.65,
  mean_days         = 13.29,
  c_hosp            = 9272,   
  c_remd            = 6188,  
  c_bari            = 3672.35,      
  c_dead            = 0,      # hazard ratio of death in sicker vs healthy
  l_sick            = 1,      # cost of remaining one cycle in the healthy state
  l_dead            = 0       # cost of remaining one cycle in the sick state
)

# Function to generate PSA parameters 
psa_params <- function(n_sim){
  df_parameters <- data.frame(
    # Hazard ratio of Remdesivir 
    hr_Remd      = rlnorm(n_sim,
                              meanlog = log(0.73), 
                              sd = (log(1.03)-log(0.52))/(2*1.96)), 

    # Hazard Ratio of Remdesivir + Ba
    hr_RemdBa_vs_Remd   = rlnorm(n_sim, 
                              meanlog = log(0.65), 
                              sd = (log(1.09)-log(0.39))/(2*1.96)),
    # Mean days
    mean_days  =  rlnorm(n_sim, 
                         meanlog = 2.4321, 
                         sdlog = 0.5153),      
    # Costs
    c_hosp     = rgamma(n_sim, shape = 8, scale = 1159),  
    # cost of remaining one cycle hospitalized 
    c_remd     = 6188,
    # cost of Remdesivir Treatment for one Cycle
    c_bari     = 3672.35,     
    # cost of Baricitinib Treatment for one Cycle
    c_dead     = rep(0, n_sim), # cost of remaining one cycle Dead
    u_sick     = rep(1, n_sim), # utility when Sick 
    u_dead     = rep(0, n_sim)  # utility when Dead
  )
  return(df_parameters)
}

df_params <- psa_params(n_sim)

# Erase data frames that won´t be used anymore
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
 
 mean_days = 13
# Run Markov model on each parameter set of PSA input dataset
 p = Sys.time()
 for(i in 1:n_sim){ #i <- 1
   l_params_psa <- updt_prm_list(l_params, df_params[i,])
   df_out_psa  <- CEA_function(l_params_psa, 
                               df_X = df_X,
                               df_h =  d_h_HD_hosp, 
                               df_pop_ch =  df_pop_ch,
                               df_QALE = df_QALEs, 
                               df_Costs = df_Costs, 
                               mean_days = mean_days)
   df_cost[i, ] <- df_out_psa$Cost
   df_effect[i, ] <- df_out_psa$Effect
   # Display simulation progress
   if(i/(n_sim/10) == round(i/(n_sim/10), 0)) { # display progress every 10%
     cat('\r', paste(i/n_sim * 100, "% done", sep = " "))
   }
 }
 comu_time = Sys.time() - p
 # Time to run
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
                                    df_h = d_h_HD_hosp,
                                    df_pop_ch =  df_pop_ch,
                                    df_QALE = df_QALE, 
                                    df_Costs = df_Costs, 
                                    mean_days = mean_days)
        df_ceas <- c(df_out_psa$Cost, df_out_psa$Effect)
        }
comu_time_par = Sys.time() - p
stopCluster(cl)

df_CEAs_hosp <- as.data.frame(df_ceas)

save(df_CEAs_hosp, file = "data/df_CEAs_hosp_07_06.Rdata")
load(file = "data/df_CEAs_hosp_07_06.Rdata")
df_CEAs_hosp_cost <- df_CEAs_hosp %>%
  select(1:3) %>% 
  rename(`No Treatment` =  V1,
         Remdesivir = V2,
         `Remdesivir & Baricitinib` = V3)

df_CEAs_hosp_cost_long <- gather(df_CEAs_hosp_cost, key = "Strategy", value = "Costs") 

df_CEAs_hosp_effect <- df_CEAs_hosp %>%
  select(4:6) %>% 
  rename(`No Treatment` =  V4,
         Remdesivir = V5,
         `Remdesivir & Baricitinib` = V6)

df_CEAs_hosp_effect_long <- gather(df_CEAs_hosp_effect, key = "Strategy", value = "Effect") 

df_CEAs_hosp_long <- bind_cols(df_CEAs_hosp_cost_long, df_CEAs_hosp_effect_long$Effect) 

df_CEAs_hosp_long <- df_CEAs_hosp_long %>% 
  rename(Effect = ...3)

save(df_CEAs_hosp_cost_long, file = "data/df_CEAs_hosp_cost_long_07_06.Rdata")
save(df_params, file = "data/df_CEAs_params_07_06.Rdata")

# load("data/df_CEAs_hosp_cost_long.Rdata")
#### PSA Object ####
df_costs <- as.data.frame(df_ceas[,1:3])
df_effects <- as.data.frame(df_ceas[,4:6])
v_strategies <- c("No Treatment", "Remdesivir", "Remdesivir.Baricitinib")

psa_obj_hosp <- make_psa_obj(cost = df_costs, 
                             effectiveness = df_effects, 
                             parameters = df_params, 
                             strategies = v_strategies)

n_strategies <- length(v_strategies)

save(df_params, df_costs, df_effects, v_strategies, n_strategies, psa_obj_hosp,
     file = "data/PSA_dataset_07_06_hosp.RData")

load("data/PSA_dataset_07_06_hosp.RData")

psa_obj_hosp$strategies <- c("No treatment", "Remdesivir", "Remdesivir and Baricitinib")

# Compute expected costs and effects for each strategy from the PSA
df_ce_psa <- summary(psa_obj_hosp)

# Calculate incremental cost-effectiveness ratios (ICERs)
df_cea_psa_Hosp <- calculate_icers(cost       = df_ce_psa$meanCost, 
                                   effect     = df_ce_psa$meanEffect,
                                   strategies = df_ce_psa$Strategy)

# Save CEA table with ICERs
save(df_cea_psa_Hosp, file = "data/ICER_results_07_06.RData")

load("data/ICER_results_07_06.RData")
# load("data/ICER_results_29_05.RData")

df_cea_psa_Hosp[1,1] <- "No Treatment"
df_cea_psa_Hosp[2,1] <- "Remdesivir and Baricitinib"
#### CEAC Object & CEAF plot ####

PIB_pc <- 9946*22.10
v_wtp <- seq(0, PIB_pc, by = (PIB_pc)/30)
ceac_obj_hosp <- ceac(wtp = v_wtp, psa = psa_obj_hosp)
summary(ceac_obj_hosp)

plot(ceac_obj_hosp)
## 09.4.3 Expected Loss Curves (ELCs)
elc_obj_hosp <- calc_exp_loss(wtp = v_wtp, psa = psa_obj_hosp)

colnames(elc_obj_hosp)

elc_obj_hosp <- elc_obj_hosp %>% 
  rename(`No Treatment` = No.Treatment,
         `Remdesivir & Baricitinib` = Remdesivir.Baricitinib,
         `Frontier & EVPI` = Frontier_EVPI)

elc_obj_hosp_long <- gather(data = elc_obj_hosp,
                            key = "Strategy", 
                            value = "Value", -WTP) 

# ELC plot
plot(elc_obj_hosp, log_y = FALSE)

#### Deterministic sensitivity analysis
metamodel_icu <- metamodel(analysis = "oneway", psa = psa_obj_icu, wtp = v_wtp )

o_icu <- owsa(sa_obj = psa_obj_icu, outcome = "nmb", wtp = v_wtp)

plot(o_hosp,
     n_x_ticks = 2)

owsa_tornado(o)


tw_hosp <- twsa(psa_obj_hosp, 
           param1 = "hr_RemdBa_vs_Remd", param2 = "hr_Remd", outcome = "nmb", wtp = v_wtp)

tw_hosp$strategy[tw_hosp$strategy == "No.Treatment"] <- "No Treatment"
tw_hosp$strategy[tw_hosp$strategy == "Remdesivir.Baricitinib"] <- "Remdesivir and Baricitinib"

tw_hosp <- tw_hosp %>% 
  rename(`HR Remdesivir` = hr_Remd,
         `HR Remdesivir & Baricitinib` = hr_RemdBa_vs_Remd)

plot(tw_hosp) + 
  theme(plot.title = element_text(face = "bold", 
                                  size = 12,
                                  family =, 
                                  hjust = 0.5),
        axis.line.x = element_line(color="#ebe4e4", size = 1),
        axis.line.y = element_line(color="#ebe4e4", size = 1),
        plot.caption = element_text(hjust = 0,
                                    colour = "#777777",
                                    size = 10),
        panel.background = element_rect(fill = "white", 
                                        colour = "white", 
                                        size = 0.15, 
                                        linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, 
                                        linetype = 'solid',
                                        colour = "#ebe4e4"),
        legend.position = "bottom", 
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8.5)) +
  labs(title = "Two way Sensitivity Analysis",
       fill = "Strategy") +
  scale_fill_manual(values=c ("#575c59", "#2cc9de", "#eb542f"))
  
