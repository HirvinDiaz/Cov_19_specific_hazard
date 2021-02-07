#### Final project script. Author: Hirvin Azael Diaz Zepeda ####

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
library(dampack)
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

######  Load data bases  ######

# Data base of suspected people with COVID-19 published by Mexico´s Ministry of 
# Health 

load("data/cov_24_01.Rdata") 

#### Load functions ####
source("analysis/04_r_functions.R")

#### Load ratetable ####
load("data/rate_mx.Rdata")

#### Compute Hazards for hospitalized population ####

#### Covid_rr_hos_res for intubated individuals higher than 44 years in national DB ####
x <- as.Date(max(Cov$date_admission), format = "%Y-%m-%d")

Covid_rr_resp <- Cov %>% 
  # filter(type == 2) %>% 
  filter(intubated == 1) %>% 
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

#### Eight data frames for each group of study
Covid_rr_resp_male_45 <- Covid_rr_resp %>% 
  filter(sex == 1 & age_range == "45 - 54")

dim(Covid_rr_resp_male_45)

Covid_rr_resp_male_55 <- Covid_rr_resp %>% 
  filter(sex == 1 & age_range == "55 - 64")

dim(Covid_rr_resp_male_55)

Covid_rr_resp_male_70 <- Covid_rr_resp %>% 
  filter(sex == 1 & age_range == "70 +")

dim(Covid_rr_resp_male_70)

Covid_rr_resp_male_65 <- Covid_rr_resp %>% 
  filter(sex == 1 & age_range == "65 - 69")

dim(Covid_rr_resp_male_65)

### Mujeres ####
Covid_rr_resp_female_45 <- Covid_rr_resp %>% 
  filter(sex == 2 & age_range == "45 - 54")

dim(Covid_rr_resp_female_45)

Covid_rr_resp_female_55 <- Covid_rr_resp %>% 
  filter(sex == 2 & age_range == "55 - 64")

dim(Covid_rr_resp_female_55)

Covid_rr_resp_female_70 <- Covid_rr_resp %>% 
  filter(sex == 2 & age_range == "70 +")

dim(Covid_rr_resp_female_70)

Covid_rr_resp_female_65 <- Covid_rr_resp %>% 
  filter(sex == 2 & age_range == "65 - 69")

#### Apply rs.surv function
KM_surv_sex_male_45 <- rs.surv.exp(Surv(time, state) ~ 1,
                                   data = Covid_rr_resp_male_45,
                                   ratetable = rate_exp_mx, 
                                   type="kaplan-meier",
                                   conf.type="log",
                                   conf.int=0.95,
                                   rmap = list(age = age, sex = sex, year = diag))

KM_surv_sex_male_55 <- rs.surv.exp(Surv(time, state) ~ 1,
                                   data = Covid_rr_resp_male_55,
                                   ratetable = rate_exp_mx, 
                                   type="kaplan-meier",
                                   conf.type="log",
                                   conf.int=0.95,
                                   rmap = list(age = age, sex = sex, year = diag))

KM_surv_sex_male_65 <- rs.surv.exp(Surv(time, state) ~ 1,
                                   data = Covid_rr_resp_male_65,
                                   ratetable = rate_exp_mx, 
                                   type="kaplan-meier",
                                   conf.type="log",
                                   conf.int=0.95,
                                   rmap = list(age = age, sex = sex, year = diag))

KM_surv_sex_male_70 <- rs.surv.exp(Surv(time, state) ~ 1,
                                   data = Covid_rr_resp_male_70,
                                   ratetable = rate_exp_mx, 
                                   type="kaplan-meier",
                                   conf.type="log",
                                   conf.int=0.95,
                                   rmap = list(age = age, sex = sex, year = diag))

#### Female ####
KM_surv_sex_female_45 <- rs.surv.exp(Surv(time, state) ~ 1,
                                     data = Covid_rr_resp_female_45,
                                     ratetable = rate_exp_mx, 
                                     type="kaplan-meier",
                                     conf.type="log",
                                     conf.int=0.95,
                                     rmap = list(age = age, sex = sex, year = diag))

KM_surv_sex_female_55 <- rs.surv.exp(Surv(time, state) ~ 1,
                                     data = Covid_rr_resp_female_55,
                                     ratetable = rate_exp_mx, 
                                     type="kaplan-meier",
                                     conf.type="log",
                                     conf.int=0.95,
                                     rmap = list(age = age, sex = sex, year = diag))

KM_surv_sex_female_65 <- rs.surv.exp(Surv(time, state) ~ 1,
                                     data = Covid_rr_resp_female_65,
                                     ratetable = rate_exp_mx, 
                                     type="kaplan-meier",
                                     conf.type="log",
                                     conf.int=0.95,
                                     rmap = list(age = age, sex = sex, year = diag))

KM_surv_sex_female_70 <- rs.surv.exp(Surv(time, state) ~ 1,
                                     data = Covid_rr_resp_female_70,
                                     ratetable = rate_exp_mx, 
                                     type="kaplan-meier",
                                     conf.type="log",
                                     conf.int=0.95,
                                     rmap = list(age = age, sex = sex, year = diag))
#### Data Frame with hazards by each group Hospitalized ####

# Males

# 45
Covid_19 <- head(KM_surv_sex_male_45$hazard_exc, 60)
Observed <- head(KM_surv_sex_male_45$hazard_obs, 60)
Population <- head(KM_surv_sex_male_45$hazard_pop, 60)
time <- seq(1,60, by = 1)

df_hazard_exp_male_45 <- data.frame(time, 
                                    Covid_19, 
                                    Observed, 
                                    Population)

df_hazard_exp_long_45 <- gather(data = df_hazard_exp_male_45, 
                                key = "Type", 
                                value = "Hazard",
                                -time)

df_hazard_exp_long_45 <- df_hazard_exp_long_45 %>% 
  mutate(group = "45 - 54")

# 55
Covid_19 <- head(KM_surv_sex_male_55$hazard_exc, 60)
Observed <- head(KM_surv_sex_male_55$hazard_obs, 60)
Population <- head(KM_surv_sex_male_55$hazard_pop, 60)
time <- seq(1,60, by = 1)

df_hazard_exp_male_55 <- data.frame(time, 
                                    Covid_19, 
                                    Observed, 
                                    Population)

df_hazard_exp_long_55 <- gather(data = df_hazard_exp_male_55, 
                                key = "Type", 
                                value = "Hazard",
                                -time)

df_hazard_exp_long_55 <- df_hazard_exp_long_55 %>% 
  mutate(group = "55 - 64")

# 65
Covid_19 <- head(KM_surv_sex_male_65$hazard_exc, 60)
Observed <- head(KM_surv_sex_male_65$hazard_obs, 60)
Population <- head(KM_surv_sex_male_65$hazard_pop, 60)
time <- seq(1,60, by = 1)

df_hazard_exp_male_65 <- data.frame(time, 
                                    Covid_19, 
                                    Observed, 
                                    Population)

df_hazard_exp_long_65 <- gather(data = df_hazard_exp_male_65, 
                                key = "Type", 
                                value = "Hazard",
                                -time)

df_hazard_exp_long_65 <- df_hazard_exp_long_65 %>% 
  mutate(group = "65 - 69")

# 70
Covid_19 <- head(KM_surv_sex_male_70$hazard_exc, 60)
Observed <- head(KM_surv_sex_male_70$hazard_obs, 60)
Population <- head(KM_surv_sex_male_70$hazard_pop, 60)
time <- seq(1,60, by = 1)

df_hazard_exp_male_70 <- data.frame(time, 
                                    Covid_19, 
                                    Observed, 
                                    Population)

df_hazard_exp_long_70 <- gather(data = df_hazard_exp_male_70, 
                                key = "Type", 
                                value = "Hazard",
                                -time)

df_hazard_exp_long_70 <- df_hazard_exp_long_70 %>% 
  mutate(group = "70 +")

df_hazard_long_male <- rbind(df_hazard_exp_long_45, df_hazard_exp_long_55,
                             df_hazard_exp_long_65, df_hazard_exp_long_70)

df_hazard_long_male <- df_hazard_long_male %>% 
  mutate(sex = "male")

gg_males_hazard <- ggplot(df_hazard_long_male, 
                          aes(x = time, 
                              y = Hazard,
                              color = Type)) +
  geom_line(size = 1)+
  facet_wrap(~group, scales = "free") +
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
                                        colour = "gray"))+ 
  # axis.text.x = element_text(angle = 90, hjust = 0))+
  scale_x_continuous(breaks = number_ticks(4))+
  scale_y_continuous(breaks = number_ticks(3))+
  scale_color_manual(values=c ("#2e7bff", "#e6ad02", "#02ebdb")) +
  labs(title = "Intubated patients: Male age groups",
       x = " Days ",
       y = " ",
       color = "Hazard")

# Females
# 45
Covid_19 <- head(KM_surv_sex_female_45$hazard_exc, 60)
Observed <- head(KM_surv_sex_female_45$hazard_obs, 60)
Population <- head(KM_surv_sex_female_45$hazard_pop, 60)
time <- seq(1,60, by = 1)

df_hazard_exp_female_45 <- data.frame(time, 
                                      Covid_19, 
                                      Observed, 
                                      Population)

df_hazard_exp_long_45 <- gather(data = df_hazard_exp_female_45, 
                                key = "Type", 
                                value = "Hazard",
                                -time)

df_hazard_exp_long_45 <- df_hazard_exp_long_45 %>% 
  mutate(group = "45 - 54")

# 55
Covid_19 <- head(KM_surv_sex_female_55$hazard_exc, 60)
Observed <- head(KM_surv_sex_female_55$hazard_obs, 60)
Population <- head(KM_surv_sex_female_55$hazard_pop, 60)
time <- seq(1,60, by = 1)

df_hazard_exp_female_55 <- data.frame(time, 
                                      Covid_19, 
                                      Observed, 
                                      Population)

df_hazard_exp_long_55 <- gather(data = df_hazard_exp_female_55, 
                                key = "Type", 
                                value = "Hazard",
                                -time)

df_hazard_exp_long_55 <- df_hazard_exp_long_55 %>% 
  mutate(group = "55 - 64")

# 65
Covid_19 <- head(KM_surv_sex_female_65$hazard_exc, 60)
Observed <- head(KM_surv_sex_female_65$hazard_obs, 60)
Population <- head(KM_surv_sex_female_65$hazard_pop, 60)
time <- seq(1,60, by = 1)

df_hazard_exp_female_65 <- data.frame(time, 
                                      Covid_19, 
                                      Observed, 
                                      Population)

df_hazard_exp_long_65 <- gather(data = df_hazard_exp_female_65, 
                                key = "Type", 
                                value = "Hazard",
                                -time)

df_hazard_exp_long_65 <- df_hazard_exp_long_65 %>% 
  mutate(group = "65 - 69")

# 70
Covid_19 <- head(KM_surv_sex_female_70$hazard_exc, 60)
Observed <- head(KM_surv_sex_female_70$hazard_obs, 60)
Population <- head(KM_surv_sex_female_70$hazard_pop, 60)
time <- seq(1,60, by = 1)

df_hazard_exp_female_70 <- data.frame(time, 
                                      Covid_19, 
                                      Observed, 
                                      Population)

df_hazard_exp_long_70 <- gather(data = df_hazard_exp_female_70, 
                                key = "Type", 
                                value = "Hazard",
                                -time)

df_hazard_exp_long_70 <- df_hazard_exp_long_70 %>% 
  mutate(group = "70 +")

df_hazard_long_female <- rbind(df_hazard_exp_long_45, df_hazard_exp_long_55,
                               df_hazard_exp_long_65, df_hazard_exp_long_70)

df_hazard_long_female <- df_hazard_long_female %>% 
  mutate(sex = "female")

gg_females_hazard <- ggplot(df_hazard_long_female, 
                            aes(x = time, 
                                y = Hazard,
                                color = Type)) +
  geom_line(size = 1)+
  facet_wrap(~group, scales = "free") +
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
                                        colour = "gray"))+
  scale_x_continuous(breaks = number_ticks(4))+
  scale_y_continuous(breaks = number_ticks(3))+
  # axis.text.x = element_text(angle = 90, hjust = 0))+
  scale_color_manual(values=c ("#2e7bff", "#e6ad02", "#02ebdb")) +
  labs(title = "Intubated patients: Female age groups",
       x = " ",
       y = " ",
       color = "Hazard")

ggarrange(gg_females_hazard, gg_males_hazard,
          #labels        = c("Females", "Males"),
          ncol          = 1, 
          common.legend = TRUE,
          label.y = c(1,0),
          legend        = "bottom")

ggsave(paste0("figs/Hazard_Intubated_age_groups",
              format(Sys.Date(), "%F"), ".pdf"), 
       width = 7, height = 5)

df_hazards_intubated <- rbind(df_hazard_long_female, df_hazard_long_male)
save(df_hazards_intubated, file = "data/hazard_intubated.Rdata")

### Hospitalized ###

#### Covid_rr_hos_res for all individuals higher than 44 years in national DB ####
x <- as.Date(max(Cov$date_admission), format = "%Y-%m-%d")

Covid_rr_hos_res_hos_res <- Cov %>% 
  filter(type == 2) %>% 
  #  filter(intubated == 1) %>% 
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

#### Eight data frames for each group of study
Covid_rr_hos_res_male_45 <- Covid_rr_hos_res %>% 
  filter(sex == 1 & age_range == "45 - 54")

dim(Covid_rr_hos_res_male_45)

Covid_rr_hos_res_male_55 <- Covid_rr_hos_res %>% 
  filter(sex == 1 & age_range == "55 - 64")

dim(Covid_rr_hos_res_male_55)

Covid_rr_hos_res_male_70 <- Covid_rr_hos_res %>% 
  filter(sex == 1 & age_range == "70 +")

dim(Covid_rr_hos_res_male_70)

Covid_rr_hos_res_male_65 <- Covid_rr_hos_res %>% 
  filter(sex == 1 & age_range == "65 - 69")

dim(Covid_rr_hos_res_male_65)

### Mujeres ####
Covid_rr_hos_res_female_45 <- Covid_rr_hos_res %>% 
  filter(sex == 2 & age_range == "45 - 54")

dim(Covid_rr_hos_res_female_45)

Covid_rr_hos_res_female_55 <- Covid_rr_hos_res %>% 
  filter(sex == 2 & age_range == "55 - 64")

dim(Covid_rr_hos_res_female_55)

Covid_rr_hos_res_female_70 <- Covid_rr_hos_res %>% 
  filter(sex == 2 & age_range == "70 +")

dim(Covid_rr_hos_res_female_70)

Covid_rr_hos_res_female_65 <- Covid_rr_hos_res %>% 
  filter(sex == 2 & age_range == "65 - 69")

#### Apply rs.surv function
KM_surv_sex_male_45 <- rs.surv.exp(Surv(time, state) ~ 1,
                                   data = Covid_rr_hos_res_male_45,
                                   ratetable = rate_exp_mx, 
                                   type="kaplan-meier",
                                   conf.type="log",
                                   conf.int=0.95,
                                   rmap = list(age = age, sex = sex, year = diag))

KM_surv_sex_male_55 <- rs.surv.exp(Surv(time, state) ~ 1,
                                   data = Covid_rr_hos_res_male_55,
                                   ratetable = rate_exp_mx, 
                                   type="kaplan-meier",
                                   conf.type="log",
                                   conf.int=0.95,
                                   rmap = list(age = age, sex = sex, year = diag))

KM_surv_sex_male_65 <- rs.surv.exp(Surv(time, state) ~ 1,
                                   data = Covid_rr_hos_res_male_65,
                                   ratetable = rate_exp_mx, 
                                   type="kaplan-meier",
                                   conf.type="log",
                                   conf.int=0.95,
                                   rmap = list(age = age, sex = sex, year = diag))

KM_surv_sex_male_70 <- rs.surv.exp(Surv(time, state) ~ 1,
                                   data = Covid_rr_hos_res_male_70,
                                   ratetable = rate_exp_mx, 
                                   type="kaplan-meier",
                                   conf.type="log",
                                   conf.int=0.95,
                                   rmap = list(age = age, sex = sex, year = diag))

#### Female ####
KM_surv_sex_female_45 <- rs.surv.exp(Surv(time, state) ~ 1,
                                     data = Covid_rr_hos_res_female_45,
                                     ratetable = rate_exp_mx, 
                                     type="kaplan-meier",
                                     conf.type="log",
                                     conf.int=0.95,
                                     rmap = list(age = age, sex = sex, year = diag))

KM_surv_sex_female_55 <- rs.surv.exp(Surv(time, state) ~ 1,
                                     data = Covid_rr_hos_res_female_55,
                                     ratetable = rate_exp_mx, 
                                     type="kaplan-meier",
                                     conf.type="log",
                                     conf.int=0.95,
                                     rmap = list(age = age, sex = sex, year = diag))

KM_surv_sex_female_65 <- rs.surv.exp(Surv(time, state) ~ 1,
                                     data = Covid_rr_hos_res_female_65,
                                     ratetable = rate_exp_mx, 
                                     type="kaplan-meier",
                                     conf.type="log",
                                     conf.int=0.95,
                                     rmap = list(age = age, sex = sex, year = diag))

KM_surv_sex_female_70 <- rs.surv.exp(Surv(time, state) ~ 1,
                                     data = Covid_rr_hos_res_female_70,
                                     ratetable = rate_exp_mx, 
                                     type="kaplan-meier",
                                     conf.type="log",
                                     conf.int=0.95,
                                     rmap = list(age = age, sex = sex, year = diag))
#### Data Frame with hazards by each group Hospitalized ####

# Males

# 45
Covid_19 <- head(KM_surv_sex_male_45$hazard_exc, 60)
Observed <- head(KM_surv_sex_male_45$hazard_obs, 60)
Population <- head(KM_surv_sex_male_45$hazard_pop, 60)
time <- seq(1,60, by = 1)

df_hazard_exp_male_45 <- data.frame(time, 
                                    Covid_19, 
                                    Observed, 
                                    Population)

df_hazard_exp_long_45 <- gather(data = df_hazard_exp_male_45, 
                                key = "Type", 
                                value = "Hazard",
                                -time)

df_hazard_exp_long_45 <- df_hazard_exp_long_45 %>% 
  mutate(group = "45 - 54")

# 55
Covid_19 <- head(KM_surv_sex_male_55$hazard_exc, 60)
Observed <- head(KM_surv_sex_male_55$hazard_obs, 60)
Population <- head(KM_surv_sex_male_55$hazard_pop, 60)
time <- seq(1,60, by = 1)

df_hazard_exp_male_55 <- data.frame(time, 
                                    Covid_19, 
                                    Observed, 
                                    Population)

df_hazard_exp_long_55 <- gather(data = df_hazard_exp_male_55, 
                                key = "Type", 
                                value = "Hazard",
                                -time)

df_hazard_exp_long_55 <- df_hazard_exp_long_55 %>% 
  mutate(group = "55 - 64")

# 65
Covid_19 <- head(KM_surv_sex_male_65$hazard_exc, 60)
Observed <- head(KM_surv_sex_male_65$hazard_obs, 60)
Population <- head(KM_surv_sex_male_65$hazard_pop, 60)
time <- seq(1,60, by = 1)

df_hazard_exp_male_65 <- data.frame(time, 
                                    Covid_19, 
                                    Observed, 
                                    Population)

df_hazard_exp_long_65 <- gather(data = df_hazard_exp_male_65, 
                                key = "Type", 
                                value = "Hazard",
                                -time)

df_hazard_exp_long_65 <- df_hazard_exp_long_65 %>% 
  mutate(group = "65 - 69")

# 70
Covid_19 <- head(KM_surv_sex_male_70$hazard_exc, 60)
Observed <- head(KM_surv_sex_male_70$hazard_obs, 60)
Population <- head(KM_surv_sex_male_70$hazard_pop, 60)
time <- seq(1,60, by = 1)

df_hazard_exp_male_70 <- data.frame(time, 
                                    Covid_19, 
                                    Observed, 
                                    Population)

df_hazard_exp_long_70 <- gather(data = df_hazard_exp_male_70, 
                                key = "Type", 
                                value = "Hazard",
                                -time)

df_hazard_exp_long_70 <- df_hazard_exp_long_70 %>% 
  mutate(group = "70 +")

df_hazard_long_male <- rbind(df_hazard_exp_long_45, df_hazard_exp_long_55,
                             df_hazard_exp_long_65, df_hazard_exp_long_70)

df_hazard_long_male <- df_hazard_long_male %>% 
  mutate(sex = "male")

gg_males_hazard <- ggplot(df_hazard_long_male, 
                          aes(x = time, 
                              y = Hazard,
                              color = Type)) +
  geom_line(size = 1)+
  facet_wrap(~group, scales = "free") +
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
                                        colour = "gray"))+ 
  # axis.text.x = element_text(angle = 90, hjust = 0))+
  scale_x_continuous(breaks = number_ticks(4))+
  scale_y_continuous(breaks = number_ticks(3))+
  scale_color_manual(values=c ("#2e7bff", "#e6ad02", "#02ebdb")) +
  labs(title = "Hospitalized: Male age groups",
       x = " Days ",
       y = " ",
       color = "Hazard")

# Females
# 45
Covid_19 <- head(KM_surv_sex_female_45$hazard_exc, 60)
Observed <- head(KM_surv_sex_female_45$hazard_obs, 60)
Population <- head(KM_surv_sex_female_45$hazard_pop, 60)
time <- seq(1,60, by = 1)

df_hazard_exp_female_45 <- data.frame(time, 
                                      Covid_19, 
                                      Observed, 
                                      Population)

df_hazard_exp_long_45 <- gather(data = df_hazard_exp_female_45, 
                                key = "Type", 
                                value = "Hazard",
                                -time)

df_hazard_exp_long_45 <- df_hazard_exp_long_45 %>% 
  mutate(group = "45 - 54")

# 55
Covid_19 <- head(KM_surv_sex_female_55$hazard_exc, 60)
Observed <- head(KM_surv_sex_female_55$hazard_obs, 60)
Population <- head(KM_surv_sex_female_55$hazard_pop, 60)
time <- seq(1,60, by = 1)

df_hazard_exp_female_55 <- data.frame(time, 
                                      Covid_19, 
                                      Observed, 
                                      Population)

df_hazard_exp_long_55 <- gather(data = df_hazard_exp_female_55, 
                                key = "Type", 
                                value = "Hazard",
                                -time)

df_hazard_exp_long_55 <- df_hazard_exp_long_55 %>% 
  mutate(group = "55 - 64")

# 65
Covid_19 <- head(KM_surv_sex_female_65$hazard_exc, 60)
Observed <- head(KM_surv_sex_female_65$hazard_obs, 60)
Population <- head(KM_surv_sex_female_65$hazard_pop, 60)
time <- seq(1,60, by = 1)

df_hazard_exp_female_65 <- data.frame(time, 
                                      Covid_19, 
                                      Observed, 
                                      Population)

df_hazard_exp_long_65 <- gather(data = df_hazard_exp_female_65, 
                                key = "Type", 
                                value = "Hazard",
                                -time)

df_hazard_exp_long_65 <- df_hazard_exp_long_65 %>% 
  mutate(group = "65 - 69")

# 70
Covid_19 <- head(KM_surv_sex_female_70$hazard_exc, 60)
Observed <- head(KM_surv_sex_female_70$hazard_obs, 60)
Population <- head(KM_surv_sex_female_70$hazard_pop, 60)
time <- seq(1,60, by = 1)

df_hazard_exp_female_70 <- data.frame(time, 
                                      Covid_19, 
                                      Observed, 
                                      Population)

df_hazard_exp_long_70 <- gather(data = df_hazard_exp_female_70, 
                                key = "Type", 
                                value = "Hazard",
                                -time)

df_hazard_exp_long_70 <- df_hazard_exp_long_70 %>% 
  mutate(group = "70 +")

df_hazard_long_female <- rbind(df_hazard_exp_long_45, df_hazard_exp_long_55,
                               df_hazard_exp_long_65, df_hazard_exp_long_70)

df_hazard_long_female <- df_hazard_long_female %>% 
  mutate(sex = "female")

gg_females_hazard <- ggplot(df_hazard_long_female, 
                            aes(x = time, 
                                y = Hazard,
                                color = Type)) +
  geom_line(size = 1)+
  facet_wrap(~group, scales = "free") +
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
                                        colour = "gray"))+
  scale_x_continuous(breaks = number_ticks(4))+
  scale_y_continuous(breaks = number_ticks(3))+
  # axis.text.x = element_text(angle = 90, hjust = 0))+
  scale_color_manual(values=c ("#2e7bff", "#e6ad02", "#02ebdb")) +
  labs(title = "Hospitalized: Female age groups",
       x = " ",
       y = " ",
       color = "Hazard")

ggarrange(gg_females_hazard, gg_males_hazard,
          #labels        = c("Females", "Males"),
          ncol          = 1, 
          common.legend = TRUE,
          label.y = c(1,0),
          legend        = "bottom")

ggsave(paste0("figs/Hazard_hos_resitalized_age_groups",
              format(Sys.Date(), "%F"), ".pdf"), 
       width = 7, height = 5)

df_hazards_hos_res <- rbind(df_hazard_long_female, df_hazard_long_male)
save(df_hazards_hos_res, file = "data/hazard_hos_res.Rdata")


### Hospitalized and intubated ###

#### Covid_rr_hos_res DB ####
x <- as.Date(max(Cov$date_admission), format = "%Y-%m-%d")

Covid_rr_hos_res <- Cov %>% 
  filter(type == 2 |intubated == 1) %>% 
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

#### Eight data frames for each group of study
Covid_rr_hos_res_male_45 <- Covid_rr_hos_res %>% 
  filter(sex == 1 & age_range == "45 - 54")

dim(Covid_rr_hos_res_male_45)

Covid_rr_hos_res_male_55 <- Covid_rr_hos_res %>% 
  filter(sex == 1 & age_range == "55 - 64")

dim(Covid_rr_hos_res_male_55)

Covid_rr_hos_res_male_70 <- Covid_rr_hos_res %>% 
  filter(sex == 1 & age_range == "70 +")

dim(Covid_rr_hos_res_male_70)

Covid_rr_hos_res_male_65 <- Covid_rr_hos_res %>% 
  filter(sex == 1 & age_range == "65 - 69")

dim(Covid_rr_hos_res_male_65)

### Mujeres ####
Covid_rr_hos_res_female_45 <- Covid_rr_hos_res %>% 
  filter(sex == 2 & age_range == "45 - 54")

dim(Covid_rr_hos_res_female_45)

Covid_rr_hos_res_female_55 <- Covid_rr_hos_res %>% 
  filter(sex == 2 & age_range == "55 - 64")

dim(Covid_rr_hos_res_female_55)

Covid_rr_hos_res_female_70 <- Covid_rr_hos_res %>% 
  filter(sex == 2 & age_range == "70 +")

dim(Covid_rr_hos_res_female_70)

Covid_rr_hos_res_female_65 <- Covid_rr_hos_res %>% 
  filter(sex == 2 & age_range == "65 - 69")

#### Apply rs.surv function
KM_surv_sex_male_45 <- rs.surv.exp(Surv(time, state) ~ 1,
                                   data = Covid_rr_hos_res_male_45,
                                   ratetable = rate_exp_mx, 
                                   type="kaplan-meier",
                                   conf.type="log",
                                   conf.int=0.95,
                                   rmap = list(age = age, sex = sex, year = diag))

KM_surv_sex_male_55 <- rs.surv.exp(Surv(time, state) ~ 1,
                                   data = Covid_rr_hos_res_male_55,
                                   ratetable = rate_exp_mx, 
                                   type="kaplan-meier",
                                   conf.type="log",
                                   conf.int=0.95,
                                   rmap = list(age = age, sex = sex, year = diag))

KM_surv_sex_male_65 <- rs.surv.exp(Surv(time, state) ~ 1,
                                   data = Covid_rr_hos_res_male_65,
                                   ratetable = rate_exp_mx, 
                                   type="kaplan-meier",
                                   conf.type="log",
                                   conf.int=0.95,
                                   rmap = list(age = age, sex = sex, year = diag))

KM_surv_sex_male_70 <- rs.surv.exp(Surv(time, state) ~ 1,
                                   data = Covid_rr_hos_res_male_70,
                                   ratetable = rate_exp_mx, 
                                   type="kaplan-meier",
                                   conf.type="log",
                                   conf.int=0.95,
                                   rmap = list(age = age, sex = sex, year = diag))

#### Female ####
KM_surv_sex_female_45 <- rs.surv.exp(Surv(time, state) ~ 1,
                                     data = Covid_rr_hos_res_female_45,
                                     ratetable = rate_exp_mx, 
                                     type="kaplan-meier",
                                     conf.type="log",
                                     conf.int=0.95,
                                     rmap = list(age = age, sex = sex, year = diag))

KM_surv_sex_female_55 <- rs.surv.exp(Surv(time, state) ~ 1,
                                     data = Covid_rr_hos_res_female_55,
                                     ratetable = rate_exp_mx, 
                                     type="kaplan-meier",
                                     conf.type="log",
                                     conf.int=0.95,
                                     rmap = list(age = age, sex = sex, year = diag))

KM_surv_sex_female_65 <- rs.surv.exp(Surv(time, state) ~ 1,
                                     data = Covid_rr_hos_res_female_65,
                                     ratetable = rate_exp_mx, 
                                     type="kaplan-meier",
                                     conf.type="log",
                                     conf.int=0.95,
                                     rmap = list(age = age, sex = sex, year = diag))

KM_surv_sex_female_70 <- rs.surv.exp(Surv(time, state) ~ 1,
                                     data = Covid_rr_hos_res_female_70,
                                     ratetable = rate_exp_mx, 
                                     type="kaplan-meier",
                                     conf.type="log",
                                     conf.int=0.95,
                                     rmap = list(age = age, sex = sex, year = diag))
#### Data Frame with hazards by each group Hospitalized ####

# Males

# 45
Covid_19 <- head(KM_surv_sex_male_45$hazard_exc, 60)
Observed <- head(KM_surv_sex_male_45$hazard_obs, 60)
Population <- head(KM_surv_sex_male_45$hazard_pop, 60)
time <- seq(1,60, by = 1)

df_hazard_exp_male_45 <- data.frame(time, 
                                    Covid_19, 
                                    Observed, 
                                    Population)

df_hazard_exp_long_45 <- gather(data = df_hazard_exp_male_45, 
                                key = "Type", 
                                value = "Hazard",
                                -time)

df_hazard_exp_long_45 <- df_hazard_exp_long_45 %>% 
  mutate(group = "45 - 54")

# 55
Covid_19 <- head(KM_surv_sex_male_55$hazard_exc, 60)
Observed <- head(KM_surv_sex_male_55$hazard_obs, 60)
Population <- head(KM_surv_sex_male_55$hazard_pop, 60)
time <- seq(1,60, by = 1)

df_hazard_exp_male_55 <- data.frame(time, 
                                    Covid_19, 
                                    Observed, 
                                    Population)

df_hazard_exp_long_55 <- gather(data = df_hazard_exp_male_55, 
                                key = "Type", 
                                value = "Hazard",
                                -time)

df_hazard_exp_long_55 <- df_hazard_exp_long_55 %>% 
  mutate(group = "55 - 64")

# 65
Covid_19 <- head(KM_surv_sex_male_65$hazard_exc, 60)
Observed <- head(KM_surv_sex_male_65$hazard_obs, 60)
Population <- head(KM_surv_sex_male_65$hazard_pop, 60)
time <- seq(1,60, by = 1)

df_hazard_exp_male_65 <- data.frame(time, 
                                    Covid_19, 
                                    Observed, 
                                    Population)

df_hazard_exp_long_65 <- gather(data = df_hazard_exp_male_65, 
                                key = "Type", 
                                value = "Hazard",
                                -time)

df_hazard_exp_long_65 <- df_hazard_exp_long_65 %>% 
  mutate(group = "65 - 69")

# 70
Covid_19 <- head(KM_surv_sex_male_70$hazard_exc, 60)
Observed <- head(KM_surv_sex_male_70$hazard_obs, 60)
Population <- head(KM_surv_sex_male_70$hazard_pop, 60)
time <- seq(1,60, by = 1)

df_hazard_exp_male_70 <- data.frame(time, 
                                    Covid_19, 
                                    Observed, 
                                    Population)

df_hazard_exp_long_70 <- gather(data = df_hazard_exp_male_70, 
                                key = "Type", 
                                value = "Hazard",
                                -time)

df_hazard_exp_long_70 <- df_hazard_exp_long_70 %>% 
  mutate(group = "70 +")

df_hazard_long_male <- rbind(df_hazard_exp_long_45, df_hazard_exp_long_55,
                             df_hazard_exp_long_65, df_hazard_exp_long_70)

df_hazard_long_male <- df_hazard_long_male %>% 
  mutate(sex = "male")

gg_males_hazard <- ggplot(df_hazard_long_male, 
                          aes(x = time, 
                              y = Hazard,
                              color = Type)) +
  geom_line(size = 1)+
  # geom_point(size  = 0.)+
  facet_wrap(~group, scales = "free") +
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
                                        colour = "gray"))+ 
  # axis.text.x = element_text(angle = 90, hjust = 0))+
  scale_x_continuous(breaks = number_ticks(4))+
  scale_y_continuous(breaks = number_ticks(3))+
  scale_color_manual(values=c ("#2e7bff", "#e6ad02", "#02ebdb")) +
  labs(title = "Hospitalized: Male",
       x = " Days ",
       y = " ",
       color = "Hazard")

ggsave(paste0("figs/Hazard_hos_res_male",
              format(Sys.Date(), "%F"), ".pdf"), 
       width = 7, height = 5)


# Females
# 45
Covid_19 <- head(KM_surv_sex_female_45$hazard_exc, 60)
Observed <- head(KM_surv_sex_female_45$hazard_obs, 60)
Population <- head(KM_surv_sex_female_45$hazard_pop, 60)
time <- seq(1,60, by = 1)

df_hazard_exp_female_45 <- data.frame(time, 
                                      Covid_19, 
                                      Observed, 
                                      Population)

df_hazard_exp_long_45 <- gather(data = df_hazard_exp_female_45, 
                                key = "Type", 
                                value = "Hazard",
                                -time)

df_hazard_exp_long_45 <- df_hazard_exp_long_45 %>% 
  mutate(group = "45 - 54")

# 55
Covid_19 <- head(KM_surv_sex_female_55$hazard_exc, 60)
Observed <- head(KM_surv_sex_female_55$hazard_obs, 60)
Population <- head(KM_surv_sex_female_55$hazard_pop, 60)
time <- seq(1,60, by = 1)

df_hazard_exp_female_55 <- data.frame(time, 
                                      Covid_19, 
                                      Observed, 
                                      Population)

df_hazard_exp_long_55 <- gather(data = df_hazard_exp_female_55, 
                                key = "Type", 
                                value = "Hazard",
                                -time)

df_hazard_exp_long_55 <- df_hazard_exp_long_55 %>% 
  mutate(group = "55 - 64")

# 65
Covid_19 <- head(KM_surv_sex_female_65$hazard_exc, 60)
Observed <- head(KM_surv_sex_female_65$hazard_obs, 60)
Population <- head(KM_surv_sex_female_65$hazard_pop, 60)
time <- seq(1,60, by = 1)

df_hazard_exp_female_65 <- data.frame(time, 
                                      Covid_19, 
                                      Observed, 
                                      Population)

df_hazard_exp_long_65 <- gather(data = df_hazard_exp_female_65, 
                                key = "Type", 
                                value = "Hazard",
                                -time)

df_hazard_exp_long_65 <- df_hazard_exp_long_65 %>% 
  mutate(group = "65 - 69")

# 70
Covid_19 <- head(KM_surv_sex_female_70$hazard_exc, 60)
Observed <- head(KM_surv_sex_female_70$hazard_obs, 60)
Population <- head(KM_surv_sex_female_70$hazard_pop, 60)
time <- seq(1,60, by = 1)

df_hazard_exp_female_70 <- data.frame(time, 
                                      Covid_19, 
                                      Observed, 
                                      Population)

df_hazard_exp_long_70 <- gather(data = df_hazard_exp_female_70, 
                                key = "Type", 
                                value = "Hazard",
                                -time)

df_hazard_exp_long_70 <- df_hazard_exp_long_70 %>% 
  mutate(group = "70 +")

df_hazard_long_female <- rbind(df_hazard_exp_long_45, df_hazard_exp_long_55,
                               df_hazard_exp_long_65, df_hazard_exp_long_70)

df_hazard_long_female <- df_hazard_long_female %>% 
  mutate(sex = "female")

gg_females_hazard <- ggplot(df_hazard_long_female, 
                            aes(x = time, 
                                y = Hazard,
                                color = Type)) +
  geom_line(size = 1)+
  facet_wrap(~group, scales = "free") +
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
                                        colour = "gray"))+
  scale_x_continuous(breaks = number_ticks(4))+
  scale_y_continuous(breaks = number_ticks(3))+
  # axis.text.x = element_text(angle = 90, hjust = 0))+
  scale_color_manual(values=c ("#2e7bff", "#e6ad02", "#02ebdb")) +
  labs(title = "Hospitalized: Female",
       x = " ",
       y = " ",
       color = "Hazard")

ggarrange(gg_females_hazard, gg_males_hazard,
          #labels        = c("Females", "Males"),
          ncol          = 2, 
          common.legend = TRUE,
          label.y = c(1,0),
          legend        = "bottom")

ggsave(paste0("figs/Hazard_hos_resitalized_age_groups",
              format(Sys.Date(), "%F"), ".pdf"), 
       width = 7, height = 5)

df_hazards_hos_res <- rbind(df_hazard_long_female, df_hazard_long_male)
save(df_hazards_hos_res, file = "data/hazard_hos_res.Rdata")

load("data/hazard_hos_res.Rdata")

df_hazards_hos_res$Type[df_hazards_hos_res$Type == "Covid_19"] <- "COVID-19"
df_hazards_hos_res$sex[df_hazards_hos_res$sex == "male"] <- "Male"
df_hazards_hos_res$sex[df_hazards_hos_res$sex == "female"] <- "Female"

df_hazards_hos_res <- df_hazards_hos_res %>% 
  mutate(age_group = paste(group, sex))

ggplot(filter(df_hazards_hos_res,
              Type != "Observed"), 
       aes(x = time, 
           y = Hazard,
           color = Type)) +
  geom_line(size = 0.8)+
  geom_line(data = filter(df_hazards_hos_res,
                          Type == "Observed"), 
            aes(x = time, 
                y = Hazard),
            alpha = 0.4, 
            size = 2)+
  facet_wrap(~age_group, ncol = 2) +
  theme(plot.title = element_text(face = "bold", 
                                  size = 14,
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
                                        colour = "#6e6866"))+
  scale_x_continuous(breaks = number_ticks(4))+
  scale_y_continuous(breaks = number_ticks(3))+
  # axis.text.x = element_text(angle = 90, hjust = 0))+
  scale_color_manual(values=c ("#120d0c", "#e32402", "#02ebdb")) +
  labs(title = "Specific Covid-19, background population and observed hazard",
       x = " ",
       y = " ",
       color = "Hazard")

ggsave(paste0("figs/Hazard_hos_resitalized_age_groups_n",
              format(Sys.Date(), "%F"), ".png"), 
       width = 7, height = 5)
