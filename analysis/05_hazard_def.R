#### Final project script. Author: Hirvin Azael Diaz Zepeda ####

####  Clean Work - space  ####
rm(list = ls()) # to clean the work - space

####  Load Packages  ####
library(dplyr)
library(ggplot2)
library(utils)
library(scales)
library(tidyr)
library(data.table)
library(reshape2)
library(survival)
library(survPen)
library(relsurv)
library(readr)
library(readxl)
library(lubridate)

######  Load data bases  ######

# Data base of suspected people with COVID-19 published by Mexico´s Ministry of 
# Health 

load("data/cov_09_02.Rdata") 

#### Load functions ####
source("analysis/04_r_functions.R")

#### Load ratetable ####
load("data/rate_mx.Rdata")

#### Compute Hazards for hospitalized population ####
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


ggplot(filter(df_hazards_hos_res,
              Type == "COVID-19"), 
       aes(x = time, 
           y = log(Hazard),
           color = Type)) +
  geom_line(size = 0.8)+
  geom_line(data = filter(df_hazards_hos_res,
                          Type == "Observed"), 
            aes(x = time, 
                y = log(Hazard)),
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

ggsave(paste0("figs/Hazard_hos_resitalized_age_groups_log",
              format(Sys.Date(), "%F"), ".png"), 
       width = 7, height = 5)

#### By Month ####

Covid_rr <- Cov %>% 
  filter(intubated == 1 | type == 2) %>% 
  filter(age >= 45) %>% 
  select(c("sex", "age", "date_admission", "date_death", "date_symptoms", "state"))%>% 
  mutate(sex = ifelse(sex == 1, 2, 1))%>% # recode 1 is male and 2 female
  mutate(time = ifelse(is.na(date_death),x - date_symptoms, 
                       date_death - date_symptoms)) %>% # time of observation is 
  mutate(diag = date_symptoms) %>%                 # date of death minus date of 
  mutate(diag = as.numeric(diag, origin = "1960-01-01")) %>%
  filter(time >= 1) %>% # eliminate negative times
  mutate(original_age = age) %>% 
  mutate(age = original_age*365.25) %>%  # age needs to be in days for relsurv
  mutate(age_range = ifelse(original_age >= 45 & original_age <= 54, "45 - 54", 
                            ifelse(original_age >= 55 & original_age <= 64, "55 - 64",
                                   ifelse(original_age >= 65 & original_age <= 69, 
                                          "65 - 69", "70 +")))) %>% 
  mutate(month = ifelse(date_symptoms > "2020-01-01" & date_symptoms < "2020-04-01", "March", 
                        ifelse(date_symptoms >= "2020-04-01" & date_symptoms <= "2020-04-30", "April", 
                               ifelse(date_symptoms >= "2020-05-01" & date_symptoms <= "2020-05-31", "May",
                                      ifelse(date_symptoms >= "2020-06-01" & date_symptoms <= "2020-06-30", "June",
                                             ifelse(date_symptoms >= "2020-07-01" & date_symptoms <= "2020-07-31", "July",
                                                    ifelse(date_symptoms >= "2020-08-01" & date_symptoms <= "2020-08-31", "August",
                                                           ifelse(date_symptoms >= "2020-09-01" & date_symptoms <= "2020-09-30", "September",
                                                                  ifelse(date_symptoms >= "2020-10-01" & date_symptoms <= "2020-10-31", "October",
                                                                         ifelse(date_symptoms >= "2020-11-01" & date_symptoms <= "2020-11-30", "November",
                                                                                ifelse(date_symptoms >= "2020-12-01" & date_symptoms <= "2020-12-31", "December", "January 21")))))))))))



#### Data_frame by month
Covid_rr_March_20 <- Covid_rr %>% 
  filter(month == "March")

Covid_rr_April_20 <- Covid_rr %>% 
  filter(month == "April")

Covid_rr_May_20 <- Covid_rr %>% 
  filter(month == "May")

Covid_rr_June_20 <- Covid_rr %>% 
  filter(month == "June")

Covid_rr_July_20 <- Covid_rr %>% 
  filter(month == "July")

Covid_rr_August_20 <- Covid_rr %>% 
  filter(month == "August")

Covid_rr_September_20 <- Covid_rr %>% 
  filter(month == "September")

Covid_rr_October_20 <- Covid_rr %>% 
  filter(month == "October")

Covid_rr_November_20 <- Covid_rr %>% 
  filter(month == "November")

Covid_rr_December_20 <- Covid_rr %>% 
  filter(month == "December")

Covid_rr_January_21 <- Covid_rr %>% 
  filter(month == "January 21")

#### KM.Surv functions ####

KM_surv_March_20 <- rs.surv.exp(Surv(time, state) ~ 1,
                               data = Covid_rr_March_20,
                               ratetable = rate_exp_mx, 
                               type="kaplan-meier",
                               conf.type="log",
                               conf.int=0.95,
                               rmap = list(age = age, 
                                           sex = sex, 
                                           year = diag))

KM_surv_April_20 <- rs.surv.exp(Surv(time, state) ~ 1,
                                data = Covid_rr_April_20,
                                ratetable = rate_exp_mx, 
                                type="kaplan-meier",
                                conf.type="log",
                                conf.int=0.95,
                                rmap = list(age = age, 
                                            sex = sex, 
                                            year = diag))

KM_surv_May_20 <- rs.surv.exp(Surv(time, state) ~ 1,
                                data = Covid_rr_May_20,
                                ratetable = rate_exp_mx, 
                                type="kaplan-meier",
                                conf.type="log",
                                conf.int=0.95,
                                rmap = list(age = age, 
                                            sex = sex, 
                                            year = diag))

KM_surv_June_20 <- rs.surv.exp(Surv(time, state) ~ 1,
                                data = Covid_rr_June_20,
                                ratetable = rate_exp_mx, 
                                type="kaplan-meier",
                                conf.type="log",
                                conf.int=0.95,
                                rmap = list(age = age, 
                                            sex = sex, 
                                            year = diag))

KM_surv_July_20 <- rs.surv.exp(Surv(time, state) ~ 1,
                                data = Covid_rr_July_20,
                                ratetable = rate_exp_mx, 
                                type="kaplan-meier",
                                conf.type="log",
                                conf.int=0.95,
                                rmap = list(age = age, 
                                            sex = sex, 
                                            year = diag))

KM_surv_August_20 <- rs.surv.exp(Surv(time, state) ~ 1,
                                data = Covid_rr_August_20,
                                ratetable = rate_exp_mx, 
                                type="kaplan-meier",
                                conf.type="log",
                                conf.int=0.95,
                                rmap = list(age = age, 
                                            sex = sex, 
                                            year = diag))

KM_surv_September_20 <- rs.surv.exp(Surv(time, state) ~ 1,
                                data = Covid_rr_September_20,
                                ratetable = rate_exp_mx, 
                                type="kaplan-meier",
                                conf.type="log",
                                conf.int=0.95,
                                rmap = list(age = age, 
                                            sex = sex, 
                                            year = diag))

KM_surv_October_20 <- rs.surv.exp(Surv(time, state) ~ 1,
                                data = Covid_rr_October_20,
                                ratetable = rate_exp_mx, 
                                type="kaplan-meier",
                                conf.type="log",
                                conf.int=0.95,
                                rmap = list(age = age, 
                                            sex = sex, 
                                            year = diag))

KM_surv_November_20 <- rs.surv.exp(Surv(time, state) ~ 1,
                                data = Covid_rr_November_20,
                                ratetable = rate_exp_mx, 
                                type="kaplan-meier",
                                conf.type="log",
                                conf.int=0.95,
                                rmap = list(age = age, 
                                            sex = sex, 
                                            year = diag))

KM_surv_December_20 <- rs.surv.exp(Surv(time, state) ~ 1,
                                data = Covid_rr_December_20,
                                ratetable = rate_exp_mx, 
                                type="kaplan-meier",
                                conf.type="log",
                                conf.int=0.95,
                                rmap = list(age = age, 
                                            sex = sex, 
                                            year = diag))

plot(KM_surv_December_20)
# Create Data frame by month
March_20  <- head(KM_surv_March_20$hazard_exc, 60)
April_20  <- head(KM_surv_April_20$hazard_exc, 60)
May_20   <- head(KM_surv_May_20$hazard_exc, 60)
June_20 <- head(KM_surv_June_20$hazard_exc, 60)
July_20 <- head(KM_surv_July_20$hazard_exc, 60)
August_20 <- head(KM_surv_August_20$hazard_exc, 60)
September_20 <- head(KM_surv_September_20$hazard_exc, 60)
October_20 <- head(KM_surv_October_20$hazard_exc, 60)
November_20 <- head(KM_surv_November_20$hazard_exc, 60)
December_20 <- head(KM_surv_December_20$hazard_exc, 60)

time <- seq(1,60, by = 1)
March_df <- as.data.frame(cbind(time, March_20, "March"))
March_df <- March_df %>% 
  rename(Month = March_20)
April_df <- as.data.frame(cbind(time, April_20, "April")) 
April_df <- April_df %>% 
  rename(Month = April_20)
May_df <- as.data.frame(cbind(time, May_20, "May"))
May_df <- May_df %>% 
  rename(Month = May_20)
June_df <- as.data.frame(cbind(time, June_20, "June"))
June_df <- June_df %>% 
  rename(Month = June_20)
July_df <- as.data.frame(cbind(time, July_20, "July")) 
July_df <- July_df %>% 
  rename(Month = July_20)
August_df <- as.data.frame(cbind(time, August_20, "August"))
August_df <- August_df %>% 
  rename(Month = August_20)
September_df <- as.data.frame(cbind(time, September_20, "September")) 
September_df <- September_df %>% 
  rename(Month = September_20)
October_df <- as.data.frame(cbind(time, October_20, "October")) 
October_df <- October_df %>% 
  rename(Month = October_20)
November_df <- as.data.frame(cbind(time, November_20, "November"))
November_df <- November_df %>% 
  rename(Month = November_20)
December_df <- as.data.frame(cbind(time, December_20, "December")) 
December_df <- December_df %>% 
  rename(Month = December_20)

Months_DF <- rbind(March_df, April_df, May_df, June_df, July_df, August_df, September_df,
                   October_df, November_df, December_df)
Months_DF$Month <- as.numeric(Months_DF$Month)
Months_DF$time <- as.numeric(Months_DF$time)

unique(Months_DF$V3)
Months_DF$label <- factor(Months_DF$V3, 
                          levels = c("March", "April", "May", "June", 
                                                   "July", "August", "September", 
                                     "October","November", "December"))
save(Months_DF, file = "Months_DF.Rdata")

ggplot(Months_DF, aes(x = time, 
                      y = Month, 
                      color = V3))+
  geom_line(size = 1)+
  facet_wrap(~label, ncol = 2)+
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
                                        colour = "#6e6866"),
        legend.position = "none")+
  scale_x_continuous(breaks = number_ticks(4))+
  scale_y_continuous(breaks = number_ticks(2))+
  # axis.text.x = element_text(angle = 90, hjust = 0))+
  scale_color_manual(values=c ("#120d0c", "#e32402", "#eaacf2", "#2d0da1", 
                               "#d19f08", "#156e1c", "#d615c9", "#1fd10f", 
                               "#37edf0", "#753731")) +
  labs(title = "Covid-19 specific hazard by month",
       x = "Days",
       y = " ",
       label = "")

ggsave(paste0("figs/sphazard by month",
              format(Sys.Date(), "%F"), ".png"), 
       width = 7, height = 5)


#### Compute Hazards just for hospitalized population ####
### Hospitalized and no intubated ###

#### Covid_rr_hosp DB ####
x <- as.Date(max(Cov$date_admission), format = "%Y-%m-%d")

Covid_rr_hosp <- Cov %>% 
  filter(type == 2 & intubated != 1) %>% 
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
Covid_rr_hosp_male_45 <- Covid_rr_hosp %>% 
  filter(sex == 1 & age_range == "45 - 54")

dim(Covid_rr_hosp_male_45)

Covid_rr_hosp_male_55 <- Covid_rr_hosp %>% 
  filter(sex == 1 & age_range == "55 - 64")

dim(Covid_rr_hosp_male_55)

Covid_rr_hosp_male_70 <- Covid_rr_hosp %>% 
  filter(sex == 1 & age_range == "70 +")

dim(Covid_rr_hosp_male_70)

Covid_rr_hosp_male_65 <- Covid_rr_hosp %>% 
  filter(sex == 1 & age_range == "65 - 69")

dim(Covid_rr_hosp_male_65)

### Mujeres ####
Covid_rr_hosp_female_45 <- Covid_rr_hosp %>% 
  filter(sex == 2 & age_range == "45 - 54")

dim(Covid_rr_hosp_female_45)

Covid_rr_hosp_female_55 <- Covid_rr_hosp %>% 
  filter(sex == 2 & age_range == "55 - 64")

dim(Covid_rr_hosp_female_55)

Covid_rr_hosp_female_70 <- Covid_rr_hosp %>% 
  filter(sex == 2 & age_range == "70 +")

dim(Covid_rr_hosp_female_70)

Covid_rr_hosp_female_65 <- Covid_rr_hosp %>% 
  filter(sex == 2 & age_range == "65 - 69")

#### Apply rs.surv function
KM_surv_sex_male_45 <- rs.surv.exp(Surv(time, state) ~ 1,
                                   data = Covid_rr_hosp_male_45,
                                   ratetable = rate_exp_mx, 
                                   type="kaplan-meier",
                                   conf.type="log",
                                   conf.int=0.95,
                                   rmap = list(age = age, sex = sex, year = diag))

KM_surv_sex_male_55 <- rs.surv.exp(Surv(time, state) ~ 1,
                                   data = Covid_rr_hosp_male_55,
                                   ratetable = rate_exp_mx, 
                                   type="kaplan-meier",
                                   conf.type="log",
                                   conf.int=0.95,
                                   rmap = list(age = age, sex = sex, year = diag))

KM_surv_sex_male_65 <- rs.surv.exp(Surv(time, state) ~ 1,
                                   data = Covid_rr_hosp_male_65,
                                   ratetable = rate_exp_mx, 
                                   type="kaplan-meier",
                                   conf.type="log",
                                   conf.int=0.95,
                                   rmap = list(age = age, sex = sex, year = diag))

KM_surv_sex_male_70 <- rs.surv.exp(Surv(time, state) ~ 1,
                                   data = Covid_rr_hosp_male_70,
                                   ratetable = rate_exp_mx, 
                                   type="kaplan-meier",
                                   conf.type="log",
                                   conf.int=0.95,
                                   rmap = list(age = age, sex = sex, year = diag))

#### Female ####
KM_surv_sex_female_45 <- rs.surv.exp(Surv(time, state) ~ 1,
                                     data = Covid_rr_hosp_female_45,
                                     ratetable = rate_exp_mx, 
                                     type="kaplan-meier",
                                     conf.type="log",
                                     conf.int=0.95,
                                     rmap = list(age = age, sex = sex, year = diag))

KM_surv_sex_female_55 <- rs.surv.exp(Surv(time, state) ~ 1,
                                     data = Covid_rr_hosp_female_55,
                                     ratetable = rate_exp_mx, 
                                     type="kaplan-meier",
                                     conf.type="log",
                                     conf.int=0.95,
                                     rmap = list(age = age, sex = sex, year = diag))

KM_surv_sex_female_65 <- rs.surv.exp(Surv(time, state) ~ 1,
                                     data = Covid_rr_hosp_female_65,
                                     ratetable = rate_exp_mx, 
                                     type="kaplan-meier",
                                     conf.type="log",
                                     conf.int=0.95,
                                     rmap = list(age = age, sex = sex, year = diag))

KM_surv_sex_female_70 <- rs.surv.exp(Surv(time, state) ~ 1,
                                     data = Covid_rr_hosp_female_70,
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

df_hazards_hosp <- rbind(df_hazard_long_female, df_hazard_long_male)

df_hazards_hosp <- df_hazards_hosp %>% 
  mutate(state = "Hosp")


#### Compute Hazards just for ICU population ####
### Intubated patients ###

#### Covid_rr_ICU DB ####
x <- as.Date(max(Cov$date_admission), format = "%Y-%m-%d")

Covid_rr_ICU <- Cov %>% 
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
Covid_rr_ICU_male_45 <- Covid_rr_ICU %>% 
  filter(sex == 1 & age_range == "45 - 54")

dim(Covid_rr_ICU_male_45)

Covid_rr_ICU_male_55 <- Covid_rr_ICU %>% 
  filter(sex == 1 & age_range == "55 - 64")

dim(Covid_rr_ICU_male_55)

Covid_rr_ICU_male_70 <- Covid_rr_ICU %>% 
  filter(sex == 1 & age_range == "70 +")

dim(Covid_rr_ICU_male_70)

Covid_rr_ICU_male_65 <- Covid_rr_ICU %>% 
  filter(sex == 1 & age_range == "65 - 69")

dim(Covid_rr_ICU_male_65)

### Mujeres ####
Covid_rr_ICU_female_45 <- Covid_rr_ICU %>% 
  filter(sex == 2 & age_range == "45 - 54")

dim(Covid_rr_ICU_female_45)

Covid_rr_ICU_female_55 <- Covid_rr_ICU %>% 
  filter(sex == 2 & age_range == "55 - 64")

dim(Covid_rr_ICU_female_55)

Covid_rr_ICU_female_70 <- Covid_rr_ICU %>% 
  filter(sex == 2 & age_range == "70 +")

dim(Covid_rr_ICU_female_70)

Covid_rr_ICU_female_65 <- Covid_rr_ICU %>% 
  filter(sex == 2 & age_range == "65 - 69")

#### Apply rs.surv function
KM_surv_sex_male_45 <- rs.surv.exp(Surv(time, state) ~ 1,
                                   data = Covid_rr_ICU_male_45,
                                   ratetable = rate_exp_mx, 
                                   type="kaplan-meier",
                                   conf.type="log",
                                   conf.int=0.95,
                                   rmap = list(age = age, sex = sex, year = diag))

KM_surv_sex_male_55 <- rs.surv.exp(Surv(time, state) ~ 1,
                                   data = Covid_rr_ICU_male_55,
                                   ratetable = rate_exp_mx, 
                                   type="kaplan-meier",
                                   conf.type="log",
                                   conf.int=0.95,
                                   rmap = list(age = age, sex = sex, year = diag))

KM_surv_sex_male_65 <- rs.surv.exp(Surv(time, state) ~ 1,
                                   data = Covid_rr_ICU_male_65,
                                   ratetable = rate_exp_mx, 
                                   type="kaplan-meier",
                                   conf.type="log",
                                   conf.int=0.95,
                                   rmap = list(age = age, sex = sex, year = diag))

KM_surv_sex_male_70 <- rs.surv.exp(Surv(time, state) ~ 1,
                                   data = Covid_rr_ICU_male_70,
                                   ratetable = rate_exp_mx, 
                                   type="kaplan-meier",
                                   conf.type="log",
                                   conf.int=0.95,
                                   rmap = list(age = age, sex = sex, year = diag))

#### Female ####
KM_surv_sex_female_45 <- rs.surv.exp(Surv(time, state) ~ 1,
                                     data = Covid_rr_ICU_female_45,
                                     ratetable = rate_exp_mx, 
                                     type="kaplan-meier",
                                     conf.type="log",
                                     conf.int=0.95,
                                     rmap = list(age = age, sex = sex, year = diag))

KM_surv_sex_female_55 <- rs.surv.exp(Surv(time, state) ~ 1,
                                     data = Covid_rr_ICU_female_55,
                                     ratetable = rate_exp_mx, 
                                     type="kaplan-meier",
                                     conf.type="log",
                                     conf.int=0.95,
                                     rmap = list(age = age, sex = sex, year = diag))

KM_surv_sex_female_65 <- rs.surv.exp(Surv(time, state) ~ 1,
                                     data = Covid_rr_ICU_female_65,
                                     ratetable = rate_exp_mx, 
                                     type="kaplan-meier",
                                     conf.type="log",
                                     conf.int=0.95,
                                     rmap = list(age = age, sex = sex, year = diag))

KM_surv_sex_female_70 <- rs.surv.exp(Surv(time, state) ~ 1,
                                     data = Covid_rr_ICU_female_70,
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

df_hazards_ICU <- rbind(df_hazard_long_female, df_hazard_long_male)

df_hazards_ICU <- df_hazards_ICU %>% 
  mutate(state = "ICU")

df_hazards_ICU_hosp <- rbind(df_hazards_hosp,df_hazards_ICU)

save(df_hazards_ICU_hosp, file = "data/df_hazards_ICU_hosp.Rdata")
