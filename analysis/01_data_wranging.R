#### Data wrangling. Author: Hirvin Azael Diaz Zepeda ####
library(readr)
library(dplyr)

# Import data base 14 - 11 - 2021
Cov_14_01 <- read_csv("data-raw/210114COVID19MEXICO.csv")

Cov <- Cov_14_01 %>% 
  select(2,6,10,11,12,13,14,16,33,40)

Cov <- Cov %>% 
  filter(RESULTADO_LAB == 1) %>% 
  rename(ID = ID_REGISTRO,
         sex = SEXO,
         type = TIPO_PACIENTE,
         date_admission = FECHA_INGRESO,
         date_symptoms = FECHA_SINTOMAS,
         date_death = FECHA_DEF,
         intubated = INTUBADO,
         age = EDAD,
         icu = UCI,
         result = RESULTADO_LAB,
  ) %>% 
  filter(date_admission <= max(date_admission) - 10) %>% 
  mutate(state = ifelse(is.na(date_death), 0, 1))

save(Cov, file = "data/cov_14_01.Rdata")
