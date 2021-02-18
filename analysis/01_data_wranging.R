#### Data wrangling. Author: Hirvin Azael Diaz Zepeda ####
library(readr)
library(dplyr)

# Import data base 14 - 11 - 2021
Cov_09_02 <- read_csv("data-raw/210209COVID19MEXICO.csv")

Cov <- Cov_09_02 %>% 
  select(2,6,10,11,12,13,14,16,33,36,40)

Cov <- Cov %>% 
  filter(CLASIFICACION_FINAL == 1 | CLASIFICACION_FINAL == 2 | CLASIFICACION_FINAL == 3) %>% 
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
  filter(date_admission <= max(date_admission) - 7) %>% 
  mutate(state = ifelse(is.na(date_death), 0, 1))

save(Cov, file = "data/cov_09_02.Rdata")


