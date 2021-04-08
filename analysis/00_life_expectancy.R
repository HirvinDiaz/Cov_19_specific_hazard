library(dplyr)
library(readxl)
library(tidyverse)

# Data with life expectancy
age_life_expectancy <- read_excel("data/age_life_expectancy.xlsx")

Life_expectancy <- age_life_expectancy %>% 
  mutate(Male = male - age, Female = female -  age) %>% 
  select(c("age", "Male", "Female"))

Life_expectancy <- Life_expectancy %>% 
  rename(male = Male,
         female = Female)

Life_expectancy <- gather(data = Life_expectancy, 
                          key = "sex", 
                          value = "Years", 
                          -age)

life_expected_dis <- read_excel("data/life_expected_dis.xlsx")
View(life_expected_dis)

Life_expectancy <- Life_expectancy %>% 
  left_join(life_expected_dis, by = "Years")

save(Life_expectancy, file = "data/Life_expectancy.Rdata")
