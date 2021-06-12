library(dplyr)

df_hazards_hos_fe <- df_hazards_hos_res %>% 
  filter(sex == "Female")

mean(df_hazards_hos_fe$Hazard)*100000


df_hazards_hos_me <- df_hazards_hos_res %>% 
  filter(sex == "Male")

mean(df_hazards_hos_me$Hazard)*100000
