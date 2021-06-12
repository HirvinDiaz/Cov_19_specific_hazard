n_sim <- 1000

days_hosp  =  rlnorm(n_sim, 
                     meanlog = 2.4321, 
                     sdlog = 0.5153)

days_icu  =  rlnorm(n_sim, 
                    meanlog = 2.3105, 
                    sdlog = 0.8008)

df_hosp_time <- as.data.frame(cbind(days_hosp, days_icu))

df_hosp_time_long <- gather(data = df_hosp_time, key = "Cohort", value = "Days")

df_hosp_time_long$Cohort[df_hosp_time_long$Cohort == "days_hosp"] <- "Hospitalized, Not Intubated"
df_hosp_time_long$Cohort[df_hosp_time_long$Cohort == "days_icu"] <- "Hospitalized, Intubated"

hr_Remd      = rlnorm(n_sim,
                      meanlog = log(0.73), 
                      sd = (log(1.03)-log(0.52))/(2*1.96))
hist(hr_Remd)
hr_RemdBa_vs_Remd   = rlnorm(n_sim, 
                             meanlog = log(0.65), 
                             sd = (log(1.09)-log(0.39))/(2*1.96))
hr_Dexa    = rlnorm(n_sim,
                    meanlog = log(0.64), 
                    sd = (log(0.81)-log(0.51))/(2*1.96))

df_effects_params <- as.data.frame(cbind(hr_Remd, hr_RemdBa_vs_Remd, hr_Dexa))

df_effects_params_long <- gather(data = df_effects_params, key = "Cohort", value = "HR")

df_effects_params_long$Cohort[df_effects_params_long$Cohort == "hr_Remd"] <- "Remdesivir"
df_effects_params_long$Cohort[df_effects_params_long$Cohort == "hr_RemdBa_vs_Remd"] <- "Remdesivir & Baricitinib"
df_effects_params_long$Cohort[df_effects_params_long$Cohort == "hr_Dexa"] <- "Dexamethasone"

# Hosp Costs
c_hosp     = rgamma(n_sim, shape = 8, scale = 1159)

c_resp     = rgamma(n_sim, shape = 8, scale = 5518.875)

df_costs_params <- as.data.frame(cbind(c_hosp, c_resp))

df_costs_params_long <- gather(data = df_costs_params, key = "Cohort", value = "Costs")

df_costs_params_long$Cohort[df_costs_params_long$Cohort == "c_hosp"] <- "Hospitalized, Not Intubated"
df_costs_params_long$Cohort[df_costs_params_long$Cohort == "c_resp"] <- "Hospitalized, Intubated"

hist(c_hosp)

t_hos_parms <- lnorm_params(m = 13.89, v = (13.17)^2)

t_hosp   = rlnorm(1000, meanlog = t_hos_parms$mu, sdlog = t_hos_parms$sigma)

hist(t_hosp, 
     main = "Time of hospitalization",
     xlab = "Days",
     ylab = " ",
     col = "#eb961e")

