hr_Remd      = rlnorm(n_sim,
                      meanlog = log(0.73), 
                      sd = (log(1.03)-log(0.52))/(2*1.96))
hist(hr_Remd, 
     main = "Remdesivir Effect",
     xlab = "Hazard Ratio", 
     ylab = " ",
     col = "#48acd4")

hr_RemdBa_vs_Remd   = rlnorm(n_sim, 
                             meanlog = log(0.65), 
                             sd = (log(1.09)-log(0.39))/(2*1.96))

hist(hr_RemdBa_vs_Remd, 
     main = "Baricitinib Effect",
     xlab = "Hazard Ratio",
     ylab = " ",
     col = "#6abd6b")

c_hosp     = rgamma(n_sim, shape = 8, scale = 1159)

hist(c_hosp, 
     main = "Daily costs hospitalization, not intubation",
     xlab = "Mexican Pesos",
     ylab = " ",
     col = "#e89f72")

hr_Dexa    = rlnorm(n_sim,
                    meanlog = log(0.64), 
                    sd = (log(0.81)-log(0.51))/(2*1.96))

hist(hr_Dexa, 
     main = "Dexamethasone Effect",
     xlab = "Hazard Ratio", 
     ylab = " ",
     col = "#34ed24")

c_resp     = rgamma(n_sim, shape = 8, scale = 5518.875)

hist(c_resp, 
     main = "Daily costs hospitalization with intubation",
     xlab = "Mexican Pesos",
     ylab = " ",
     col = "#eb961e")
