set.seed(02021989)
x <- exp(rnorm(1000, log(0.73), sd = (log(1.03)-log(0.52))/(2*1.96)))
mean(x)
quantile(x, probs = c(0.025, 0.975))

# Todos los parametros deben de tener sus distribuciones
# Para sacar los costos, se suma y se restan 20% 
# Utilidades Años de vida, utilidades, hospitalizados.

y <- rlnorm(1000, log(0.73), sd = (log(1.03)-log(0.52))/(2*1.96))
y <- rlnorm(1000, meanlog = log(0.73), sd = (log(1.03)-log(0.52))/(2*1.96))
mean(y)

1.03 - 0.73
0.73 - 0.3

c_remd <- rgamma(1000, shape = 106.743, scale = 100)  
hist(c_remd)
View(c_remd)

c_Bari <- rgamma(1000, shape = 36.72, scale = 100)
hist(c_Bari)
View(c_Bari)

c_hosp     <- rgamma(1000, shape = 8, scale = 1159)
hist(c_hosp)
View(c_hosp)

gamma_params(mu = 9272, sigma = 3278.14704)

gamma_params(mu = 44151, sigma = 15609.7357)

c_hosp     <- rgamma(1000, shape = 8, scale = 5518.875)
hist(c_hosp)
View(c_hosp)



