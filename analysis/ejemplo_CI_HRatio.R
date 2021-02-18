x <- exp(rnorm(1000, log(0.73), sd = (log(1.03)-log(0.52))/(2*1.96)))
mean(x)
quantile(x, probs = c(0.025, 0.975))

# Todos los parametros deben de tener sus distribuciones
# Para sacar los costos, se suma y se restan 20% 
# Utilidades Años de vida, utilidades, hospitalizados.
# 