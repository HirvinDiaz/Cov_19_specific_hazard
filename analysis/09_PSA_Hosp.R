#### 09 Probabilistic Sensitibity Analysis  ####

# Based in the code developed by the Decision Analysis in R for Technologies in 
# Health (DARTH) workgroup:

n_sim <- 1000
# Function to generate PSA parameters ()
psa_params <- function(n_sim){
  set.seed(02021989) # set a seed to be able to reproduce the same results
  df_parameters <- data.frame(
    # Hazard ratios 
    Remd_effect      = rlnorm(n_sim,
                          meanlog = log(0.73), 
                          sd = (log(1.03)-log(0.52))/(2*1.96)), 
    # Hazard Ratio of Remdesivir
    Remd_Ba_effect   = rlnorm(n_sim, 
                          meanlog = log(0.65), 
                          sd = (log(1.09)-log(0.39))/(2*1.96)),
    # Hazard Ratio of Remdesivir + Ba 
    # Costs
    c_hosp     = rgamma(n_sim, shape = 81.4163, scale = 100),  # cost of remaining one cycle hospitalized 
    c_Remd     = rgamma(n_sim, shape = 106.743, scale = 100),
    # cost of Remdesivir Treatment for one Cycle
    c_Bari     = rgamma(n_sim, shape = 36.72, scale = 100),     
    # cost of Baricitinib Treatment for one Cycle
    c_dead     = rep(0, n_sim),        # cost of remaining one cycle Dead
    # Utilities
    l_sick     = rep(1, n_sim), # utility when Sick 
    l_dead     = rep(0, n_sim), # utility when Dead                              
    d_e        = rep(0.0082, n_sim)   # discount factor for costs
  )
  return(df_parameters)
}

# Try it
df_par <- psa_params(1000) 

l_params_all <- list( #El elemento que estamos creando es una lista
  Remd_effect    = 0.73,  # probability to die when healthy
  Remd_Ba_effect = 0.65,   # probability to become sick when healthy
  c_hosp         = 8100,    # probability to become healthy when sick
  c_Remd         = 10674.3,  # probability to become sicker when sick
  c_Bari         = 3672.35,      # hazard ratio of death in sick vs healthy
  c_dead         = 0,     # hazard ratio of death in sicker vs healthy
  l_sick         = 1,   # cost of remaining one cycle in the healthy state
  l_dead         = 0,   # cost of remaining one cycle in the sick state
  d_e            = 0.0082    # discount factor for costs
)

x <- as.list(l_params_all)

i <- 1
df <- updt_prm_list(l_params_all,df_par[i,])

# Run Markov model on each parameter set of PSA input dataset
for(i in 1:n_sim){
  l_psa_input <- update_param_list(l_params_all, df_psa_input[i,])
  df_out_psa  <- calculate_ce_out(l_psa_input)
  df_c[i, ] <- df_out_psa$Cost
  df_e[i, ] <- df_out_psa$Effect
  # Display simulation progress
  if(i/(n_sim/10) == round(i/(n_sim/10), 0)) { # display progress every 10%
    cat('\r', paste(i/n_sim * 100, "% done", sep = " "))
  }
}
