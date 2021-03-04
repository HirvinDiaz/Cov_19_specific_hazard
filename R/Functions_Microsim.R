#### Functions for Cost-Effectiveness Analysis ####

source("R/samplev.R")

updt_prm_list <- function(l_params_all, params_updated){
  if (typeof(params_updated)!="list"){
    params_updated <- split(unname(params_updated),names(params_updated)) 
  }
  l_params_all <- modifyList(l_params_all, params_updated) #update the values
  return(l_params_all)
}

Probs <- function(v_M_t, df_X, t, n_states, n_i, Trt = FALSE) { # t <- 1
  m_p_t           <- matrix(0, nrow = n_states, ncol = n_i)  
  # give the state names to the rows
  rownames(m_p_t) <-  v_names_states
  # Lookup baseline probability of dying from Covid-19 or other causes  
  if (Trt == "Remd"){
    p_die_CoV_all <- dt_p_CoV[.(df_X$group, df_X$time + t, df_X$sex), p_dCoV_Remd] 
  } else if (Trt == "Remd_Ba"){
    p_die_CoV_all <- dt_p_CoV[.(df_X$group, df_X$time + t, df_X$sex), p_dCoV_Remd_Ba] 
  } else if (Trt == "None" | Trt == FALSE){
    p_die_CoV_all <- dt_p_CoV[.(df_X$group, df_X$time + t, df_X$sex), p_dCoV] 
  } else {
    message("Choose a treatment within the options")
  }
  p_die_Pop_all <- dt_p_CoV[.(df_X$group, df_X$time + t, df_X$sex), p_dPop] 
  
  p_die_CoV     <- p_die_CoV_all[v_M_t == "Cov19+"]  
  p_die_Pop     <- p_die_Pop_all[v_M_t == "Cov19+"] 
  # transition probabilities when healthy
  m_p_t[, v_M_t == "Cov19+"] <- rbind(1 - (p_die_CoV + p_die_Pop), p_die_CoV, p_die_Pop)    
  # transition probabilities when sick 
  m_p_t[, v_M_t == "CoV19_Dead"] <- rbind(0, 1 ,0)
  # transition probabilities when sicker 
  m_p_t[, v_M_t == "O_Causes_Dead"] <- rbind(0, 0, 1)
  return(t(m_p_t))
}     

Costs <- function (l_params, v_M_t, Trt = FALSE) {
  with(as.list(l_params),{
  # v_M_t: current health state
  c_t <- c()
  if (Trt == "Remd"){
    c_t[v_M_t == "Cov19+"]  <- c_hosp + c_remd
  } else if (Trt == "Remd_Ba"){
    c_t[v_M_t == "Cov19+"]  <- c_hosp + c_remd + c_bari
  } else if (Trt == "None" | Trt == FALSE){
    c_t[v_M_t == "Cov19+"]  <- c_hosp 
  } else {
    message("Choose a treatment within the options")
  }
  # costs accrued by being Sick this cycle
  c_t[v_M_t == "CoV19_Dead"] <- c_dead  # costs accrued by being sick this cycle
  c_t[v_M_t == "O_Causes_Dead"] <- c_dead  # costs accrued by being sicker this cycle
  return(c_t)  # return costs accrued this cycle
  }
  )
}

Effs <- function (l_params, v_M_t) {
  with(as.list(l_params),{
  # v_M_t: current health state
  q_t <- c() 
  q_t[v_M_t == "Cov19+"] <- l_sick   
  q_t[v_M_t == "CoV19_Dead"] <- l_dead     
  q_t[v_M_t == "O_Causes_Dead"]  <- l_dead     
  return(q_t)  # return the LDs accrued this cycle
  }
  )
}

CEA_microsim <- function(l_params, df_X, df_h, Trt = FALSE) { 
  with(as.list(l_params), {
       set.seed(02021989)
       
       d_h_HD_hosp <- df_h # Data frame with specific hazards
       Remd_effect <- Remd_effect
       Remd_Ba_effect <- Remd_Ba_effect
       
       # Transform hazards to transition probabilities
       d_h_HD_hosp <- d_h_HD_hosp %>% 
         mutate(haz_Remd = Hazard.Covid_19*Remd_effect) %>% 
         mutate(haz_Remd_Ba = haz_Remd*Remd_Ba_effect)
       
       d_p_HD_hosp <- d_h_HD_hosp %>% 
         mutate(p_dCoV = 1 - exp(-Hazard.Covid_19)) %>% 
         mutate(p_dPop = 1 - exp(- Hazard.Population)) %>% 
         mutate(p_dCoV_Remd = 1 - exp(- haz_Remd)) %>% 
         mutate(p_dCoV_Remd_Ba = 1 - exp(- haz_Remd_Ba)) %>% 
         arrange(group, time)%>% 
         select(c("time", "group", "sex", "p_dCoV", "p_dPop", "p_dCoV_Remd", 
                  "p_dCoV_Remd_Ba"))
       
       # Convert data frame to a data table for efficiency
       dt_p_CoV <- data.table(d_p_HD_hosp)
       
       # set the data table to be indexed by age, day, sex and ICU
       setkey(dt_p_CoV, group, time, sex)
       
       n_i              <- length(df_X$sex)
       n_t              <- 50 # time horizon, 50 days
       v_names_states   <- c("Cov19+", "CoV19_Dead", "O_Causes_Dead") # state names
       n_states         <- length(v_names_states)   # number of states
       
       # vector with discount weights for costs
       d_c   <- 0.0082 # daily discount rates 
       v_dwc <- 1 / (1 + d_c) ^ (0:n_t) 
       v_dwe <- 1 / (1) ^ (0:n_t) 
       
       Probs <- function(v_M_t, df_X, t, Trt = FALSE) { # t <- 1
         m_p_t           <- matrix(0, nrow = n_states, ncol = n_i)  
         # give the state names to the rows
         rownames(m_p_t) <-  v_names_states
         # Lookup baseline probability of dying from Covid-19 or other causes  
         if (Trt == "Remd"){
           p_die_CoV_all <- dt_p_CoV[.(df_X$group, df_X$time + t, df_X$sex), p_dCoV_Remd] 
         } else if (Trt == "Remd_Ba"){
           p_die_CoV_all <- dt_p_CoV[.(df_X$group, df_X$time + t, df_X$sex), p_dCoV_Remd_Ba] 
         } else if (Trt == "None" | Trt == FALSE){
           p_die_CoV_all <- dt_p_CoV[.(df_X$group, df_X$time + t, df_X$sex), p_dCoV] 
         } else {
           message("Choose a treatment within the options")
         }
         p_die_Pop_all <- dt_p_CoV[.(df_X$group, df_X$time + t, df_X$sex), p_dPop] 
         
         p_die_CoV     <- p_die_CoV_all[v_M_t == "Cov19+"]  
         p_die_Pop     <- p_die_Pop_all[v_M_t == "Cov19+"] 
         # transition probabilities when healthy
         m_p_t[, v_M_t == "Cov19+"] <- rbind(1 - (p_die_CoV + p_die_Pop), p_die_CoV, p_die_Pop)    
         # transition probabilities when sick 
         m_p_t[, v_M_t == "CoV19_Dead"] <- rbind(0, 1 ,0)
         # transition probabilities when sicker 
         m_p_t[, v_M_t == "O_Causes_Dead"] <- rbind(0, 0, 1)
         return(t(m_p_t))
       }     
       
       m_M <- m_C <- m_E <-  matrix(NA, nrow = n_i, ncol = n_t + 1, 
                                    dimnames = list(paste("ind"  , 
                                                          1:n_i, sep = " "), 
                                                    paste("cycle",
                                                          0:n_t, sep = " "))) 
       v_M_init  <- rep("Cov19+", n_i)
       m_M[, 1] <- v_M_init          # initial health state
       m_C[, 1] <- Costs(l_params, m_M[, 1])   # costs accrued during cycle 0
       m_E[, 1] <- Effs(l_params, m_M[, 1])    # LDs accrued during cycle 0
       
       # open a loop for time running cycles 1 to n_t 
       for (t in 1:n_t) { # t <- 1
         # calculate the transition probabilities for the cycle based on health state t
         m_P <- Probs(v_M_t = m_M[, t], df_X, t = t, Trt = Trt) 
         # sample the current health state and store that state in matrix m_M
         m_M[, t + 1]  <- samplev(m_P, 1)    
         # calculate costs per individual during cycle t + 1
         m_C[, t + 1]  <- Costs(l_params, m_M[, t + 1], Trt = Trt)  
         # calculate QALYs per individual during cycle t + 1
         m_E[, t + 1]  <- Effs(l_params, m_M[, t + 1])  
         
       } # close the loop for the time points 
       
       if (Trt == "Remd"){
         m_C[, 11:n_t + 1]  <- m_C[, 11:n_t + 1] - c_remd
       } else if (Trt == "Remd_Ba"){
         m_C[, 11:n_t + 1]  <- m_C[, 11:n_t + 1] - c_remd
         m_C[, 15:n_t + 1]  <- m_C[, 15:n_t + 1] - c_bari
       } else {
         m_C <- m_C
       }
       
       m_C[m_C < 0] <- 0
       
       # calculate  
       tc <- m_C %*% v_dwc    # total (discounted) cost per individual
       te <- m_E %*% v_dwe    # total LDs per individual 
       tc_hat <- mean(tc)     # average (discounted) cost 
       te_hat <- mean(te)     # average LDs per individual
       tc_sum <- sum(tc)      # sum (discounted) cost
       te_sum <- sum(te)      # sum LDs per individual
       
       # store the results from the simulation in a list
       results <- list(m_M = m_M, 
                       tc_hat = tc_hat, 
                       te_hat = te_hat,
                       tc_sum = tc_sum,
                       te_sum = te_sum)   
                      
       return(results)  # return the results
  
     } # end of the MicroSim function  
  )
}

CEA_function <- function(l_params, df_X, df_h, n_wtp = 10267*2*22.10){ # User defined
  with(as.list(l_params), {
    ## Create discounting vectors
    v_names_str = c("No Treatment", "Remdesivir", "Remdesivir & Bariticinib")
    ## Run STM model at a parameter set for each intervention
    Res_model          <- CEA_microsim(l_params = l_params, df_X,
                                     df_h, Trt = FALSE)
    Res_model_remd     <- CEA_microsim(l_params = l_params, df_X,
                                     df_h, Trt = "Remd")
    Res_model_remd_ba  <- CEA_microsim(l_params = l_params, df_X, 
                                     df_h, Trt = "Remd_Ba")

    ## Vector with total discounted mean Costs and LDs
    v_tc_d      <- c(Res_model$tc_hat, Res_model_remd$tc_hat, 
                     Res_model_remd_ba$tc_hat)
    v_tu_d      <- c(Res_model$te_hat, Res_model_remd$te_hat, 
                     Res_model_remd_ba$te_hat)
    
    ## Vector with discounted net monetary benefits (NMB)
    v_nmb_d     <- v_tu_d * n_wtp - v_tc_d
    
    ## Dataframe with discounted costs, effectiveness and NMB
    df_cea <- data.frame(Strategy = v_names_str,
                        Cost     = v_tc_d,
                        Effect   = v_tu_d,
                        NMB      = v_nmb_d)
    return(df_cea)
  }
  )
}
