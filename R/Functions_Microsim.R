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

CEA_microsim <- function(l_params, df_X, df_h, df_pop_ch, life_expectancy, 
                         Trt = FALSE) { 
  with(as.list(l_params), {

       d_h_HD_hosp <- df_h # Data frame with specific hazards
       
       # Transform hazards to transition probabilities
       d_h_HD_hosp <- d_h_HD_hosp %>% 
         mutate(haz_Remd = exp(log(Hazard.Covid_19) + log(hr_Remd))) %>% 
         mutate(haz_Remd_Ba = exp(log(haz_Remd) + log(hr_RemdBa_vs_Remd)))
       
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
       
       # Set the data table to be indexed by age, day, sex and ICU
       setkey(dt_p_CoV, group, time, sex)

       # States
       v_names_states   <- c("Cov19+", "CoV19_Dead", "O_Causes_Dead") # state names
       n_states         <- length(v_names_states)   # number of states
       
       # Number of individuals and cycles
       n_i <- length(df_X$sex)
       n_t <- 50
       
       # vector with discount weights for costs
       v_dwc <- 1 / (1) ^ (0:n_t) 
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
       
       Final_outcome <- as.data.frame(m_E) 
       Final_outcome <- Final_outcome %>% 
         select(`cycle 50`)
       
       Effects <- cbind(Final_outcome,df_pop_ch)
       
       Effects <- Effects %>% 
         merge(Life_expectancy, 
               by = c("sex" = "sex",
                      "age" = "age")) %>% 
         mutate(years = ifelse(`cycle 50` == 0, 0, Years_dis))
       
       # calculate  
       tc <- m_C %*% v_dwc    # total (discounted) cost per individual
       te <- m_E %*% v_dwe    # total LDs per individual 
       tc_hat <- mean(tc)     # average (discounted) cost 
       te_hat <- mean(te)     # average LDs per individual
       lf_hat <- mean(Effects$years) # average Life years Gained
       
       # store the results from the simulation in a list
       results <- list(tc_hat = tc_hat, 
                       te_hat = te_hat,
                       lys_g  = lf_hat)   
                      
       return(results)  # return the results
  
     } # end of the MicroSim function  
  )
}

CEA_function <- function(l_params, df_X, df_h, df_pop_ch, life_expectancy,
                         n_wtp = 10267*2*22.10){ # User defined
  with(as.list(l_params), {
    ## Create discounting vectors
    v_names_str = c("No Treatment", "Remdesivir", "Remdesivir & Bariticinib")
    ## Run STM model at a parameter set for each intervention
    model_no_trt  <- CEA_microsim(l_params = l_params, df_X,
                             df_h, df_pop_ch, life_expectancy, Trt = FALSE)
    model_remd    <- CEA_microsim(l_params = l_params, df_X,
                             df_h, df_pop_ch, life_expectancy, Trt = "Remd")
    model_remd_ba <- CEA_microsim(l_params = l_params, df_X, 
                             df_h, df_pop_ch, life_expectancy, Trt = "Remd_Ba")

    ## Vector with total discounted mean Costs and Lys
    v_tc_d      <- c(model_no_trt$tc_hat, model_remd$tc_hat, 
                     model_remd_ba$tc_hat)
    v_tu_d      <- c(model_no_trt$lys_g, model_remd$lys_g, 
                     model_remd_ba$lys_g)
    
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

#### Functions for ICU ####

Costs_icu <- function (l_params, v_M_t, Trt = FALSE) {
  with(as.list(l_params),{
  # v_M_t: current health state
  c_t <- c()
  if (Trt == "Dexa"){
    c_t[v_M_t == "Cov19+"]  <- c_resp + c_dexa
  } else if (Trt == "None" | Trt == FALSE){
    c_t[v_M_t == "Cov19+"]  <- c_resp 
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

CEA_microsim_icu <- function(l_params, df_X, df_h, df_pop_ch, life_expectancy,
                             Trt = FALSE) { 
  with(as.list(l_params), {
   
    d_h_HD_resp <- df_h # Data frame with specific hazards
    
    d_h_HD_resp <- d_h_HD_resp %>% 
      mutate(haz_Dexa = exp(log(Hazard.Covid_19) + log(hr_Dexa)))
    
    d_p_HD_resp <- d_h_HD_resp %>% 
      mutate(p_dCoV = 1 - exp(-Hazard.Covid_19)) %>% 
      mutate(p_dPop = 1 - exp(- Hazard.Population)) %>% 
      mutate(p_dCoV_Dexa = 1 - exp(- haz_Dexa)) %>%  
      arrange(group, time)%>% 
      select(c("time", "group", "sex", "p_dCoV", "p_dPop", "p_dCoV_Dexa"))
    
    # Convert data frame to a data table for efficiency
    dt_p_CoV <- data.table(d_p_HD_resp)
    
    # set the data table to be indexed by age, day, sex and ICU
    setkey(dt_p_CoV, group, time, sex)
    
    # States
    v_names_states   <- c("Cov19+", "CoV19_Dead", "O_Causes_Dead") # state names
    n_states         <- length(v_names_states)   # number of states
    
    # Number of individuals and cycles
    n_i <- length(df_X$sex)
    n_t <- 50
    
    # vector with discount weights for costs
    v_dwc <- 1 / (1) ^ (0:n_t) 
    v_dwe <- 1 / (1) ^ (0:n_t) 
    
    Probs <- function(v_M_t, df_X, t, Trt = FALSE) { # t <- 1
      m_p_t           <- matrix(0, nrow = n_states, ncol = n_i)  
      # give the state names to the rows
      rownames(m_p_t) <-  v_names_states
      # Lookup baseline probability of dying from Covid-19 or other causes  
      if (Trt == "Dexa"){
        p_die_CoV_all <- dt_p_CoV[.(df_X$group, df_X$time + t, df_X$sex), p_dCoV_Dexa] 
      } else if (Trt == "None" | Trt == FALSE){
        p_die_CoV_all <- dt_p_CoV[.(df_X$group, df_X$time + t, df_X$sex), p_dCoV] 
      } else {
        message("Choose a treatment within the options")
      }
      p_die_Pop_all <- dt_p_CoV[.(df_X$group, df_X$time + t, df_X$sex), p_dPop] 
      
      p_die_CoV     <- p_die_CoV_all[v_M_t == "Cov19+"]  
      p_die_Pop     <- p_die_Pop_all[v_M_t == "Cov19+"] 
      # update m_p_t with the appropriate probabilities   
      # transition probabilities when healthy
      m_p_t[, v_M_t == "Cov19+"] <- rbind(1 - (p_die_CoV + p_die_Pop), p_die_CoV, p_die_Pop)    
      # transition probabilities when sick 
      m_p_t[, v_M_t == "CoV19_Dead"] <- rbind(0, 1 ,0)
      # transition probabilities when sicker 
      m_p_t[, v_M_t == "O_Causes_Dead"] <- rbind(0, 0, 1)
      return(t(m_p_t))
    }     
    # Matrix
    m_M <- m_C <- m_E <-  matrix(NA, nrow = n_i, ncol = n_t + 1, 
                                 dimnames = list(paste("ind"  , 
                                                       1:n_i, sep = " "), 
                                                 paste("cycle",
                                                       0:n_t, sep = " "))) 
    v_M_init  <- rep("Cov19+", n_i)
    m_M[, 1] <- v_M_init          # initial health state
    m_C[, 1] <- Costs_icu(l_params, m_M[, 1])   # costs accrued during cycle 0
    m_E[, 1] <- Effs(l_params, m_M[, 1])    # LDs accrued during cycle 0
    
    # open a loop for time running cycles 1 to n_t 
    for (t in 1:n_t) { # t <- 1
      # calculate the transition probabilities for the cycle based on health state t
      m_P <- Probs(v_M_t = m_M[, t], df_X, t = t, Trt = Trt) 
      # sample the current health state and store that state in matrix m_M
      m_M[, t + 1]  <- samplev(m_P, 1)    
      # calculate costs per individual during cycle t + 1
      m_C[, t + 1]  <- Costs_icu(l_params, m_M[, t + 1], Trt = Trt)  
      # calculate QALYs per individual during cycle t + 1
      m_E[, t + 1]  <- Effs(l_params, m_M[, t + 1])  
      
    } # close the loop for the time points 
  
  if (Trt == "Dexa"){
    m_C[, 11:n_t + 1]  <- m_C[, 11:n_t + 1] - c_dexa
  } else {
    m_C <- m_C
  }
  
  m_C[m_C < 0] <- 0
  
  Final_outcome <- as.data.frame(m_E) 
  Final_outcome <- Final_outcome %>% 
    select(`cycle 50`)
  
  Effects <- cbind(Final_outcome,df_pop_ch)
  
  Effects <- Effects %>% 
    merge(Life_expectancy, 
          by = c("sex" = "sex",
                 "age" = "age")) %>% 
    mutate(years = ifelse(`cycle 50` == 0, 0, Years_dis))
    
    # calculate  
    tc <- m_C %*% v_dwc    # total (discounted) cost per individual
    te <- m_E %*% v_dwe    # total LDs per individual 
    tc_hat <- mean(tc)     # average (discounted) cost 
    te_hat <- mean(te)     # average LDs per individual
    tc_sum <- sum(tc)      # sum (discounted) cost
    te_sum <- sum(te)      # sum LDs per individual
    lf_hat <- mean(Effects$years) # average Life years Gained
    
    # store the results from the simulation in a list
    results <- list(tc_hat = tc_hat, 
                    te_hat = te_hat,
                    lys_g  = lf_hat)   
    
    return(results)  # return the results
    
  } # end of the MicroSim function  
  )
}

CEA_function_icu <- function(l_params, df_X, df_h, df_pop_ch, life_expectancy,
                             n_wtp = 10267*2*22.10){ # User defined
  with(as.list(l_params), {
    ## Create discounting vectors
    v_names_str = c("No Treatment", "Dexamethasone")
    ## Run STM model at a parameter set for each intervention
    model_n_trt    <- CEA_microsim_icu(l_params = l_params, df_X, df_h,  
                                       df_pop_ch, life_expectancy, Trt = FALSE)
    model_Dexa     <- CEA_microsim_icu(l_params = l_params, df_X, df_h, 
                                       df_pop_ch, life_expectancy, Trt = "Dexa")

    ## Vector with total discounted mean Costs and LDs
    v_tc_d      <- c(model_n_trt$tc_hat, model_Dexa$tc_hat)
    v_tu_d      <- c(model_n_trt$lys_g, model_Dexa$lys_g)
    
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

