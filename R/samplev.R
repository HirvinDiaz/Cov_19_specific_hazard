#---------------------------------------------------------------------------#
#### R function to sample states for multiple individuals simultaneously ####
#---------------------------------------------------------------------------#

# Developed by the Decision Analysis in R for Technologies in Health (DARTH) workgroup
# Fernando Alarid-Escudero, PhD (1) 
# Eva A. Enns, MS, PhD (2)	
# M.G. Myriam Hunink, MD, PhD (3,4)
# Hawre J. Jalal, MD, PhD (5) 
# Eline M. Krijkamp, MSc (3)	
# Petros Pechlivanoglou, PhD (6) 

# In collaboration of: 		
# 1 Drug Policy Program, Center for Research and Teaching in Economics (CIDE) - CONACyT, 
#   Aguascalientes, Mexico
# 2 University of Minnesota School of Public Health, Minneapolis, MN, USA
# 3 Erasmus MC, Rotterdam, The Netherlands
# 4 Harvard T.H. Chan School of Public Health, Boston, USA
# 5 University of Pittsburgh Graduate School of Public Health, Pittsburgh, PA, USA
# 6 The Hospital for Sick Children, Toronto and University of Toronto, Toronto ON, Canada

#####################################################################################
# Please cite our publications when using this code
# - Jalal H, Pechlivanoglou P, Krijkamp E, Alarid-Escudero F, Enns E, Hunink MG. 
# An Overview of R in Health Decision Sciences. Med Decis Making. 2017; 37(3): 735-746. 
# https://journals.sagepub.com/doi/abs/10.1177/0272989X16686559
# - Krijkamp EM, Alarid-Escudero F, Enns EA, Jalal HJ, Hunink MGM, Pechlivanoglou P. 
# Microsimulation modeling for health decision sciences using R: A tutorial. 
# Med Decis Making. 2018;38(3):400–22. 
# https://journals.sagepub.com/doi/abs/10.1177/0272989X18754513

#####################################################################################
# Copyright 2017, THE HOSPITAL FOR SICK CHILDREN AND THE COLLABORATING INSTITUTIONS. 
# All rights reserved in Canada, the United States and worldwide. Copyright, 
# trademarks, trade names and any and all associated intellectual property are 
# exclusively owned by THE HOSPITAL FOR Sick CHILDREN and the collaborating 
# institutions. These materials may be used, reproduced, modified, distributed 
# and adapted with proper attribution.
#####################################################################################

samplev <- function(m_Probs, m) {
  # Arguments
  # m_Probs: matrix with probabilities (n.i * n.s)
  # m:       number of states than need to be sampled per individual  
  # Return
  # ran:    n.i x m matrix filled with sampled health state(s) per individual
  
  d <- dim(m_Probs)  # dimensions of the matrix filled with the multinomical probabilities for the health states 
  n <- d[1]          # first dimension - number of rows (number of individuals to sample for)
  k <- d[2]          # second dimension - number of columns (number of health states considered)
  lev <- dimnames(m_Probs)[[2]]  # extract the names of the health states considered for sampling
  if (!length(lev))  # in case names for the health states are missing, use numbers to specify the health states
    lev <- 1:k       # create a sequence from 1:k (number of health states considered)
  # create a matrix 
  ran <- matrix(lev[1], ncol = m, nrow = n) # create the matrix ran, filled with the first health state of the levels 
  U <- t(m_Probs)    # transposed m.Probs matrix n.i x n.s --> n.s x n.i 
  
  for(i in 2:k) {    # start loop, from the 2nd health states
    U[i, ] <- U[i, ] + U[i - 1, ] # start summing the probabilities of the different health states per individual 
  }
  if (any((U[k, ] - 1) > 1e-05))  # sum of all probs per individual - 1 should be 0 (use 1e-05 for rounding issues), else print the error statement
    stop("error in multinom: probabilities do not sum to 1")
  
  for (j in 1:m) {   # start loop of the state that needs to be sampled (m)
    un <- rep(runif(n), rep(k, n))       # sample from a uniform distribution of length n*k
    ran[, j] <- lev[1 + colSums(un > U)] # store the health state at the jth column of the U matrix
  }
  ran # return the new health state per individual n.i x m
} # close the function 
