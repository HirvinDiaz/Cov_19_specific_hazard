#####################################################################################
## Function to make a simple plot from the data generated in a Microsimulation model ##
#####################################################################################

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

# Histogram showing variability in individual total costs
plot(density(tc), main = paste("Total cost per person"), xlab = "Cost ($)")

# Histogram showing variability in individual total QALYs
plot(density(te), main = paste("Total QALYs per person"), xlab = "QALYs")


# PLOT THE DISTRIBUTION OF THE POPULATON ACROSS HEALTH STATES OVER TIME (TRACE)
# count the number of individuals in each health state at each cycle
m_TR <- t(apply(m_M, 2, function(x) table(factor(x, levels = v_n, ordered = TRUE)))) 
m_TR <- m_TR / n_i                                       # calculate the proportion of individuals 
colnames(m_TR) <- v_n                                    # name the rows of the matrix
rownames(m_TR) <- paste("Cycle", 0:n_t, sep = " ")       # name the columns of the matrix

# Plot trace of first health state
plot(0:n_t, m_TR[, 1], type = "l", main = "Health state trace", 
     ylim = c(0, 1), ylab = "Proportion of cohort", xlab = "Cycle")
# add a line for each additional state
for (n_s in 2:length(v_n)) {
  lines(m_TR[, n_s], col = n_s)       # adds a line to current plot
}
legend("topright", v_n, col = 1:3,    # add a legend to current plot
       lty = rep(1, 3), bty = "n", cex = 0.65)


