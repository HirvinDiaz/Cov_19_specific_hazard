#### Final project script. Author: Hirvin Azael Diaz Zepeda ####

####  Clean Work - space  ####
rm(list = ls()) # to clean the work - space

####  Load Packages  ####
library(dplyr)
library(ggplot2)
library(utils)
library(scales)
library(dampack)
library(chron)
library(tibble)
library(knitr)
library(pander)
library(kableExtra)
library(CEAutil)
library(tidyr)
library(dampack)
library(data.table)
library(reshape2)
library(survival)
library(survPen)
library(relsurv)
library(readr)
library(readxl)
library(lubridate)
library(muhaz)
library(ggpubr)

######  Load data bases  ######

# Data base of suspected people with COVID-19 published by Mexico´s Ministry of 
# Health 

load("data/cov_14_01.Rdata") # data base of november 15th






