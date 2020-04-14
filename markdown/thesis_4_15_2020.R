
# Last edited: 04/14/2020
# This script creates figures and tables for thesis submission for committee meeting by April 17, 2020

#### Outlines ================================== ####

# 1. Control Model:
#    a. raw data plot vs model plot
#    b. tables for coefficient estimates and 95% intervals
#    c. hazard ratio between control and uninfected 
# 2. ACO Model:
#    a. raw data plot vs model plot
#    b. tables for coefficient estimates and 95% intervals
#    c. hazard ratio between males and females 
#    d. immunity measurements: scale parameters
# 3. CO Model:
#    a. raw data plot vs model plot
#    b. tables for coefficient estimates and 95% intervals
#    c. hazard ratio between males and females 
#    d. immunity measurements: scale parameters

library(tidyverse)
library(data.table)
library(magrittr)
library(ggsci)
library(kableExtra)

source("./scripts/functions.R")
raw_data <- readxl::read_excel("./data/IS_data.xlsx") %>% clean_data(c("Day after spray","Deaths"))

#### Control Model ============================== ####

#### a. raw data plot vs model plot ####

raw_control <- raw_data %>% filter(Age==14 & Treatment!="Infected")
raw_control %>% get_survival_percent() %>% plot_survival_percent()

posterior_dist <- 
  read.csv("./data/posterior_control_2020-04-13.csv") %>% 
  burn_n_thin_draws(jump = 80, burn_at = 1)



#### ACO Model ================================ ####








#### CO Model ================================= ####












