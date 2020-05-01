# - - - - - - - - - - - - - - - - - 
# Analysis functions for model of isolation, contact tracing, and physical distancing
# Adam Kucharski (2020)
# https://github.com/adamkucharski/2020-cov-tracing
# - - - - - - - - - - - - - - - - - 

library(dplyr)
library(tibble)
library(readr)
library(magrittr)
library(igraph)
library(doMC)

registerDoMC(4)  #change the 4 to your number of CPU cores

setwd("~/Documents/GitHub/2020-cov-tracing/")

# Set local path ----------------------------------------------

wdir <- getwd()
out_dir <- paste0(wdir,"/outputs/") # Set data output path

# Create output directories if don't exist
subdir_list <- c("runs","runs2","sensitivity","rr_out")
for(ii in 1:length(subdir_list)){
  ifelse(!dir.exists(paste0(out_dir, subdir_list[ii])), dir.create(file.path(out_dir, subdir_list[ii])), FALSE)
}

# Load contact data ----------------------------------------------

data_user_col_red_u18 <- read_csv("data/contact_distributions_u18.csv")
data_user_col_red_o18 <- read_csv("data/contact_distributions_o18.csv")

n_user_u18 <- nrow(data_user_col_red_u18)
n_user_o18 <- nrow(data_user_col_red_o18)


# Run simulation model ----------------------------------------------
# For each function, outputs are saved in 'dir_pick' directory.

n_run_pick <- 1e4 # model iterations

set.seed(201)

# Load model functions
source("R/model_functions.R")

scen_choose <- c(2,3,5,7,8)
isolate_choose <- c(0,0,0.25,0.25,0.25,0.25)
pre_inf_choose <- 2
inf_period_choose <- 5
non_risk_pick <- 0.002


# - - - - - - 
# PCR test probability 

pos_test <- function(x){
  if(x==0){y=0}
  if(x==1){y=0.05}
  if(x==2){y=0.2}
  if(x==3){y=0.4}
  if(x==4){y=0.5}
  if(x>=5){y=0.6}
  y
}

# pos_test <- function(x){
#   if(x==0){y=0.05}
#   if(x==1){y=0.2}
#   if(x==2){y=0.5}
#   if(x==3){y=0.6}
#   if(x==4){y=0.9}
#   if(x>=5){y=1.0}
#   y
# }

prob_test <- 0.9*pos_test(5)

# - - - - - - 
# Baseline case:
offspring_model(n_run = n_run_pick, range_n = scen_choose,
                isolate_distn = isolate_choose,
                inf_period = inf_period_choose, sample_delay = 0, test_delay = 2, pre_inf = pre_inf_choose,
                non_risk = non_risk_pick,
                dir_pick = paste0(out_dir,"sensitivity/delay0_"),
                p_tested = prob_test
                )

# - - - - - - 
# Sensitivity case:

offspring_model(n_run = n_run_pick, range_n = scen_choose,
                isolate_distn = isolate_choose,
                inf_period = inf_period_choose, sample_delay = 3, test_delay = 2, pre_inf = pre_inf_choose,
                non_risk = non_risk_pick,
                dir_pick = paste0(out_dir,"sensitivity/delay3_"),
                p_tested = prob_test
)


offspring_model(n_run = n_run_pick, range_n = scen_choose,
                isolate_distn = isolate_choose,
                inf_period = inf_period_choose, sample_delay = 6, test_delay = 2, pre_inf = pre_inf_choose,
                non_risk = non_risk_pick,
                dir_pick = paste0(out_dir,"sensitivity/delay6_"),
                p_tested = prob_test
)

offspring_model(n_run = n_run_pick, range_n = scen_choose,
                isolate_distn = isolate_choose,
                inf_period = inf_period_choose, sample_delay = 12, test_delay = 2, pre_inf = pre_inf_choose,
                non_risk = non_risk_pick,
                dir_pick = paste0(out_dir,"sensitivity/delay14_"),
                p_tested = prob_test
)



# Plot outputs ------------------------------------------------------------

# - - - - - - 
# Output combined supplementary tables (Table S1-S2)
table_outputs(dir_pick = out_dir)


# Calculate mean
sum(c(0,0,0.25,0.25,0.25,0.25)*(0:5))-2


# - - - - - - 
# Plot R distribution (Fig S1)
plot_R_distribution(dir_pick = out_dir)




# - - - - - - 
# Plot contact networks (Fig 1)
plot_networks(dir_pick = out_dir)



# - - - - - - 
# Plot outputs for different combinations (Fig 2)
plot_contacts(dir_pick = out_dir)



# - - - - - - 
# Plot outputs for different asymptomatic (Fig S2)
plot_symptom_reduction(dir_pick = out_dir)




# Deprectated ------------------------------------------------------------

# # - - - - - - 
# Optional scenarios:
# offspring_model(n_run = n_run_pick,range_n = c(1,2,3,11),pt_extra = 0.5,pt_extra_reduce = 0.3,dir_pick = "extra1_")





