# - - - - - - - - - - - - - - - - - 
# Analysis functions for model of isolation, contact tracing, and physical distancing
# Adam Kucharski (2020)
# https://github.com/adamkucharski/2020-cov-tracing
# - - - - - - - - - - - - - - - - - 

library(tidyverse)
library(igraph)
library(doMC)

registerDoMC(4)  #change the 4 to your number of CPU cores


# Set local path ----------------------------------------------

wdir <- "~/GitHub/2020-cov-tracing/" # Set git path

setwd(wdir)
out_dir <- paste0(wdir,"outputs/") # Set data output path


# Create output directories if don't exist
subdir_list <- c("runs","runs2","sensitivity")
for(ii in 1:length(subdir_list)){
  ifelse(!dir.exists(file.path(out_dir, subdir_list[ii])), dir.create(file.path(out_dir, subdir_list[ii])), FALSE)
}

# Load contact data ----------------------------------------------

data_user_col_red_u18 <- read_csv("data/contact_distributions_u18.csv")
data_user_col_red_o18 <- read_csv("data/contact_distributions_o18.csv")


# Run simulation model ----------------------------------------------

hh_risk <- 0.2 # HH risk
cc_risk <- 0.06 # contact risk
baseline_incident_cases <- 1e5 # baseline incidence symptomatic cases

n_run_pick <- 1e4 # model iterations

# Load model functions
source("R/model_functions.R")

# - - - - - - 
# Baseline case:
offspring_model(n_run = n_run_pick)

# - - - - - - 
# Sensitivity analysis on presymptomatic and isolation:
offspring_model(n_run = n_run_pick,isolate_distn = c(0.25,0.25,0.2,0.3,0,0),dir_pick = "sensitivity/no_presym_")
offspring_model(n_run = n_run_pick,isolate_distn = c(0,0,0.25,0.25,0.2,0.3),dir_pick = "sensitivity/late_detection_")

# Output mean delay from onset-to-isolation in baseline and delayed scenario
sum(c(0,0.25,0.25,0.2,0.3,0)*(0:5))
sum(c(0,0,0.25,0.25,0.2,0.3)*(0:5))

# - - - - - - 
# Iterate over different proportions traced
trace_range <- seq(0,1,0.1)

#for(ii in trace_range){
foreach(ii = trace_range) %dopar% {
  offspring_model(trace_prop = ii, range_n = c(4:7),wfh_prob=0,n_run = n_run_pick,dir_pick = "runs/")
}

# - - - - - - 
# Iterate over different limits on contacts and WFH
limit_range <- round(10^seq(0,2,0.25))

foreach(ii = limit_range) %dopar% {
  offspring_model(max_low_fix = ii,range_n = c(5,6,8,9),wfh_prob=0,n_run = n_run_pick,dir_pick = "runs/")
  offspring_model(max_low_fix = ii,range_n = c(5,6,8,9),wfh_prob=0.5,n_run = n_run_pick,dir_pick = "runs/")
}

# - - - - - - 
# Iterate over different WFH range

wfh_range <- seq(0,0.8,0.1)

foreach(ii = wfh_range) %dopar% {
  offspring_model(range_n = c(5,7,8),wfh_prob=ii,n_run = n_run_pick,dir_pick = "runs/")
}

# - - - - - - 
# Iterate over cell-phone range

app_range <- seq(0,1,0.1)

foreach(ii = app_range) %dopar% {
  offspring_model(range_n = c(8,9),wfh_prob=0,n_run = n_run_pick,app_cov = ii, dir_pick = "runs/")
}


# - - - - - - 
# Iterate over prop symptomatic

prop_asymp <- seq(0.2,0.8,0.1)

foreach(ii = prop_asymp) %dopar% {
  offspring_model(range_n = c(1,5:9),n_run = n_run_pick,prob_symp = ii, dir_pick = "runs2/")
}

# - - - - - - 
# Iterate over asymptomatic contribution

t_asymp <- seq(0,1,0.1)

foreach(ii = t_asymp) %dopar% {
  offspring_model(range_n = c(1,5:9),n_run = n_run_pick,prob_t_asymp = ii, dir_pick = "runs2/")
}



# Plot outputs ------------------------------------------------------------

# - - - - - - 
# Plot contact networks
plot_networks()

# - - - - - - 
# Plot outputs for different combinations
plot_contacts()

# - - - - - - 
# Plot outputs for different combinations
plot_symptom_reduction()


# Deprectated ------------------------------------------------------------

# # - - - - - - 
# Optional case:
# offspring_model(n_run = n_run_pick,range_n = c(1,2,3,11),pt_extra = 0.5,pt_extra_reduce = 0.3,dir_pick = "extra1_")





