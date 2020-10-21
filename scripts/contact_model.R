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

# Set local path ----------------------------------------------

# setwd("~/Documents/GitHub/2020-cov-tracing/")

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

# - - - - - - 
# MAIN FIGURES
# - - - - - - 

# - - - - - - 
# Baseline case (Table 3 and 4):
#offspring_model(n_run = n_run_pick, range_n = c(1:13), dir_pick = out_dir,output_r = T)
#table_outputs_1(dir_pick = out_dir)


# - - - - - - 
# Baseline with reduced contacts

covid_s <- 0

offspring_model(n_run = n_run_pick, 
                wfh_prob = 0,
                wfh_probC = 0,
                other_prob = 0,
                cc_risk = 0.06*(1-covid_s), # Outside HH contact risk
                trace_prop = 0.9, # Proportion of contacts traced
                pre_inf = 2, # Pre-infectious period
                test_delay = 2, # Delay to tracing
                trace_delay = 2, # Delay to tracing
                p_tested = 0.8, # Probability index case isolated/tested
                trace_adherence = 0.8, # Adherence of traced contacts to quarantine
                range_n = c(1,2,3,7), dir_pick = out_dir,output_r = T)


# - - - - - - 
# Iterate over different limits on contacts and with/without limit on children
limit_range <- round(10^seq(0,2,0.2)) %>% unique()
max_lim_children_vals <- c(T,F)

foreach(ii = limit_range) %dopar% {
  for(jj in max_lim_children_vals){
    offspring_model(max_low_fix = ii,
                    max_lim_children=jj,
                    range_n = c(6),
                    wfh_prob=0.5,
                    n_run = 2e4,
                    pre_inf = 2, # Pre-infectious period
                    trace_prop = 0.5, # Proportion of contacts traced
                    test_delay = 3, # Delay to testing
                    trace_delay = 2, # Delay to tracing
                    p_tested = 0.8, # Probability index case isolated/tested
                    trace_adherence = 0.5, # Adherence of traced contacts to quarantine
                    dir_pick = paste0(out_dir,"runs_contacts/",jj,"_"))
  }
}

plot_contacts()


# - - - - - - 
# Iterate over different tracing speeds
tracing_range <- seq(0,1,0.2)
speed_range_test <- c(0,1,2,3,5)
speed_range_trace <- c(1:3)
t_asymp_range <- c(0.25,0.5,0.75)
p_test_range <- c(0.4,0.6,0.8)

isolate_choose <- c(0,0,0.25,0.25,0.25,0.25)

#for(ii in trace_range){
foreach(jj = speed_range_test) %dopar% {
  for(pp in tracing_range){
      for(kk in speed_range_trace){
        for(mm in t_asymp_range){
          for(tt in p_test_range){
              offspring_model(range_n = c(7),
                            isolate_distn = isolate_choose, # onset to isolation
                            trace_prop = pp, # Proportion of contacts traced
                            test_delay = jj, # Delay to test results that trigger tracing of contacts
                            trace_delay = kk, # Delay to tracing
                            p_tested = tt, # Probability index case isolated/tested
                            trace_adherence = 1, # Adherence of traced contacts to quarantine (wrap into trace parameter)
                            prob_t_asymp = mm, # Asymptomatic contribution
                            n_run = n_run_pick,
                            dir_pick = paste0(out_dir,"runs_delay/out",jj,"_",pp,"_",kk,"_",mm,"_",tt))
          }
        }
      }
  }
}

# - - - - - - 
# Output combined table with values
table_outputs_delay_range(dir_pick = out_dir)


# - - - - - - 
# Sensitivity analysis on presymptomatic and isolation:
offspring_model(n_run = n_run_pick, range_n = c(1:13),isolate_distn = c(0,0,0.25,0.25,0.2,0.3),dir_pick = paste0(out_dir,"sensitivity/late_detection_"))
offspring_model(n_run = n_run_pick, range_n = c(1:13),isolate_distn = c(0,0.8,0.2,0,0,0),dir_pick = paste0(out_dir,"sensitivity/fast_isolation_"))


# Output mean delay from onset-to-isolation in baseline and delayed scenario
sum(c(0,0.25,0.25,0.2,0.3,0)*(0:5))
sum(c(0,0.8,0.2,0,0,0)*(0:5))
sum(c(0,0,0.25,0.25,0.2,0.3)*(0:5))

# - - - - - - 
# Sensitivity analysis on higher non-household contact SAR
offspring_model(n_run = n_run_pick, range_n = c(1:13), cc_risk = 0.07, dir_pick = paste0(out_dir,"sensitivity/CC_SAR_higher_"))
offspring_model(n_run = n_run_pick, range_n = c(1:13), hh_risk = 0.5, cc_risk = 0.05, dir_pick = paste0(out_dir,"sensitivity/HH_SAR_higher_"))


# - - - - - - 
# Iterate over different proportions traced
trace_range <- seq(0,1,0.1)

#for(ii in trace_range){
foreach(ii = trace_range) %dopar% {
  offspring_model(trace_prop = ii, range_n = c(4:7),wfh_prob=0,n_run = n_run_pick,dir_pick = paste0(out_dir,"runs/"))
}

# - - - - - - 
# Iterate over different limits on contacts and WFH
limit_range <- round(10^seq(0,2,0.25))

foreach(ii = limit_range) %dopar% {
  offspring_model(max_low_fix = ii,range_n = c(5,6,8,9),wfh_prob=0,n_run = n_run_pick,dir_pick = paste0(out_dir,"runs/"))
  offspring_model(max_low_fix = ii,range_n = c(5,6,8,9),wfh_prob=0.5,n_run = n_run_pick,dir_pick = paste0(out_dir,"runs/"))
}

# - - - - - - 
# Iterate over different WFH range

wfh_range <- seq(0,0.6,0.1)

foreach(ii = wfh_range) %dopar% {
  offspring_model(range_n = c(5,7,8),wfh_prob=ii,n_run = n_run_pick,dir_pick = paste0(out_dir,"runs/"))
}

# With schools closed:

foreach(ii = wfh_range) %dopar% {
  offspring_model(range_n = c(5,7,8),wfh_prob=ii,wfh_probC=1, n_run = n_run_pick,dir_pick = paste0(out_dir,"runs/SC_"))
}

# - - - - - - 
# Iterate over cell-phone range

app_range <- seq(0,1,0.1)

foreach(ii = app_range) %dopar% {
  offspring_model(range_n = c(8,9),wfh_prob=0,n_run = n_run_pick,app_cov = ii, dir_pick = paste0(out_dir,"runs/"))
}

# - - - - - - 
# SUPPLEMENTARY FIGURES
# - - - - - - 

# Iterate over prop symptomatic
# Note: children scaled according to adult proportion (25%/60%)

prop_asymp <- seq(0.2,0.8,0.1)

foreach(ii = prop_asymp) %dopar% {
  offspring_model(range_n = c(1,5:9),n_run = n_run_pick,p_symptomaticC = ii*(0.3/0.7), p_symptomaticA = ii, dir_pick = paste0(out_dir,"runs2/"))
}

# - - - - - - 
# Iterate over asymptomatic contribution

t_asymp <- seq(0,1,0.1)

foreach(ii = t_asymp) %dopar% {
  offspring_model(range_n = c(1,5:9),n_run = n_run_pick,prob_t_asymp = ii, dir_pick = paste0(out_dir,"runs2/"))
}



# Plot outputs ------------------------------------------------------------

# - - - - - - 
# Plot contact networks (Fig 1)
plot_networks(dir_pick = out_dir)

# - - - - - - 
# Plot outputs for different combinations (Fig 2)
plot_contacts(dir_pick = out_dir)

# - - - - - - 
# Plot R distribution (Fig S1)
plot_R_distribution(dir_pick = out_dir)

# - - - - - - 
# Plot outputs for different asymptomatic (Fig S2)
plot_symptom_reduction(dir_pick = out_dir)

# - - - - - - 
# Output combined main table
table_outputs_1(dir_pick = out_dir)

# - - - - - - 
# Output combined supplementary tables (Table S1-S2)
table_outputs(dir_pick = out_dir)


# Deprectated ------------------------------------------------------------

# # - - - - - - 
# Optional scenarios:
# offspring_model(n_run = n_run_pick,range_n = c(1,2,3,11),pt_extra = 0.5,pt_extra_reduce = 0.3,dir_pick = "extra1_")

# - - - - - - 
# 2020_06_22

# - - - - - - 
# Iterate over different proportions of contacts
contact_range <- seq(0,0.8,0.1)
reduction_range <- c(0,0.05,0.25,0.5)

#for(ii in trace_range){
foreach(ii = contact_range) %dopar% {
  for(pp in c(0.2,0.35,0.8)){
    for(kk in c(0)){
      for(jj in reduction_range){
        offspring_model(range_n = c(1,3,7),
                        wfh_prob = ii,
                        wfh_probC = kk, # school contacts
                        other_prob = ii,
                        cc_risk = 0.06*(1-jj), # Outside HH contact risk
                        trace_prop = pp, # Proportion of contacts traced
                        trace_delay = 1, # Delay to tracing
                        p_tested = 0.5, # Probability index case isolated/tested
                        trace_adherence = 0.8, # Adherence of traced contacts to quarantine
                        n_run = n_run_pick,dir_pick = paste0(out_dir,"runs/"))
      }
    }
  }
}

# - - - - - - 
# Output combined table with interactions
table_outputs_contact_range(dir_pick = out_dir)


