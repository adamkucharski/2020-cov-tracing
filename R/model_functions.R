# - - - - - - - - - - - - - - - - - 
# Model of isolation, contact tracing, and physical distancing
# Adam Kucharski (2020)
# https://github.com/adamkucharski/2020-cov-tracing
#
# Libraries used: dplyr, tibble, readr, magrittr, igraph, doMC
# - - - - - - - - - - - - - - - - - 


# Helper functions ----------------------------------------------

c.text<-function(x,sigF=3){
  bp1=signif(c(median(x),quantile(x,0.025),quantile(x,0.975)),sigF)
  paste(bp1[1]," (",bp1[2],"-",bp1[3],")",sep="")
}

c.nume<-function(x){
  x=x[!is.na(x)] # remove NA
  bp1=c(median(x),quantile(x,0.025),quantile(x,0.975))
  as.numeric(bp1)
}


col_def <-   list(col1="dark grey",col2=rgb(0.9,0.7,0),col3=rgb(0,0,0.8),col4=rgb(0.1,0.4,0.1),col5=rgb(1,0.4,1),col6=rgb(0.2,0,0.8),rgb(0,0.7,0.7))
col_def_F <- list(col1="grey",col2=rgb(0.9,0.7,0,0.5),col3=rgb(0,0,0.8,0.5),col5=rgb(0.1,0.4,0.1,0.5),col6=rgb(1,0.4,1,0.5),col7=rgb(0.2,0,0.8,0.5),rgb(0,0.7,0.7,0.2))

# Offspring simulation model ----------------------------------------------

offspring_model <- function(max_low_fix = 4, # Social distancing limit in these scenarios
                            max_lim_children = T, # Impose limit on children as well
                            wfh_probC = 0, # Probability children have no school contacts
                            wfh_prob = 0, # Probability adults have no work contacts
                            other_prob = 0, # Probability people have no other contacts
                            range_n = NULL, # Pick specific scenarios to run
                            trace_prop = 0.95, # Proportion of contacts traced
                            n_run = 5e3, # Number of simualtions
                            app_cov = 0.53, # App coverage
                            p_symptomaticC = 0.3, # Proportion children symptomatic
                            p_symptomaticA = 0.7, # Proportion adults symptomatic
                            prob_t_asymp = 0.5, # Relative transmission from asymptomatic
                            isolate_distn = c(0,0.25,0.25,0.2,0.3,0), # distribution of time to isolate (1st day presymp)
                            dir_pick = "", # Output directory
                            pt_extra = 0, # Optional extra transmission intervention
                            pt_extra_reduce = 0, # Reduction from extra intervention
                            output_r = F,
                            hh_risk = 0.2, # HH risk
                            cc_risk = 0.06, # Outside HH contact risk
                            inf_period = 5, # Infectious period - default 5 days
                            pre_inf = 2, # Pre infectious period - default 1 day
                            #sample_delay = 0, # Delay to testing of contacts
                            test_delay = 2, # Delay to test results that trigger tracing of contacts
                            trace_delay = 2, # Delay to tracing of contacts in manual tracing
                            non_risk = 0.002, # Background risk
                            trace_adherence = 0.9, # Adherence of traced contacts to quarantine
                            p_tested = 0.9 # Probability index case tested
                            ){
  
  
  # DEBUG
  # max_low_fix = 4; wfh_prob = 0; range_n = NULL; trace_prop = 0.95; n_run = 5e2; app_cov = 0.53; prob_symp = 0.6; prob_t_asymp = 0.5; dir_pick = ""; pt_extra = 0.95; pt_extra_reduce = 0; output_r = F
  
  # p_symptomaticA=0.7; p_symptomaticC=0.7; trace_adherence=0.8; wfh_probC=0; p_tested=0.8; other_prob = 0; range_n = c(2,3,5,7,8); isolate_distn = c(0,0,0.25,0.25,0.25,0.25); inf_period = 5; sample_delay=0; trace_delay=2; test_delay = 5; pre_inf = 1; hh_risk = 0.2; cc_risk = 0.06; non_risk=0.001
  
  
  # - - - - - - - - - - - - - - - - - - - - 
  # Define parameters across scenarios (note: some will be redefined for specific scenarios)

  # Demographic parameters
  under_18_prob <- 0.21 # Probability under 18
  
  # Transmission and baseline risk
  baseline_incident_cases <- 20e3 # baseline incidence symptomatic cases (suspected or confirmed)
  
  # Symptomatic and proportion getting tested
  # Adherence to testing/quarantine

  time_isolate <- isolate_distn # Distribution over symptomatic period
  transmission_asymp <- prob_t_asymp
  phone_coverage <- 1 # App coverage in non-app scenarios
  p_pop_test <- 0.05 # Proportion mass tested (5% per week)
  
  
  # Define default scenarios
  scenario_list <- c("no_measures","isolation_only","hh_quaratine_only","hh_work_only",
                     "isolation_manual_tracing_met_only","isolation_manual_tracing_met_limit",
                     "isolation_manual_tracing","cell_phone","cell_phone_met_limit",
                     "isolation_manual_tracing_met_only_cell",
                     "isolation_manual_tracing_met_only_cell_met_limit",
                     "pop_testing",
                     "isolation_outside") 
  
  if(is.null(range_n)){nn_choose <- 1:length(scenario_list)}else{nn_choose <- range_n}
  
  store_table_scenario <- NULL
  
  # - - - - - 
  # RUN MODEL 
  # Iterate over scenarios
  
  for(kk in nn_choose){
    
    
    # Proportion of contacts met before
    met_before_w <- 0.79 # At work. 
    met_before_s <- 0.9
    met_before_h <- 1 # Within HH
    met_before_o <- 0.52 # In other settings
    
    # Set contact limit default high to avoid censoring in default scenarios
    max_contacts <- 2e3 
    
    # Tracing parameters - need to reset
    hh_trace <- 1 # Tracing in HH
    ww_trace <- trace_prop # Tracing at work
    other_trace <- trace_prop # Tracing others
  
    # Scenario set up
    scenario_pick <- scenario_list[kk]

    # Define scenario parameters
    if(scenario_pick=="no_measures"){
      do_isolation <- F
      do_tracing <- F
      cell_yes <- F
    }
    
    if(scenario_pick=="isolation_outside"){
      do_isolation <- T
      do_tracing <- F
      cell_yes <- F
    }
  
    if(scenario_pick=="isolation_only"){
      do_isolation <- T
      do_tracing <- T
      hh_trace <- 0 # Tracing at home
      ww_trace <- 0 # Tracing at work
      other_trace <- 0 # Tracing others
      cell_yes <- F
    }
    
    if(scenario_pick=="hh_quaratine_only"){
      do_isolation <- T
      do_tracing <- T
      ww_trace <- 0 # Tracing at work
      other_trace <- 0 # Tracing others
      cell_yes <- F
    }
    
    if(scenario_pick=="hh_work_only"){
      do_isolation <- T
      do_tracing <- T
      met_before_s <- 1 # Met before
      met_before_w <- 1 # Met before
      met_before_o <- 1 # Met before
      other_trace <- 0 # Tracing others
      cell_yes <- F
    }
    
    if(scenario_pick=="isolation_manual_tracing_met_limit"){
      do_isolation <- T
      do_tracing <- T
      max_contacts <- max_low_fix
      cell_yes <- F
    }

    if(scenario_pick=="isolation_manual_tracing_met_only"){
      do_isolation <- T
      do_tracing <- T
      cell_yes <- F
    }
    
    if(scenario_pick=="isolation_manual_tracing"){
      do_isolation <- T
      do_tracing <- T
      met_before_s <- 1 # Met before
      met_before_w <- 1 # Met before
      met_before_o <- 1 # Met before
      cell_yes <- F
    }
    
    if(scenario_pick=="cell_phone"){
      do_isolation <- T
      do_tracing <- T
      ww_trace <- 0 # Tracing at work
      other_trace <- 0 # Tracing others
      phone_coverage <- app_cov
      cell_yes <- T
    }
    
    if(scenario_pick=="cell_phone_met_limit"){ # Note this includes manual tracing
      do_isolation <- T
      do_tracing <- T
      ww_trace <- 0 # Tracing at work
      other_trace <- 0 # Tracing others
      phone_coverage <- app_cov
      max_contacts <- max_low_fix
      cell_yes <- T
    }
    
    if(scenario_pick=="isolation_manual_tracing_met_only_cell"){
      do_isolation <- T
      do_tracing <- T
      phone_coverage <- app_cov
      cell_yes <- T
    }
    
    if(scenario_pick=="isolation_manual_tracing_met_only_cell_met_limit"){
      do_isolation <- T
      do_tracing <- T
      phone_coverage <- app_cov
      max_contacts <- max_low_fix
      cell_yes <- T
    }

    if(scenario_pick=="pop_testing"){
      do_isolation <- F
      do_tracing <- F
      cell_yes <- F
    }
    
    # if(scenario_pick=="pop_testing_tracing_met_only"){
    #   do_isolation <- F
    #   do_tracing <- F
    #   do_isolation <- T
    #   do_tracing <- T
    #   cell_yes <- F
    #   cell_yes <- F
    # }
    
    if(scenario_pick=="pt_extra"){
      do_isolation <- T
      do_tracing <- T
      ww_trace <- 0 # Tracing at work
      other_trace <- 0 # Tracing others
      cell_yes <- F
    }
    
    
    # - - -
    # Iterate over users
    store_r <- NULL
    
    for(ii in 1:n_run){

      # Decide if child or adult
      if(runif(1)<under_18_prob){
        adult_T <- F
        pick_user <- sample(1:n_user_u18,1)
        data_ii <- data_user_col_red_u18[pick_user,]
        met_before_w <- met_before_s
        wfh_t <- runif(1)< wfh_probC # schools closed?
        other_t <- runif(1)< other_prob  # other contacts
        
        symp_T <- runif(1)<p_symptomaticC # proportion symptomatic
        
      }else{
        adult_T <- T
        pick_user <- sample(1:n_user_o18,1)
        data_ii <- data_user_col_red_o18[pick_user,]
        wfh_t <- runif(1)< wfh_prob # work contacts?
        other_t <- runif(1)< other_prob  # other contacts
        
        symp_T <- runif(1)<p_symptomaticA # proportion symptomatic
      }

      # Simulate infectious period
      # Decide if symptomatic & tested
      phone_T <- runif(1)<phone_coverage # has phone app?
      tested_T <- runif(1)<p_tested # would be tested
      inf_scale_mass_pop <- 1
      extra_red <- 1

      # Set infectious period based on symptomatic & tested
      if(tested_T==T & symp_T==T & do_isolation==T ){
        inf_period_ii <- sample(0:inf_period,1,prob=time_isolate)
      }else{
        inf_period_ii <- inf_period
      }
      
      # Set infectious period for population testing
      if(scenario_pick=="pop_testing" & runif(1)<p_pop_test){
        inf_period_ii <- sample(c(0:inf_period,5),1) # Include extra day to reflect week
        tested_T <- T
      }

      # Set relative transmission of asymptomatics
      if(symp_T==T){
        inf_propn <- 1
      }else{
        inf_propn <- transmission_asymp
      }
      
      # Define amalgamated tracing probilities
      ww_trace_all <- ww_trace*met_before_w
      other_trace_all <- other_trace*met_before_o
      
      # Check if contacts phone traced (in cell phone scenario):
      if((cell_yes==T) & phone_T==T){
          tracePhone <- phone_coverage # Tracing at work

          ww_trace_all <- 1-(1-ww_trace_all)*(1-tracePhone) # product of both reductions
          other_trace_all <- 1-(1-ww_trace_all)*(1-tracePhone) # product of both reductions

      }else{
        tracePhone <- 0
      }
      
      # Extra transmission reduction
      if(scenario_pick=="pt_extra" & runif(1)<pt_extra){
        tested_T <- T
        extra_red <- (1-pt_extra_reduce) # Multiply by both
      }

      # Proportion infectious
      inf_ratio <- inf_period_ii/inf_period

      # Tally contacts
      home_c <- data_ii$e_home
      work_c <- data_ii$e_work*inf_period # Limit other contacts if needed
      #school_c <- (data_ii$e_school)*inf_period # deprecated as merged with work
      other_c <- data_ii$e_other*inf_period

      # Fix NA entries
      if(is.na(home_c)){home_c <- 0}
      if(is.na(work_c)){work_c <- 0}
      #if(is.na(school_c)){school_c <- 0}
      if(is.na(other_c)){other_c <- 0}
      
      # Add max limit on other settings
      scale_other <- 1
      
      if(max_lim_children==T | adult_T==T){ # cgeck if children included
        scale_other <- min(1,(max_contacts*inf_period)/other_c) # scale down based on max other contacts
      }
      
      # Generate basic infections
      home_inf_basic <- rbinom(1,home_c,prob=hh_risk*inf_propn)
      work_inf_basic <- rbinom(1,work_c,prob=cc_risk*inf_propn)
      other_inf_basic <- rbinom(1,other_c,prob=cc_risk*inf_propn)
      rr_basic_ii <- home_inf_basic + work_inf_basic + other_inf_basic
      
      # Generate infections with interventions in place
      inf_ratio_h <- ifelse(scenario_pick=="isolation_outside",inf_ratio,1)
      inf_ratio_w <-  ifelse(wfh_t,0,inf_ratio) # check if reduced work contacts
      inf_ratio_o <-  ifelse(other_t,0,inf_ratio) # check if reduced work contacts
      
      home_infect <- rbinom(1,home_inf_basic,prob=inf_ratio_h)
      work_infect <- rbinom(1,work_inf_basic,prob=inf_ratio_w*extra_red)
      other_infect <- rbinom(1,other_inf_basic,prob=inf_ratio_o*scale_other*extra_red) # scale by maximum
      rr_ii <- home_infect + work_infect + other_infect
      
      # Contact tracing - tally contacts
      home_traced <- rbinom(1,home_c,prob=hh_trace)
      work_traced <- rbinom(1,work_c,prob=ww_trace_all)
      other_traced <- rbinom(1,other_c,prob=other_trace_all*scale_other)
      
      # Infections averted
      home_averted <- rbinom(1,home_infect,prob=hh_trace*trace_adherence)
      work_averted <- rbinom(1,work_infect,prob=ww_trace_all*trace_adherence)
      other_averted <- rbinom(1,other_infect,prob=other_trace_all*trace_adherence)
      
      # Calculate secondary cases averted, adjusting for delay to isolation/tracing
      # This adjusment accounts for only proportion of onward transmission potentially being averted
      
      if(tested_T==T & symp_T==T & do_tracing==T ){
        
        # Total infections eventually averted
        total_averted1 <- home_averted + work_averted + other_averted
        
        # Sample infection times of contacts (relative to onset of infectiousness in primary)
        sample_inf_contacts <- sample(inf_period_ii,total_averted1,replace = T)
        total_averted <- sample_inf_contacts
        
        # Calculate which contacts would not be traced before onset of infectiousness
        # (Delay onset to infectious secondary) <= isolated + test delay
        
        if(scenario_pick=="cell_phone" | scenario_pick=="cell_phone_met_limit"){
            delay_to_trace <- rep(test_delay,total_averted1)
          }else{
            delay_to_trace <- c(rep(test_delay,home_averted), rep(test_delay+ trace_delay,work_averted + other_averted) ) # Instant for phone, X days for manual
        } 
        
        # Identify source of tracing if using both app and manual
        
        if(scenario_pick=="isolation_manual_tracing_met_only_cell" | scenario_pick=="isolation_manual_tracing_met_only_cell_met_limit"){
          delay_to_trace <- rep(test_delay,total_averted1)
          delay_to_trace <- (runif(total_averted1)>tracePhone) # decide if located by phone
          delay_to_trace[delay_to_trace] <- test_delay + trace_delay # set delay based on phone
        }
        
        upper_timing <- inf_period_ii +  delay_to_trace # time of quarantine of contacts
        pick_missed <- ( (sample_inf_contacts + (inf_period-pre_inf) ) <= upper_timing ) # select those quarantined late

        # Scale by proportion of secondary infections while awaiting testing
        delay_to_inf <- (5-pre_inf) # how long from exposure to infectiousness
        
        tally_n <- rep(1,total_averted1)

        contact_inf_duration <- pmin(inf_period, pmax(0,upper_timing-(sample_inf_contacts+delay_to_inf)) ) # days infectious (max = inf_period)
        tally_n[pick_missed] <- (1-(contact_inf_duration[pick_missed])/inf_period) # calculate proportion of onward transmisison pre-quarantine
        
        # Calculate expected secondary transmission averted
        total_averted <- sum(tally_n)
        
      }else{
        total_averted <- 0
      }
  
      rr_reduced <- rr_ii - total_averted
      
      rr_reduced
      
      # Tally up contacts traced
      
      if(scenario_pick=="no_measures" | scenario_pick=="pop_testing" | 
         scenario_pick=="isolation_only" | scenario_pick=="isolation_outside"){
        total_traced <- 0
      }
      
      if(scenario_pick=="hh_quaratine_only"){
        total_traced <- home_traced 
      }
      if(scenario_pick=="hh_work_only"){
        total_traced <- home_traced + work_traced 
      }
      if(scenario_pick=="isolation_manual_tracing_met_only" | scenario_pick=="isolation_manual_tracing_met_limit" | 
         scenario_pick== "isolation_manual_tracing" | scenario_pick== "cell_phone" | 
         scenario_pick=="cell_phone_met_limit" | scenario_pick=="isolation_manual_tracing_met_only_cell" |
         scenario_pick=="isolation_manual_tracing_met_only_cell_met_limit" | scenario_pick=="pop_testing_tracing_met_only"){
        total_traced <- home_traced + work_traced + other_traced 
      }
        
      # store outputs
      store_r <- rbind(store_r,c(rr_basic_ii,rr_reduced,total_traced,inf_period_ii))
    
    } # end individual loop
    
    # Convert
    store_r <- as_tibble(store_r)
    names(store_r) <- c("rr","rr_reduced","total_traced","inf_p")
    
    diff_r <- store_r$rr - store_r$rr_reduced
    
    # Output R and store table
    if(output_r == T){write_csv(store_r,paste0(dir_pick,"rr_out/RR_",scenario_pick,".csv"))}

    store_table_scenario <- rbind(store_table_scenario,c(scenario_pick,mean(store_r$rr),mean(diff_r),
                                                         mean(store_r$rr_reduced),median(store_r$total_traced),mean(store_r$total_traced),quantile(store_r$total_traced,c(0.05,0.95)),
                                                         sum(store_r$rr_reduced>1)/n_run,sum(store_r$rr_reduced>3)/n_run
                                                         ))
    
  } # end scenario loop
  
  # Convert
  store_table_scenario <- as_tibble(store_table_scenario)
  names(store_table_scenario) <- c("scenario","basic","r_diff","r_eff","contacts_med","contacts","traced_90_1","traced_90_2","above1","above5")
  store_table_scenario$r_eff <- as.numeric(store_table_scenario$r_eff)
  store_table_scenario$basic <- as.numeric(store_table_scenario$basic)
  store_table_scenario$above1 <- as.numeric(store_table_scenario$above1)
  store_table_scenario$above5 <- as.numeric(store_table_scenario$above5)
  store_table_scenario$traced_90_1 <- as.numeric(store_table_scenario$traced_90_1)
  store_table_scenario$traced_90_2 <- as.numeric(store_table_scenario$traced_90_2)
  store_table_scenario$contacts <- signif(as.numeric(store_table_scenario$contacts),2)
  store_table_scenario$contacts_med <- as.numeric(store_table_scenario$contacts_med)
  
  max_val <- store_table_scenario$r_eff[1]
  
  pop_UK <- 67820904
  non_risk <- 0.001
  nn_COVID <- round(non_risk*pop_UK)
  
  store_table_scenarioA <- store_table_scenario %>% mutate(reduction_raw = 1-r_eff/basic,
                                                           aboveX = paste0(100* signif(above1,2),"%"), #,100* signif(above5,2),"%"),
                                                           r_effective = signif(max_val*(1-reduction_raw),2), # Normalise for consistency
                                                           reduction = paste0(100* signif(reduction_raw,2),"%"),
                                                           t_contacts = paste0(signif(contacts_med,2)," (",signif(traced_90_1,2),"-",signif(traced_90_2,2),")"),
                                                           total_contacts100 =  signif(p_tested*contacts*baseline_incident_cases,2),
                                                           total_contacts50 =   signif(p_tested*contacts*0.25*baseline_incident_cases,2),
                                                           total_contacts10 =   signif(p_tested*contacts*0.05*baseline_incident_cases,2),
                                                           contactsR1 = signif(p_tested*(nn_COVID)*contacts,3),
                                                           wfh = wfh_prob,
                                                           limit_other = max_low_fix,
                                                           trace_p = trace_prop,
                                                           app_cov,
                                                           p_symptomaticA,
                                                           prob_t_asymp,
                                                           wfh_probC,
                                                           cc_risk,
                                                           test_delay,
                                                           trace_delay,
                                                           max_lim_children,
                                                           p_tested
                                                           )
  

  write_csv(store_table_scenarioA,paste0(dir_pick,"table",hh_risk,"_",cc_risk,"_minother_",max_low_fix,"_wfh_",
                                         wfh_prob,"_trace_",trace_prop,"_symp_",p_symptomaticA,"_app_",app_cov,"_tasymp_",prob_t_asymp,".csv"))


}

# Plot R distribtion -----------------------------------------------------------

plot_R_distribution <- function(dir_pick){
  
  file_names <- list.files(paste0(dir_pick,"rr_out/"))
  n_files <- length(file_names)

  par(mfrow=c(1,2),mar=c(3,3,1,1),mgp=c(2,0.6,0))
  yy_lim <- c(1,5e3)
  
  plot(1,1,log="xy",xlab="secondary cases",ylab="frequency",ylim=yy_lim,xlim=c(0.5,100),yaxt="n",xaxt="n",col="white",pch=19,cex=0.8)
  xticks <- c(0.5,10^seq(0,3,1)); xtick_lab <- c(0,10^seq(0,3,1))
  axis(1, at = xticks,labels = xtick_lab,col = "black") 
  axis(2, at = 10^c(0:3),labels = 10^c(0:3),col = "black") 
  
  for(ii in 9){
    
    store_r <- read_csv(paste0(dir_pick,"rr_out/",file_names[ii]))
  
    dataRR <- table(store_r$rr)
    xx <- names(dataRR); xx[xx==0] <- 0.5; yy <- dataRR
    
    points(xx,yy,col=rgb(0,0,0,0.5),cex=0.8,pch=19)
    
    mean_yy <- mean(store_r$rr)
    
    lines(c(mean_yy,mean_yy),c(0.5,1e4),lty=2)
  
  }
  
  title(LETTERS[1],adj=0)
  
  # Plot distribution of scenarios
  
  label_list <- c("baseline","shorter delay","longer delay","no pre-symptomatic")
  
  plot(-1,-1,ylab="probability",xlab="days from infectious-to-isolation",xlim=c(0,6),ylim=c(0,1))
  lines(c(1,1),c(-1,2),col="grey",lty=2)
  
  store_disnt <- list(c(0,0.25,0.25,0.2,0.3,0),
                      c(0,0.8,0.2,0,0,0),
                      c(0,0,0.25,0.25,0.2,0.3),
                      c(0.25,0.25,0.2,0.3,0,0)
  )

  for(ii in 1:3){
    data_ii <- c(store_disnt[[ii]],0)
    points(c(0:6),data_ii,col=col_def[[ii]],pch=19)
    lines(c(0:6),data_ii,col=col_def[[ii]])
    text(x=2,y=1-0.1*ii,labels=label_list[ii],adj=0,col=col_def[[ii]],cex=0.8)
  }
  
  title(LETTERS[2],adj=0)
  
  #dev.copy(png,paste0(dir_pick,"rr_supp_plot.png"),units="cm",width=20,height=8,res=150)
  dev.copy(pdf,paste0(dir_pick,"Figure_S1.pdf"),width=8,height=3)
  dev.off()
  
  
}

# Plot networks -----------------------------------------------------------

plot_networks <- function(dir_pick){
  
  par(mfrow=c(3,2),mar=c(3,3,1,1),mgp=c(2,0.6,0),las=0)
  
  label_list <- c("home contacts","work or school contacts","other contacts")
  
  # Plot distributions
  for(ii in 1:2){

    if(ii==1){dataH1 <- table(data_user_col_red_u18$e_home)
              dataH2 <- table(data_user_col_red_u18$e_work)
              dataH3 <- table(data_user_col_red_u18$e_other)
    }
    if(ii==2){dataH1 <- table(data_user_col_red_o18$e_home)
              dataH2 <- table(data_user_col_red_o18$e_work)
              dataH3 <- table(data_user_col_red_o18$e_other)
    }


    xx <- names(dataH1) %>% as.numeric(); xx[xx==0] <- 0.5; yy <- dataH1 %>% as.numeric(); #yy <- yy/sum(yy);
    yy_lim <- if(ii==1){c(1,1.4e3)}else{c(1,1.5e4)}  #1.4e4 #ifelse(ii==1,)
    plot(xx,yy,log="xy",xlab="contacts",ylab="participants",ylim=yy_lim,xlim=c(0.5,1000),yaxt="n",xaxt="n",col=col_def_F[[1]],pch=19,cex=0.8)
    
    xx <- names(dataH2) %>% as.numeric(); xx[xx==0] <- 0.5; yy <- dataH2 %>% as.numeric(); #yy <- yy/sum(yy);
    points(xx,yy,col=col_def_F[[2]],cex=0.8,pch=19)
    
    xx <- names(dataH3) %>% as.numeric(); xx[xx==0] <- 0.5; yy <- dataH3 %>% as.numeric(); #yy <- yy/sum(yy);
    points(xx,yy,col=col_def_F[[3]],cex=0.8,pch=19)

    xticks <- c(0.5,10^seq(0,3,1)); xtick_lab <- c(0,10^seq(0,3,1))
    axis(1, at = xticks,labels = xtick_lab,col = "black") 
    axis(2, at = 10^c(0:3),labels = 10^c(0:3),col = "black") 
    
    if(ii==2){
      for(kk in 1:3){
        text(x=20,y=exp(6+kk), labels=label_list[kk],cex=0.9,col=col_def[[kk]],adj=0)
      }
    }
    
    #title(LETTERS[ii],adj=0)

  }
  
  
  # - - - - 
  # Network plot
  for(ii in 1:4){
    #pick_user <- sample(1:n_user,1)
    #data_ii <- data_user_col_red2[pick_user,]
    
    if(ii==1){contact_val <- c(7,38,4) }
    if(ii==2){contact_val <- c(1,0,2) }
    if(ii==3){contact_val <- c(3,1,3) }
    if(ii==4){contact_val <- c(1,20,0) }
    if(ii==5){contact_val <- c(0,1,6) }
    if(ii==6){contact_val <- c(4,8,12) }

    home_c <- contact_val[1]
    work_c <- contact_val[2] # Limit other contacts if needed
    other_c <- contact_val[3]
    
    total_contacts <- home_c + work_c + other_c
  
    # Network plotting
    edge_mat <- matrix(c(rep(1,total_contacts),c(2:(total_contacts+1))),ncol=2,byrow=F)
    g2 <- graph_from_edgelist(el=edge_mat,directed=F) # define directed graph
  
    # Extract MF colours and match so in correct sequence
    pickshape <- c(rep(col_def[[1]],home_c),rep(col_def[[2]],work_c),rep(col_def[[3]],other_c))
    pickcol <- c('black',rep(col_def[[1]],home_c),rep(col_def[[2]],work_c),rep(col_def[[3]],other_c))
    
    # Plot graph
    par(mar=c(1,1,1,1))
    layout_1 <- layout_nicely #layout_in_circle(g2, order = V(g2))
    plot(g2,layout=layout_1,vertex.size=10,vertex.shape='circle',vertex.label=NA,vertex.color=pickcol,
         vertex.label.cex=0.5,vertex.label.family="",edge.arrow.width=2,edge.arrow.size=0.1,edge.color=rgb(0.6,0.6,0.6),edge.width=0.5) 
    
    #title(LETTERS[ii+2],adj=0)
    
  }

  #dev.copy(png,paste0(dir_pick,"plot1.png"),units="cm",width=14,height=20,res=150)
  dev.copy(pdf,paste0(dir_pick,"Figure_1.pdf"),width=6,height=8)
  dev.off()
  
  
}

# Plot contact limit against effectiveness -----------------------------------------------------------

plot_contacts <- function(dir_pick){
  
  file_names <- list.files(paste0(dir_pick,"runs/"))
  n_files <- length(file_names)
  
  par(mfrow=c(2,2),mar=c(3,3,1,1),mgp=c(2,0.6,0),las=1)
  
  # Compile data
  store_dat <- NULL
  for(ii in 1:n_files){
    
    input_ii <- read_csv(paste0(dir_pick,"runs/",file_names[ii]))
    store_dat <- rbind(store_dat,input_ii)
    
  }
  

  # - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Plot by tracing success
  
  label_list <- c("Self-isolation + HQ + school/work tracing",
                  "SI + manual tracing (acquaintance only)",
                  "SI + manual (acquaintance, max 4 other contacts)",
                  "SI + manual tracing (all)",
                  "SI + app-based tracing",
                  "SI + app-based (max 4 other contacts)")
  
  store_dat0 <- store_dat %>% filter(limit_other==4 & wfh==0 & app_cov==0.53 & trace_p!=0.95) %>% arrange(trace_p)
  
  for(ii in 1:4){
    
    if(ii==1){
      xx <- store_dat0 %>% filter(wfh==0 & scenario=="hh_work_only")
      plot(xx$trace_p,xx$r_eff,xlab="proportion successfully traced outside HH",ylab="R",ylim=c(0.5,2.5),xlim=c(0,1),
           col="white",type="l",lwd=2)
      lines(c(0,1e3),c(1,1),col="dark grey")
      #grid(ny = NULL, nx=NA, col = "lightgray")
      kk <- ii+4
    }

    if(ii==2){
      xx <- store_dat0 %>% filter(wfh==0 & scenario=="isolation_manual_tracing_met_only") 
    }
    if(ii==3){
      xx <- store_dat0 %>% filter(wfh==0 & scenario=="isolation_manual_tracing_met_limit") 
    }
    if(ii==4){
      xx <- store_dat0 %>% filter(wfh==0 & scenario=="isolation_manual_tracing") 
    }

    
    lines(xx$trace_p,xx$r_eff,col=col_def[[ii]],lwd=2) 
    
    text(x=0,y=2.5-0.1*(ii-1), labels=label_list[ii],cex=0.7,col=col_def[[ii]],adj=0)
    
  }
  
  
  title(LETTERS[1],adj=0)
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Plot by other contact limit scenarios
  label_list <- c("Self-isolation + app-based tracing",
                  "SI + manual tracing (acquaintance only)",
                  "- - 50% of adults with no work contacts")
  
  store_dat0 <- store_dat %>% filter(trace_p==0.95 & app_cov==0.53 & limit_other!=4) %>% arrange(limit_other)
  
  for(ii in 1:4){
  
    if(ii==1){
      xx <- store_dat0 %>% filter(wfh==0  & scenario=="cell_phone_met_limit") 
      plot(xx$limit_other,xx$r_eff,log="x",xlab="maximum daily other contacts",ylab="R",ylim=c(0.5,2.5),xlim=c(1,100),
           xaxt="n",col="white",type="l",lwd=2)
      #grid(ny = NULL, nx=NA, col = "lightgray")
      lines(c(1e-5,1e3),c(1,1),col="dark grey")
      xticks <- c(10^seq(0,3,1)); xtick_lab <- c(10^seq(0,3,1))
      axis(1, at = xticks,labels = xtick_lab,col = "black") 
      kk <- 7; l_type <- 1
    }
      
    if(ii==2){
      xx <- store_dat0 %>% filter(wfh==0.5 & scenario=="cell_phone_met_limit") 
      kk <- 7; l_type <- 2
    }
    if(ii==3){
      xx <- store_dat0 %>% filter(wfh==0 & scenario=="isolation_manual_tracing_met_limit") 
      kk <- 3; l_type <- 1
    }
    if(ii==4){
      xx <- store_dat0 %>% filter(wfh==0.5 & scenario=="isolation_manual_tracing_met_limit") 
      kk <- 3; l_type <- 2
    }
    
    lines(xx$limit_other,xx$r_eff,col=col_def[[kk]],lwd=2,lty=l_type) 
  }
  
  text(x=1,y=2.5-0.1*(1-1), labels=label_list[1],cex=0.7,col=col_def[[7]],adj=0)
  text(x=1,y=2.5-0.1*(2-1), labels=label_list[2],cex=0.7,col=col_def[[3]],adj=0)
  text(x=1,y=2.5-0.1*(3-1), labels=label_list[3],cex=0.7,col="black",adj=0)
  
  
  title(LETTERS[2],adj=0)
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Plot by WFH scenarios
  
  label_list <- c("SI + manual tracing (acquaintance only)",
                  "SI + manual tracing (all)",
                  "SI + app-based tracing")
  
  store_dat0 <- store_dat %>% filter(limit_other==4 & wfh_probC==0 & app_cov==0.53 & trace_p==0.95) %>% arrange(wfh)
  
  for(ii in 1:3){
    
    if(ii==1){
      xx <- store_dat0 %>% filter(scenario=="isolation_manual_tracing_met_only")
      plot(xx$wfh,xx$r_eff,xlab="proportion with no work contacts",ylab="R",ylim=c(0.5,2.5),xlim=c(0,0.6),col="white",type="l",lwd=2)
      lines(c(0,1e3),c(1,1),col="dark grey")
      #grid(ny = NULL, nx=NA, col = "lightgray")
      kk <- 2
    }
    # if(ii==2){
    #   xx <- store_dat %>% filter(wfh==0 & limit_other==6 & scenario=="hh_quaratine_only") 
    # }

    if(ii==2){
      xx <- store_dat0 %>% filter(scenario=="isolation_manual_tracing") 
      kk <- 4
    }
    if(ii==3){
      xx <- store_dat0 %>% filter(scenario=="cell_phone") 
      kk <- 7
    }
    

    lines(xx$wfh,xx$r_eff,col=col_def[[kk]],lwd=2) 
    
    text(x=0,y=2.5-0.1*(ii-1), labels=label_list[ii],cex=0.7,col=col_def[[kk]],adj=0)
    
    if(ii==3){
    text(x=0,y=2.5-0.1*(ii), labels="- - with school closure",cex=0.7,col="black",adj=0)
    }
    
  }
  
  # With school closure
  store_dat0 <- store_dat %>% filter(limit_other==4 & wfh_probC==1 & app_cov==0.53 & trace_p==0.95) %>% arrange(wfh)
  
  for(ii in 1:3){
    
    if(ii==1){
      xx <- store_dat0 %>% filter(scenario=="isolation_manual_tracing_met_only")
      kk <- 2
    }
    
    if(ii==2){
      xx <- store_dat0 %>% filter(scenario=="isolation_manual_tracing") 
      kk <- 4
    }
    if(ii==3){
      xx <- store_dat0 %>% filter(scenario=="cell_phone") 
      kk <- 7
    }

    lines(xx$wfh,xx$r_eff,col=col_def[[kk]],lwd=2,lty=2) 

  }
  
  
  title(LETTERS[3],adj=0)
  
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Plot by app coverage
  
  label_list <- c("SI + app-based tracing",
                  "- - with max 4 other contacts")
  
  store_dat0 <- store_dat %>% filter(limit_other==4 & wfh==0 & app_cov!=0.53 & trace_p==0.95) %>% arrange(app_cov)
  
  for(ii in 1:2){
    
    if(ii==1){
      xx <- store_dat0 %>% filter(scenario=="cell_phone")
      plot(xx$app_cov,xx$r_eff,xlab="app coverage",ylab="R",ylim=c(0.5,2.5),xlim=c(0,1),col="white",type="l",lwd=2)
      lines(c(0,1e3),c(1,1),col="dark grey")
      #grid(ny = NULL, nx=NA, col = "lightgray")
      kk <- 7; l_type=1
    }
    # if(ii==2){
    #   xx <- store_dat %>% filter(wfh==0 & limit_other==6 & scenario=="hh_quaratine_only") 
    # }
    
    if(ii==2){
      xx <- store_dat0 %>% filter(scenario=="cell_phone_met_limit") 
      kk <- 7; l_type=2
    }
    
    lines(xx$app_cov,xx$r_eff,col=col_def[[kk]],lwd=2,lty=l_type) 
    
    if(ii==1){
      text(x=0,y=2.5-0.1*(ii-1), labels=label_list[ii],cex=0.7,col=col_def[[kk]],adj=0)
    }
    if(ii==2){
      text(x=0,y=2.5-0.1*(ii-1), labels=label_list[ii],cex=0.7,col="black",adj=0)
    }

    
  }
  
  
  title(LETTERS[4],adj=0)
  
  #dev.copy(png,paste0(dir_pick,"plot2.png"),units="cm",width=20,height=16,res=150)
  dev.copy(pdf,paste0(dir_pick,"Figure_2.pdf"),width=8,height=6)
  dev.off()
  

  
}

# Table of contact range -----------------------------------------------------------


table_outputs_contact_range <- function(dir_pick){
  
  file_names <- list.files(paste0(dir_pick,"runs/"))
  n_files <- length(file_names)
  
  # Compile data
  store_dat <- NULL
  for(ii in 1:n_files){
    
    input_ii <- read_csv(paste0(dir_pick,"runs/",file_names[ii]))
    store_dat <- rbind(store_dat,input_ii)
    
  }
  
  store_dat$wfh <- 1-store_dat$wfh # convert to % active
  store_dat <- store_dat %>% mutate(c_safe = 1-cc_risk/0.06)
  
  # Extract ranges
  
  trace_range <- unique(store_dat$trace_p)
  c_safe_range <- unique(store_dat$c_safe) %>% sort()
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Plot by WFH scenarios
  
  col_def <- list("red","orange","cyan","grey")
  
  label_list <- c("- Self-isolation and HH quarantine",
                  "- - SI + HHQ + 35% contact tracing")

  par(mfrow=c(length(c_safe_range),length(trace_range)),mar=c(3,3,1,1),mgp=c(2,0.6,0),las=1)
  
  store_dat0 <- store_dat %>% arrange(wfh)

  for(ii in 1:length(c_safe_range)){
    for(kk in 1:length(trace_range)){
      
      for(jj in 1:2){
          
          if(jj==1){scenario_pick <- "hh_quaratine_only"}
          if(jj==2){scenario_pick <- "isolation_manual_tracing"}
          
          xx <- store_dat0 %>% filter(scenario==scenario_pick) %>% filter(c_safe==c_safe_range[ii] & trace_p==trace_range[kk])
          
          if(jj==1){
            plot(xx$wfh,xx$r_eff,xlab="% active work and other contacts",
                 main = paste0("C-S: ",100*c_safe_range[ii],"% , trace:", 100*trace_range[kk],"%"),
                 ylab="R",ylim=c(0,2),xlim=c(0.2,1),col="white",type="l",lwd=2)
            lines(c(0.2,1),c(1,1),col="dark grey",lty=2)
            #grid(ny = NULL, nx=NA, col = "lightgray")
            kk <- 2
          }
          if(ii==1){
            text(x=0.2,y=2-0.1*jj, labels=label_list[jj],cex=1,col="black",adj=0)
          }
      
          lines(xx$wfh,xx$r_eff,col=col_def[[ii]],lwd=2,lty=jj) 
          
          if(jj==1){
            text(x=0.8,0.5-0.1*ii,labels=paste0(round(40),"% school contacts"),col=col_def[[ii]],adj=0)
          }
      }
    }
  }
  

  
  dev.copy(pdf,paste0(dir_pick,"Figure_contacts.pdf"),width=6,height=5)
  dev.off()
  
  # Output csv
  output_plot <- store_dat %>% filter(scenario=="hh_quaratine_only" | scenario=="isolation_manual_tracing") %>%
                           mutate(prop_active_contacts=signif(wfh,3),prop_covid_safe=signif(1-cc_risk/0.06),3) %>%
                           select(scenario,prop_active_contacts,prop_covid_safe,r_eff)
   
  write_csv(output_plot,paste0(dir_pick,"output_estimates.csv"))
  
  
  
}

# Table of delay range -----------------------------------------------------------


table_outputs_delay_range <- function(dir_pick){
  
  file_names <- list.files(paste0(dir_pick,"runs_delay/"))
  n_files <- length(file_names)
  
  # Compile data
  store_dat <- NULL
  for(ii in 1:n_files){
    
    input_ii <- read_csv(paste0(dir_pick,"runs_delay/",file_names[ii]))
    store_dat <- rbind(store_dat,input_ii)
    
  }

  # Extract ranges
  trace_range <- unique(store_dat$trace_p)
  test_delay_range <- unique(store_dat$test_delay) %>% sort()
  trace_delay_range <- unique(store_dat$trace_delay) %>% sort()
  
  ptest_delay_range <- unique(store_dat$p_tested) %>% sort()
  t_aymp_delay_range <- unique(store_dat$prob_t_asymp) %>% sort()
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Plot by scenarios
  
  col_def <- list("blue","orange","grey")
  label_list <- c("Delay test-to-trace: 1d",
                  "Delay test-to-trace: 2d",
                  "Delay test-to-trace: 3d")
  
  

  # - - -
  # Loop over options
  
  for(pp1 in 1:length(ptest_delay_range)){
  for(pp2 in 1:length(t_aymp_delay_range)){
    
    store_dat0 <- store_dat %>% filter(p_tested==ptest_delay_range[pp1] & prob_t_asymp==t_aymp_delay_range[pp2]) %>% arrange(trace_p)
    
    par(mfrow=c(1,5),mar=c(3,3,1,1),mgp=c(2,0.6,0),las=1)
    
  
    for(ii in 1:length(test_delay_range)){
      
      for(kk in 1:length(trace_delay_range)){
        
        for(jj in 2){
          
          if(jj==1){scenario_pick <- "hh_quaratine_only"}
          if(jj==2){scenario_pick <- "isolation_manual_tracing"}
          
          xx <- store_dat0 %>% filter(scenario==scenario_pick) %>% filter(test_delay==test_delay_range[ii] & trace_delay==trace_delay_range[kk])
          
          if(jj==2 & kk==1){
            plot(xx$wfh,xx$r_eff,xlab="% contacts traced and quarantined",
                 main = paste0("Delay isolation-to-test: ",test_delay_range[ii]," days"),
                 ylab="reduction in transmission (%)",ylim=c(0,100),xlim=c(0,1),col="white",type="l",lwd=2)
            #lines(c(0.2,1),c(1,1),col="dark grey",lty=2)
            grid(ny = NULL, nx=NA, col = "lightgray")
          }
          
          if(jj==2){
            #text(x=0.2,y=100*(1-0.1*jj), labels=label_list[jj],cex=1,col="black",adj=0)
            text(x=0.2,y=80*(1-0.05*kk), labels=label_list[kk],cex=1,col=col_def[[kk]],adj=0)
          }
          
          lines(xx$trace_p,100*xx$reduction_raw,col=col_def[[kk]],lwd=2,lty=1) 
          
          # if(jj==1){
          #   text(x=0.8,0.5-0.1*ii,labels=paste0(round(40),"% school contacts"),col=col_def[[ii]],adj=0)
          # }
        } # end scenarios
      } # end trace delay
      
    } # end test delay
  
    dev.copy(pdf,paste0(dir_pick,"Figure_delay_ptest_",pp1,"_tasymp",pp2,".pdf"),width=12,height=4)
    dev.off()
  
  }
  }
  
  # Output csv
  # output_plot <- store_dat %>% filter(scenario=="hh_quaratine_only" | scenario=="isolation_manual_tracing") %>%
  #   mutate(prop_active_contacts=signif(wfh,3),prop_covid_safe=signif(1-cc_risk/0.06),3) %>%
  #   select(scenario,prop_active_contacts,prop_covid_safe,r_eff)
  # 
  # write_csv(output_plot,paste0(dir_pick,"output_estimates.csv"))
  
  
  output_plot <- store_dat %>% select(p_tested,prob_t_asymp,test_delay,trace_delay,trace_p,reduction_raw)
  output_plot <- output_plot %>% rename(proportion_symptomatic_isolate_and_tested = p_tested,
                                        relative_transmission_from_asymptomatics = prob_t_asymp,
                                        delay_isolation_to_test = test_delay,
                                        delay_test_to_trace_and_quarantine = trace_delay,
                                        proportion_contacts_traced_and_quarantine = trace_p,
                                        estimated_transmission_reduction = reduction_raw)
  
  write_csv(output_plot,paste0(dir_pick,"output_values_plot.csv"))
  
  
  
}


# Plot symptomatic sensitivity -----------------------------------------------------------

plot_symptom_reduction <- function(dir_pick){
  
  file_names <- list.files(paste0(dir_pick,"runs2/"))
  n_files <- length(file_names)
  
  par(mfrow=c(1,2),mar=c(3,3,1,1),mgp=c(2,0.6,0),las=1)
  
  # Compile data
  store_dat <- NULL
  for(ii in 1:n_files){
    
    input_ii <- read_csv(paste0(out_dir,"runs2/",file_names[ii]))
    store_dat <- rbind(store_dat,input_ii)
    
  }
  
  # Swap ordering
  store_dat$reduction_raw <- 1-store_dat$reduction_raw
  

  # - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Plot different symptomatic assumptions
  label_list <- c("SI + manual tracing (acquaintance only)",
                  "SI + manual (acquaintance, max 4 other contacts)",
                  "SI + manual tracing (all)",
                  "SI + app-based tracing",
                  "SI + app-based (max 4 other contacts)")
  
  store_dat0 <- store_dat %>% filter(trace_p==0.95 & app_cov==0.53 & prob_t_asymp==0.5 & limit_other==4) %>% arrange(p_symptomaticA)
  
  for(ii in 1:5){
    
    if(ii==1){
      xx <- store_dat0 %>% filter(wfh==0  & scenario=="isolation_manual_tracing_met_only") 
      plot(xx$p_symptomaticA,xx$reduction_raw,xlab="proportion adults symptomatic",ylab="relative R",ylim=c(0,1),xlim=c(0.2,0.8),col="white",type="l",lwd=2)
      lines(c(0.6,0.6),c(-1,2),lty=2,col="grey")
      #lines(c(1e-5,1e3),c(1,1),col="dark grey")
      #xticks <- c(10^seq(0,3,1)); xtick_lab <- c(10^seq(0,3,1))
      #axis(1, at = xticks,labels = xtick_lab,col = "black") 
      kk <- 2
    }
    
    if(ii==2){
      xx <- store_dat0 %>% filter(wfh==0 & scenario=="isolation_manual_tracing_met_limit") 
      kk <- 3; l_type <- 1
    }
    if(ii==3){
      xx <- store_dat0 %>% filter(wfh==0 & scenario=="isolation_manual_tracing") 
      kk <- 4; l_type <- 1
    }
    if(ii==4){
      xx <- store_dat0 %>% filter(wfh==0 & scenario=="cell_phone") 
      kk <- 7; l_type <- 1
    }
    if(ii==5){
      xx <- store_dat0 %>% filter(wfh==0 & scenario=="cell_phone_met_limit") 
      kk <- 5; l_type <- 1
    }
    
    lines(xx$p_symptomaticA,xx$reduction_raw,col=col_def[[kk]],lwd=2) 

    text(x=0.2,y=1-0.05*(ii-1), labels=label_list[ii],cex=0.7,col=col_def[[kk]],adj=0)
    
  }
  
  title(LETTERS[1],adj=0)
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Plot different asymptomatic transmission assumptions
  label_list <- c("SI + manual tracing (acquaintance only)",
                  "SI + manual (acquaintance, max 4 other contacts)",
                  "SI + manual tracing (all)",
                  "SI + app-based tracing",
                  "SI + app-based (max 4 other contacts)")
  
  store_dat0 <- store_dat %>% filter(trace_p==0.95 & app_cov==0.53 & p_symptomaticA==0.7 & limit_other==4) %>% arrange(prob_t_asymp)
  
  for(ii in 1:5){
    
    if(ii==1){
      xx <- store_dat0 %>% filter(wfh==0  & scenario=="isolation_manual_tracing_met_only") 
      plot(xx$prob_t_asymp,xx$reduction_raw,xlab="relative asymptomatic transmission",ylab="relative R",ylim=c(0,1),
           xlim=c(0,1),col="white",type="l",lwd=2)
      lines(c(0.5,0.5),c(-1,2),lty=2,col="grey")
      #lines(c(1e-5,1e3),c(1,1),col="dark grey")
      #xticks <- c(10^seq(0,3,1)); xtick_lab <- c(10^seq(0,3,1))
      #axis(1, at = xticks,labels = xtick_lab,col = "black") 
      kk <- 2
    }
    
    if(ii==2){
      xx <- store_dat0 %>% filter(wfh==0 & scenario=="isolation_manual_tracing_met_limit") 
      kk <- 3
    }
    if(ii==3){
      xx <- store_dat0 %>% filter(wfh==0 & scenario=="isolation_manual_tracing") 
      kk <- 4
    }
    if(ii==4){
      xx <- store_dat0 %>% filter(wfh==0 & scenario=="cell_phone") 
      kk <- 7
    }
    if(ii==5){
      xx <- store_dat0 %>% filter(wfh==0 & scenario=="cell_phone_met_limit") 
      kk <- 5
    }
    
    lines(xx$prob_t_asymp,xx$reduction_raw,col=col_def[[kk]],lwd=2) 

    text(x=0.05,y=1-0.05*(ii-1), labels=label_list[ii],cex=0.7,col=col_def[[kk]],adj=0)
    
  }
  
  title(LETTERS[2],adj=0)
  
  #dev.copy(png,paste0(dir_pick,"plot_symp.png"),units="cm",width=20,height=8,res=150)
  dev.copy(pdf,paste0(dir_pick,"Figure_S2.pdf"),width=8,height=3)
  dev.off()
  
}

# Compile table -----------------------------------------------------------

table_outputs <- function(dir_pick){
  
  # Compile delay sensitivity analysis
  
  matchLL <- c("no_measures","isolation_only","isolation_outside",
                     "hh_quaratine_only","hh_work_only",
                     "isolation_manual_tracing_met_only",
                     "isolation_manual_tracing",
                     "cell_phone","isolation_manual_tracing_met_only_cell",
                     "isolation_manual_tracing_met_limit",
                     "isolation_manual_tracing_met_only_cell_met_limit",
                     "pop_testing")
  
  input_base <- read_csv(paste0(out_dir,"table0.2_minother_4_wfh_0_trace_0.95_symp_0.7_app_0.53_tasymp_0.5.csv"))
  input_longer <- read_csv(paste0(out_dir,"sensitivity/","late_detection_table0.2_minother_4_wfh_0_trace_0.95_symp_0.7_app_0.53_tasymp_0.5.csv"))
  input_fast <- read_csv(paste0(out_dir,"sensitivity/","fast_isolation_table0.2_minother_4_wfh_0_trace_0.95_symp_0.7_app_0.53_tasymp_0.5.csv"))
  
  
  out_tab <- cbind(input_base$scenario,input_base$reduction,input_base$t_contacts,
        input_fast$reduction,input_fast$t_contacts,
        input_longer$reduction,input_longer$t_contacts
        )
  
  # Edit ordering
  out_tab <-  out_tab[match(matchLL, out_tab[,1]),] 
  
  write_csv(as_tibble(out_tab),paste0(out_dir,"sensitivity/compile_table_duration.csv"))
  
  # Compile SAR sensitivity analysis
  
  input_base <- read_csv(paste0(out_dir,"table0.2_minother_4_wfh_0_trace_0.95_symp_0.7_app_0.53_tasymp_0.5.csv"))
  input_cc_7 <- read_csv(paste0(out_dir,"sensitivity/","CC_SAR_higher_table0.2_minother_4_wfh_0_trace_0.95_symp_0.7_app_0.53_tasymp_0.5.csv"))
  input_hh_40 <- read_csv(paste0(out_dir,"sensitivity/","HH_SAR_higher_table0.5_minother_4_wfh_0_trace_0.95_symp_0.7_app_0.53_tasymp_0.5.csv"))

  
  out_tab <- cbind(input_base$scenario,input_base$reduction,input_base$t_contacts,
                   input_cc_7$reduction,input_cc_7$t_contacts,
                   input_hh_40$reduction,input_hh_40$t_contacts,
                   input_cc_7$r_eff,input_hh_40$r_eff
  )
  
  # Edit ordering
  out_tab <-  out_tab[match(matchLL, out_tab[,1]),] 
  
  write_csv(as_tibble(out_tab),paste0(out_dir,"sensitivity/compile_table_SAR.csv"))
  
  
}



# Compile table (V0) -----------------------------------------------------------


table_outputs_1 <- function(dir_pick){
  
  input_base <- read_csv(paste0(out_dir,"table0.2_minother_4_wfh_0_trace_0.95_symp_0.7_app_0.53_tasymp_0.5.csv"))
  
  input_base2 <- input_base[match(c("no_measures","isolation_only","isolation_outside",
                                                   "hh_quaratine_only","hh_work_only",
                     "isolation_manual_tracing_met_only",
                     "isolation_manual_tracing",
                     "cell_phone","isolation_manual_tracing_met_only_cell",
                     "isolation_manual_tracing_met_limit",
                     "isolation_manual_tracing_met_only_cell_met_limit",
                     "pop_testing"),input_base$scenario
                      ),] %>% dplyr::select(scenario,aboveX,r_effective,reduction,
                                     t_contacts,total_contacts100,total_contacts50,total_contacts10)
  
  write_csv(input_base2,paste0(out_dir,"compile_table_base.csv"))
  
  
}
  
# Plot contact limits -----------------------------------------------------------

plot_contacts <- function(dir_pick){
  
  dir_pick <- out_dir
  
  file_names <- list.files(paste0(dir_pick,"runs_contacts/"))
  n_files <- length(file_names)
  
  par(mfrow=c(1,1),mar=c(3,3,1,1),mgp=c(2,0.6,0),las=1)
  
  # Compile data
  store_dat <- NULL
  for(ii in 1:n_files){
    
    input_ii <- read_csv(paste0(dir_pick,"runs_contacts/",file_names[ii]))
    store_dat <- rbind(store_dat,input_ii)
    
  }
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Plot by other contact limit scenarios
  # label_list <- c("SI + manual tracing (acquaintance only)",
  #                 "SI + manual tracing (acquaintance only) + limit on contacts")
  # 
  label_list <- c("Limit for everyone",
                  "Limit for adults only")
  
  
  store_dat0 <- store_dat %>% filter(wfh==0.25) %>% arrange(limit_other)
  
  for(ii in 1:2){
    
    if(ii==1){
      xx <- store_dat0 %>% filter(max_lim_children==T  & scenario=="isolation_manual_tracing_met_limit") 
      plot(xx$limit_other,xx$r_eff,log="x",xlab="maximum daily other contacts",ylab="R",ylim=c(0.5,2),xlim=c(1,100),
           xaxt="n",col="white",type="l",lwd=2)
      #grid(ny = NULL, nx=NA, col = "lightgray")
      lines(c(1e-5,1e3),c(1,1),col="dark grey")
      xticks <- c(10^seq(0,3,1)); xtick_lab <- c(10^seq(0,3,1))
      axis(1, at = xticks,labels = xtick_lab,col = "black") 
      kk <- 7; l_type <- 1
    }
    
    if(ii==2){
      xx <- store_dat0 %>% filter(max_lim_children==F  & scenario=="isolation_manual_tracing_met_limit") 
      kk <- 3; l_type <- 1
    }
    lines(xx$limit_other,xx$r_eff,col=col_def[[kk]],lwd=2,lty=l_type)

  }
  

  lines(c(5,5),c(0,3),col="grey",lwd=2,lty=2) 
  text(x=6,y=0.6, labels="Rule of 6",cex=0.7,col="grey",adj=0)
  
  text(x=1,y=2-0.1*(1-1), labels=label_list[1],cex=0.7,col=col_def[[7]],adj=0)
  text(x=1,y=2-0.1*(2-1), labels=label_list[2],cex=0.7,col=col_def[[3]],adj=0)
  

  
  dev.copy(png,paste0(dir_pick,"plot_contacts.png"),units="cm",width=12,height=12,res=150)
  #dev.copy(pdf,paste0(dir_pick,"Figure_contacts.pdf"),width=8,height=6)
  dev.off()
  
  
  
}
