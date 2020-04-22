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
                            wfh_prob = 0, # Probability people are working from home
                            range_n = NULL, # Pick specific scenarios to run
                            trace_prop = 0.95, # Proportion of contacts traced
                            n_run = 5e3, # Number of simualtions
                            app_cov = 0.53, # App coverage
                            prob_symp = 0.6, # Proportion symptomatic
                            prob_t_asymp = 0.5, # Proportion symptomatic
                            isolate_distn = c(0,0.25,0.25,0.2,0.3,0), # distribution of time to isolate (1st day presymp)
                            dir_pick = "", # Output directory
                            pt_extra = 0, # Optional extra transmission intervention
                            pt_extra_reduce = 0, # Reduction from extra intervention
                            output_r = F,
                            hh_risk = 0.2, # HH risk
                            cc_risk = 0.065 # Outside HH contact risk
                            ){
  
  # DEBUG
  # max_low_fix = 4; wfh_prob = 0; range_n = NULL; trace_prop = 0.95; n_run = 5e3; app_cov = 0.53; prob_symp = 0.6; prob_t_asymp = 0.5; isolate_distn = c(0,0.25,0.25,0.2,0.3,0); dir_pick = ""; pt_extra = 0.95; pt_extra_reduce = 0; output_r = F
  
  # - - - - - - - - - - - - - - - - - - - - 
  # Define parameters across scenarios (note: some will be redefined for specific scenarios)

  # Demographic parameters
  under_18_prob <- 0.21 # Probability under 18
  
  # Transmission and baseline risk
  baseline_incident_cases <- 20e4 # baseline incidence symptomatic cases
  
  # Symptomatic and proportion getting tested
  trace_adherence <- 0.9 # Adherence to testing/quarantine
  p_tested <- trace_adherence # Proportion who get tested
  time_isolate <- isolate_distn # Distribution over symptomatic period
  p_symptomatic <- prob_symp
  transmission_asymp <- prob_t_asymp
  phone_coverage <- 1 # App coverage in non-app scenarios
  p_pop_test <- 0.05 # Proportion mass tested (5% per week)
  inf_period <- 5 # Infectious period
  
  # Define default scenarios
  scenario_list <- c("no_measures","isolation_only","hh_quaratine_only","hh_work_only",
                     "isolation_manual_tracing_met_only","isolation_manual_tracing_met_limit",
                     "isolation_manual_tracing","cell_phone","cell_phone_met_limit",
                     "pop_testing","pt_extra")
  
  if(is.null(range_n)){nn_choose <- 1:length(scenario_list)}else{nn_choose <- range_n}
  
  store_table_scenario <- NULL
  
  # - - - - - 
  # RUN MODEL 
  # Iterate over scenarios
  
  for(kk in nn_choose){
    
    
    # Proportion of contacts met before
    met_before_w <- 0.79 # At work. At school = 90%, which is defined in function later on
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
    }
  
    if(scenario_pick=="isolation_only"){
      do_isolation <- T
      do_tracing <- F
    }
    
    if(scenario_pick=="hh_quaratine_only"){
      do_isolation <- T
      do_tracing <- T
      ww_trace <- 0 # Tracing at work
      other_trace <- 0 # Tracing others
    }
    
    if(scenario_pick=="hh_work_only"){
      do_isolation <- T
      do_tracing <- T
      met_before_w <- 1 # Met before
      met_before_o <- 1 # Met before
      other_trace <- 0 # Tracing others
    }
    
    if(scenario_pick=="isolation_manual_tracing_met_limit"){
      do_isolation <- T
      do_tracing <- T
      max_contacts <- max_low_fix
    }

    if(scenario_pick=="isolation_manual_tracing_met_only"){
      do_isolation <- T
      do_tracing <- T
    }
    
    if(scenario_pick=="isolation_manual_tracing"){
      do_isolation <- T
      do_tracing <- T
      met_before_w <- 1 # Met before
      met_before_o <- 1 # Met before
    }
    
    if(scenario_pick=="cell_phone"){
      do_isolation <- T
      do_tracing <- T
      met_before_w <- 1 # Met before
      met_before_o <- 1 # Met before
      phone_coverage <- app_cov
    }
    
    if(scenario_pick=="cell_phone_met_limit"){
      do_isolation <- T
      do_tracing <- T
      met_before_w <- 1 # Met before
      met_before_o <- 1 # Met before
      phone_coverage <- app_cov
      max_contacts <- max_low_fix
    }

    if(scenario_pick=="pop_testing"){
      do_isolation <- F
      do_tracing <- F
    }
    
    if(scenario_pick=="pt_extra"){
      do_isolation <- T
      do_tracing <- T
      ww_trace <- 0 # Tracing at work
      other_trace <- 0 # Tracing others
    }
    
    
    # - - -
    # Iterate over users
    store_r <- NULL
    
    for(ii in 1:n_run){

      # Decide if child or adult
      if(runif(1)<under_18_prob){
        pick_user <- sample(1:n_user_u18,1)
        data_ii <- data_user_col_red_u18[pick_user,]
        met_before_w <- 0.9
        wfh_t <- F # working from home?
      }else{
        pick_user <- sample(1:n_user_o18,1)
        data_ii <- data_user_col_red_o18[pick_user,]
        wfh_t <- runif(1)< wfh_prob # working from home?
      }

      # Simulate infectious period
      # Decide if symptomatic & tested
      phone_T <- runif(1)<phone_coverage # has phone app?
      symp_T <- runif(1)<p_symptomatic # syptomatic?
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
      
      # Check if contacts phone traced (in cell phone scenario):
      if(scenario_pick=="cell_phone" | scenario_pick=="cell_phone_met_limit"){
        if(phone_T==T){
          ww_trace <- phone_coverage # Tracing at work
          other_trace <- phone_coverage # Tracing others
        }else{
          ww_trace <- 0 # Tracing at work
          other_trace <- 0 # Tracing others
        }
      }
      
      # Extra transmission reduction
      if(scenario_pick=="pt_extra" & runif(1)<pt_extra){
        tested_T <- T
        extra_red <- (1-pt_extra_reduce) # Multiple by both
      }

      # Proportion infectious
      inf_ratio <- inf_period_ii/inf_period

      # Tally contacts
      home_c <- data_ii$e_home
      work_c <- data_ii$e_work*inf_period # Limit other contacts if needed
      #school_c <- (data_ii$e_school)*inf_period # merged with work
      other_c <- data_ii$e_other*inf_period

      # Fix NA entries
      if(is.na(home_c)){home_c <- 0}
      if(is.na(work_c)){work_c <- 0}
      #if(is.na(school_c)){school_c <- 0}
      if(is.na(other_c)){other_c <- 0}
      
      # Collate work/school
      work_c <- work_c    
      scale_other <- min(1,(max_contacts*inf_period)/other_c) # scale down based on max other contacts

      # Generate basic infections
      home_inf_basic <- rbinom(1,home_c,prob=hh_risk*inf_propn)
      work_inf_basic <- rbinom(1,work_c,prob=cc_risk*inf_propn)
      other_inf_basic <- rbinom(1,other_c,prob=cc_risk*inf_propn)
      rr_basic_ii <- home_inf_basic + work_inf_basic + other_inf_basic
      
      # Generate infections
      inf_ratio_w <- ifelse(wfh_t==T,0,inf_ratio) # check if working from home
      
      home_infect <- rbinom(1,home_inf_basic,prob=inf_ratio)
      work_infect <- rbinom(1,work_inf_basic,prob=inf_ratio_w*extra_red)
      other_infect <- rbinom(1,other_inf_basic,prob=inf_ratio*scale_other*extra_red) # scale by maximum
      rr_ii <- home_infect + work_infect + other_infect
      
      # Contact tracing - tally contacts
      home_traced <- rbinom(1,home_c,prob=hh_trace)
      work_traced <- rbinom(1,work_c,prob=ww_trace*met_before_w)
      other_traced <- rbinom(1,other_c,prob=ww_trace*met_before_o*scale_other)
      
      # Infections averted
      home_averted <- rbinom(1,home_infect,prob=hh_trace*trace_adherence)
      work_averted <- rbinom(1,work_infect,prob=ww_trace*met_before_w*trace_adherence)
      other_averted <- rbinom(1,other_infect,prob=met_before_o*other_trace*trace_adherence)
      
      if(tested_T==T & symp_T==T & do_tracing==T ){
        total_averted <- home_averted+work_averted+other_averted
      }else{
        total_averted <- 0
      }
  
      rr_reduced <- rr_ii - total_averted
      
      # Trace contacts-of-contacts (i.e. people who were infected before detection)
      
      if(scenario_pick=="no_measures" | scenario_pick=="pop_testing" | scenario_pick=="isolation_only"){
        total_traced <- 1
      }
      
      if(scenario_pick=="hh_quaratine_only"){
        total_traced <- home_traced + 1
      }
      if(scenario_pick=="hh_work_only"){
        total_traced <- home_traced + work_traced + 1
      }
      if(scenario_pick=="isolation_manual_tracing_met_only" | scenario_pick=="isolation_manual_tracing_met_limit" | 
         scenario_pick== "isolation_manual_tracing" | scenario_pick== "cell_phone" | 
         scenario_pick=="cell_phone_met_limit"){
        total_traced <- home_traced + work_traced + other_traced + 1
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
                                                         mean(store_r$rr_reduced),mean(store_r$total_traced),quantile(store_r$total_traced,c(0.05,0.95)),
                                                         sum(store_r$rr_reduced>1)/n_run,sum(store_r$rr_reduced>3)/n_run
                                                         ))
    
  } # end scenario loop
  
  # Convert
  store_table_scenario <- as_tibble(store_table_scenario)
  names(store_table_scenario) <- c("scenario","basic","r_diff","r_eff","contacts","traced_90_1","traced_90_2","above1","above5")
  store_table_scenario$r_eff <- as.numeric(store_table_scenario$r_eff)
  store_table_scenario$basic <- as.numeric(store_table_scenario$basic)
  store_table_scenario$above1 <- as.numeric(store_table_scenario$above1)
  store_table_scenario$above5 <- as.numeric(store_table_scenario$above5)
  store_table_scenario$contacts <- signif(as.numeric(store_table_scenario$contacts),2)
  
  
  max_val <- store_table_scenario$r_eff[1]
  
  store_table_scenarioA <- store_table_scenario %>% mutate(reduction_raw = 1-r_eff/basic,
                                                           r_effective = signif(max_val*(1-reduction_raw),3), # Normalise for consistency
                                                           aboveX = paste0(100* signif(above1,2),"%, ",100* signif(above5,2),"%"),
                                                           reduction = paste0(100* signif(reduction_raw,2),"%"),
                                                           t_contacts = contacts,
                                                           total_contacts100 =  signif(1e-3*p_tested*contacts*baseline_incident_cases,2),
                                                           total_contacts50 =   signif(1e-3*p_tested*contacts*0.5*baseline_incident_cases,2),
                                                           total_contacts10 =   signif(1e-3*p_tested*contacts*0.25*baseline_incident_cases,2),
                                                           wfh = wfh_prob,
                                                           limit_other = max_low_fix,
                                                           trace_p = trace_prop,
                                                           app_cov,
                                                           prob_symp,
                                                           prob_t_asymp
                                                           )
  

  write_csv(store_table_scenarioA,paste0(dir_pick,"table",hh_risk,"_minother_",max_low_fix,"_wfh_",
                                         wfh_prob,"_trace_",trace_prop,"_symp_",prob_symp,"_app_",app_cov,"_tasymp_",prob_t_asymp,".csv"))


}

# Plot R distribtion -----------------------------------------------------------

plot_R_distribution <- function(dir_pick){
  
  file_names <- list.files(paste0(dir_pick,"rr_out/"))
  n_files <- length(file_names)

  par(mfrow=c(1,1),mar=c(3,3,1,1),mgp=c(2,0.6,0))
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
  
  
  dev.copy(png,paste0(dir_pick,"rr_plot.png"),units="cm",width=10,height=8,res=150)
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

  dev.copy(png,paste0(dir_pick,"plot1.png"),units="cm",width=14,height=20,res=150)
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
  
  label_list <- c("Self-isolation + HQ + work tracing",
                  "SI + manual tracing (familiar only)",
                  "SI + manual (familiar, max 4 other contacts)",
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
  label_list <- c("— Self-isolation + App-based tracing, 0% WFH",
                  "- - SI + App-based tracing, 50% WFH",
                  "— SI + manual tracing (familiar only), 0% WFH",
                  "- - SI + Manual tracing (familiar only), 50% WFH")
  
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
    
    text(x=1,y=2.5-0.1*(ii-1), labels=label_list[ii],cex=0.7,col=col_def[[kk]],adj=0)
  
  }
  
  title(LETTERS[2],adj=0)
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Plot by WFH scenarios
  
  label_list <- c("SI + manual tracing (familiar only)",
                  "SI + manual tracing (all)",
                  "SI + app-based tracing")
  
  store_dat0 <- store_dat %>% filter(limit_other==4 & app_cov==0.53 & trace_p==0.95) %>% arrange(wfh)
  
  for(ii in 1:3){
    
    if(ii==1){
      xx <- store_dat0 %>% filter(scenario=="isolation_manual_tracing_met_only")
      plot(xx$wfh,xx$r_eff,xlab="proportion with no work contacts",ylab="R",ylim=c(0.5,2.5),xlim=c(0,0.8),col="white",type="l",lwd=2)
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
    
  }
  
  
  title(LETTERS[3],adj=0)
  
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Plot by app coverage
  
  label_list <- c("— SI + app-based tracing",
                  "- - SI + app-based (max 4 other contacts)")
  
  store_dat0 <- store_dat %>% filter(limit_other==4 & wfh==0 & trace_p==0.95) %>% arrange(app_cov)
  
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
    
    text(x=0,y=2.5-0.1*(ii-1), labels=label_list[ii],cex=0.7,col=col_def[[kk]],adj=0)
    
  }
  
  
  title(LETTERS[4],adj=0)
  
  dev.copy(png,paste0(dir_pick,"plot2.png"),units="cm",width=20,height=16,res=150)
  dev.off()
  
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
  label_list <- c("SI + manual tracing (familiar only)",
                  "SI + manual (familiar, max 4 other contacts)",
                  "SI + manual tracing (all)",
                  "SI + app-based tracing",
                  "SI + app-based (max 4 other contacts)")
  
  store_dat0 <- store_dat %>% filter(trace_p==0.95 & app_cov==0.53 & prob_t_asymp==0.5 & limit_other==4) %>% arrange(prob_symp)
  
  for(ii in 1:5){
    
    if(ii==1){
      xx <- store_dat0 %>% filter(wfh==0  & scenario=="isolation_manual_tracing_met_only") 
      plot(xx$prob_symp,xx$reduction_raw,xlab="proportion symptomatic",ylab="relative R",ylim=c(0,1),xlim=c(0.2,0.8),col="white",type="l",lwd=2)
      #grid(ny = NULL, nx=NA, col = "lightgray")
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
    
    lines(xx$prob_symp,xx$reduction_raw,col=col_def[[kk]],lwd=2) 
    
    text(x=0.2,y=1-0.05*(ii-1), labels=label_list[ii],cex=0.7,col=col_def[[kk]],adj=0)
    
  }
  
  title(LETTERS[1],adj=0)
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Plot different asymptomatic transmission assumptions
  label_list <- c("SI + manual tracing (familiar only)",
                  "SI + manual (familiar, max 4 other contacts)",
                  "SI + manual tracing (all)",
                  "SI + app-based tracing",
                  "SI + app-based (max 4 other contacts)")
  
  store_dat0 <- store_dat %>% filter(trace_p==0.95 & app_cov==0.53 & prob_symp==0.6 & limit_other==4) %>% arrange(prob_t_asymp)
  
  for(ii in 1:5){
    
    if(ii==1){
      xx <- store_dat0 %>% filter(wfh==0  & scenario=="isolation_manual_tracing_met_only") 
      plot(xx$prob_t_asymp,xx$reduction_raw,xlab="relative asymptomatic transmission",ylab="relative R",ylim=c(0,1),
           xlim=c(0,1),col="white",type="l",lwd=2)
      #grid(ny = NULL, nx=NA, col = "lightgray")
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
    
    text(x=0.2,y=1-0.05*(ii-1), labels=label_list[ii],cex=0.7,col=col_def[[kk]],adj=0)
    
  }
  
  title(LETTERS[2],adj=0)
  
  dev.copy(png,paste0(dir_pick,"plot_symp.png"),units="cm",width=20,height=8,res=150)
  dev.off()
  
}

