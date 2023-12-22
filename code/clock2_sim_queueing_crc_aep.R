# R script for handling cluster queuing for clock2 simulations
library(tidyverse)
library(BAMBI)
if (sum(stringr::str_detect(Sys.info(), "Alex|alexdombrovski"))>1) {
  basedir <- "~/code/clock2/code/"
  output_dir <- "~/code/clock2/simulations"
  sbatch_dir <- "~/code/clock2/code/sbatch/"
} else if (sum(stringr::str_detect(Sys.info(), "andypapale"))>1)  {
  basedir <- "~/clock2/code/"
  output_dir <- "~/clock2/simulations"
  sbatch_dir <- "~/clock2/code/sbatch/"
} else {
  basedir <- "~/code/clock2/code/"
  output_dir <- "~/code/clock2/simulations"
  sbatch_dir <- "~/code/clock2/code/sbatch/"
}

test <- F
test_on_mac <- F
animate <- F
# if (test && test_on_mac) {
#   basedir <- "~/OneDrive - University of Pittsburgh/Momentum_EMA/eeg_data_t_split/"
#   output_dir <- "~/OneDrive - University of Pittsburgh/Momentum_EMA/"
# }
setwd(output_dir)
silent <- F
generate_inquisit_lists <- T

# set up simulation grid, write files
if (sum(stringr::str_detect(Sys.info(), "andypapale"))>1)  {
  rob_grid <- expand.grid(alpha = c(0.2), gamma = c(0.1),                 # model params
                          beta = c(1), # at very high betas, h and u are decorrelated, no need to test
                          epsilon_u = c(0.9999), # 0.0833 is at chance, low correlation -- not worth testing
                          block_length = c(10), # block length > 15 had higher correlations, not worth testing
                          low_avg = c(10),
                          iteration = c(6520),
                          #drift = c(1, 2, 4), bump_prom = c(8, 10, 15),
                          seed = 1)
} else {
  rob_grid <- expand.grid(alpha = c(0.2, 0.5), gamma = c(0.1, 0.5, 0.9),                 # model params
                          beta = c(1, 5), # at very high betas, h and u are decorrelated, no need to test
                          epsilon_u = c(0.3, 0.9), # 0.0833 is at chance, low correlation -- not worth testing
                          block_length = c(10), # block length > 15 had higher correlations, not worth testing
                          low_avg = c(10, 20),
                          iteration = c(223:1222),
                          #drift = c(1, 2, 4), bump_prom = c(8, 10, 15),
                          seed = 1:100)
}
if (!sum(stringr::str_detect(Sys.info(), "andypapale"))>1)  {
  idf_list <- rob_grid %>% group_split(iteration)
  niterations <- length(idf_list)
  for (f in 1:niterations) {data.table::fwrite(idf_list[[f]], file = paste0("grid_", f, ".csv"))}
}
if (sum(stringr::str_detect(Sys.info(), "andypapale"))>1)  {
  source('~/clock2/code/clock2_sim_crc_aep.R')
  setwd(output_dir)
  i = 1
  j = 1
  #files <- list.files(pattern = "grid_")
  #df <- data.table::fread(file = files[j])
  set.seed(rob_grid$iteration[i])
  ncenters <- 9 # how many gaussians there are
  mean_val <- 10 # mean reward rate
  sd_val <- 2 # standard deviation of reward / range of rewards
  centers <- sample(seq(0, 2*pi, by = pi/20), ncenters, replace = FALSE) # line up gaussians here
  values <- sample(truncnorm::rtruncnorm(ncenters, a = 0, mean = mean_val, sd = sd_val))
  width_sd <- 0.349 # fixed, how wide are the underlying Gaussians
  sanity_checks = T # diagnostic plots inside simulation loop
  ntrials = 300
  cat(sprintf("In loop i: %d, j: %d\n", i, j), file = "run_log.txt", append=T)
  bump_prominence <- 10
  bump_value <- mean_val * bump_prominence
  bump_center <- sample(seq(0, 2*pi, by = pi/20), 1, replace = FALSE)
  setwd(base_dir)
  dum_run <- iterate_sim(rob_grid, bump_prominence, ncenters, centers, values, width_sd, i, j)
  rm(dum_run)
  tt <- iterate_sim(rob_grid, bump_prominence, ncenters, centers, values, width_sd, i, j)
  cc <- round(tt$get_values_matrix(type = 'objective'),0) 
  inq_tri <- data.frame(cc)
  inq_tri <- inq_tri %>% mutate(trial = row_number()) %>% rowwise() %>% pivot_longer(cols = starts_with("X"), names_to = "RT") %>% mutate(RT = readr::parse_number(RT))
  inq_tri <- inq_tri %>% arrange(trial,RT)
  write.csv(inq_tri,paste0('Design-Matrix-with-Erasures-',as.character(rob_grid$iteration),'.csv'))
  values1 <- data.frame(t(cc))
  values1 <- values1 %>% mutate(timepoint = row_number()) %>% rowwise() %>% pivot_longer(cols = starts_with("X"), names_to = "trial") %>%
    mutate(trial = extract_numeric(trial)) %>% group_by(trial) %>% summarise(vmax = max(value),
                                                                             vmax_location = timepoint[which.max(value)])
  plot(values1$vmax_location)
  
  set.seed(rob_grid$iteration[i])
  ncenters <- 9 # how many gaussians there are
  mean_val <- 10 # mean reward rate
  sd_val <- 2 # standard deviation of reward / range of rewards
  centers <- sample(seq(0, 2*pi, by = pi/20), ncenters, replace = FALSE) # line up gaussians here
  values <- sample(truncnorm::rtruncnorm(ncenters, a = 0, mean = mean_val, sd = sd_val))
  width_sd <- 0.349 # fixed, how wide are the underlying Gaussians
  sanity_checks = F # diagnostic plots inside simulation loop
  ntrials = 300
  i = 1
  j = 1
  cat(sprintf("In loop i: %d, j: %d\n", i, j), file = "run_log.txt", append=T)
  # set up contingency
  bump_prominence <- 10
  bump_value <- mean_val * bump_prominence
  bump_center <- sample(seq(0, 2*pi, by = pi/20), 1, replace = FALSE)
  setwd(base_dir)
  dum_run <- iterate_sim(rob_grid, bump_prominence, ncenters, centers, values, width_sd, i, j)
  rm(dum_run)
  contingency <- vm_circle_contingency(centers = c(centers, bump_center), weights = c(values, bump_value), widths = rep(width_sd, ncenters + 1), units = "radians")
  qq <- troll_world$new(n_trials=ntrials, values=contingency$get_wfunc(), drift_sd=1)
  qq$apply_flex(high_avg = 1, high_spread = 0, low_avg = rob_grid$low_avg[i], spread_max = 100, jump_high = T)
  values2 <- data.frame(round(t(qq$get_values_matrix())),0)
  #values <- values %>% mutate(timepoint = row_number()) %>% rowwise() %>% pivot_longer(cols = starts_with("X"), names_to = "trial")
  values2 <- values2 %>% mutate(timepoint = row_number()) %>% rowwise() %>% pivot_longer(cols = starts_with("X"), names_to = "trial") %>%
    mutate(trial = extract_numeric(trial)) %>% group_by(trial) %>% summarise(vmax = max(value),
                                                                             vmax_location = timepoint[which.max(value)])
  plot(values2$vmax_location)
  aa <- round(qq$get_values_matrix(type = 'objective',quiet = F),0) 
  
  # for (r in 1:nrow(aa)) {
  #   plot(aa[r,])
  #   print(r)
  #   Sys.sleep(.2)
  # }
  
  if (generate_inquisit_lists){
    
    setwd('/Users/andypapale/clock2/Inquisit_design_files')
    
    aa <- data.frame(aa)
    aa <- aa %>% mutate(trial = row_number()) %>% rowwise() %>% pivot_longer(cols = starts_with("X"), names_to = "RT") %>% mutate(RT = extract_numeric(RT))
    aa <- aa %>% arrange(trial,RT)
    write.csv(aa,paste0('Design-Matrix-',as.character(rob_grid$iteration),'.csv'))
    epoch <- data.frame(tt$epoch)
    write.csv(epoch,paste0('epoch-',as.character(rob_grid$iteration),'.csv'))
    
    #generate value, RT and trial lists as 1 x (nT x nRT) inquisit lists
    options("encoding" = "UTF-8") # encode in UTF-8 as suggested here https://forums.millisecond.com/Topic15777.aspx#15778
    df0 <- NULL;
    dq0 <- NULL;
    dz0 <- NULL;
    nR <- nrow(aa);
    for (iR in 1:nR){
      if (iR==1){
        df0 <- paste0('<list values>\n/ items = (',as.character(aa$value[iR]),',');
        dq0 <- paste0('<list RT>\n/ items = (',as.character(aa$RT[iR]),',');
        dz0 <- paste0('<list trial>\n/ items = (',as.character(aa$trial[iR]),',');
      } else if (iR > 1 && iR < nR){
        df0 <- paste0(df0,as.character(aa$value[iR]),',');
        dq0 <- paste0(dq0,as.character(aa$RT[iR]),',');
        dz0 <- paste0(dz0,as.character(aa$trial[iR]),',');
      } else if (iR==nR){
        df0 <- paste0(df0,as.character(aa$value[iR]),')\n/ selectionrate = always\n/ selectionmode = values.master_idx;\n</list>')
        dq0 <- paste0(dq0,as.character(aa$RT[iR]),')\n/ selectionrate = always\n/ selectionmode = values.master_idx;\n</list>')
        dz0 <- paste0(dz0,as.character(aa$trial[iR]),')\n/ selectionrate = always\n/ selectionmode = values.master_idx;\n</list>')
      }
      if ((iR %% 1000)==0){
        print(iR/nR);
      }
    }
    write.table(df0,paste0('values-',as.character(rob_grid$iteration),'.txt'),row.names=F,col.names=F,quote=F)
    write.table(dq0,paste0('RTs-',as.character(rob_grid$iteration),'.txt'),row.names=F,col.names=F,quote=F)
    write.table(dz0,paste0('trials-',as.character(rob_grid$iteration),'.txt'),row.names=F,col.names=F,quote=F)
    options("encoding" = "native.enc") # change encoding back to native
    
    
    # write erasure schedule
    seg_min <- tt$erasure_segments$segment_min*180/pi;
    seg_max <- tt$erasure_segments$segment_max*180/pi;
    era_loc <- NULL;
    for (iT in 1:300){
      if (!is.na(seg_max[iT])){
      if (seg_max[iT] < 31){
        # then seg_min >= 360
        if (seg_max[iT] < 16){
          # then era_loc > 15
          era_loc[iT] <- seg_min[iT] + 15;
        } else if (seg_max[iT] >= 16){
          # then era_loc <= 15
          era_loc[iT] <- seg_max[iT] - 15;
        }
      } else {
        era_loc[iT] <- seg_min[iT] + 15;
      }
      } else {
        era_loc[iT] <- NA;
      }
    }
    trial_type <- tt$erasure_segments$trial_type
    drift <- tt$get_drift_vec()
    spread <- tt$spread
    options("encoding" = "UTF-8")
    df0 <- NULL;
    dq0 <- NULL;
    dz0 <- NULL;
    dh0 <- NULL;
    nR <- length(era_loc);
    for (iR in 1:nR){
      
      #if (is.na(era_loc[iR])){
      #  temp <- 'NULL';
      #} else {
      #  temp <- era_loc[iR];
      #}
      
      if (iR==1){
        #df0 <- paste0('<list erasure_locations>\n/ items = (',as.character(temp),',');
        dq0 <- paste0('<list block_type>\n/ items = ("',as.character(trial_type[iR]),'",');
        dz0 <- paste0('<list drift>\n/ items = (',as.character(drift[iR]),',');
        dh0 <- paste0('<list spread>\n/ items = (',as.character(spread[iR]),',');
      } else if (iR > 1 && iR < nR){
        #df0 <- paste0(df0,as.character(temp),',');
        dq0 <- paste0(dq0,'"',as.character(trial_type[iR]),'",');
        dz0 <- paste0(dz0,as.character(drift[iR]),',');
        dh0 <- paste0(dh0,as.character(spread[iR]),',');
      } else if (iR==nR){
        #df0 <- paste0(df0,as.character(temp),')\n/ selectionrate = always\n/ selectionmode = values.era_loc_index;\n</list>')
        dq0 <- paste0(dq0,'"',as.character(trial_type[iR]),'")\n/ selectionrate = always\n/ selectionmode = values.trial;\n</list>')
        dz0 <- paste0(dz0,as.character(drift[iR]),')\n/ selectionrate = always\n/ selectionmode = values.trial; \n</list>')
        dh0 <- paste0(dh0,as.character(spread[iR]),')\n/ selectionrate = always\n/ selectionmode = values.trial; \n</list>')
      }
    }
    #write.table(df0,paste0('era_loc-',as.character(rob_grid$iteration),'.txt'),row.names=F,col.names=F,quote=F)
    write.table(dq0,paste0('block_type-',as.character(rob_grid$iteration),'.txt'),row.names=F,col.names=F,quote=F)
    write.table(dz0,paste0('drift-vector-',as.character(rob_grid$iteration),'.txt'),row.names=F,col.names=F,quote=F)
    write.table(dh0,paste0('spread-',as.character(rob_grid$iteration),'.txt'),row.names=F,col.names=F,quote=F)
    options("encoding" = "native.enc") # change encoding back to native
    
    era_val = NULL
    era_loc1 = NULL
    att_loc = NULL
    iC = 1;
    iD = 1;
    for (iT in 1:300){
      if (tt$erasure_segments$trial_type[iT]=='erasure' && tt$erasure_segments$clicks_remain[iT]==2){
        era_val[iC] <- cc[iT,era_loc[iT]];
        era_loc1[iC] <- era_loc[iT];
        iC <- iC + 1;
      }  
      if (tt$erasure_segments$trial_type[iT]=='attention' && (tt$erasure_segments$trial[iT] %% 10) == 0){
        att_loc[iD] <- era_loc[iT];
        iD <- iD + 1;
      }
    }
    
    options("encoding" = "UTF-8")
    dq0 <- NULL;
    dh0 <- NULL;
    nR <- length(era_loc1);
    for (iR in 1:nR){
      if (iR==1){
        dq0 <- paste0('<list era_loc>\n/ items = (',as.integer(era_loc1[iR]),',');
        dh0 <- paste0('<list era_val>\n/ items = (',as.character(era_val[iR]),',');
      } else if (iR > 1 && iR < nR){
        dq0 <- paste0(dq0,as.integer(era_loc1[iR]),',');
        dh0 <- paste0(dh0,as.character(era_val[iR]),',');
      } else if (iR==nR){
        dq0 <- paste0(dq0,as.integer(era_loc1[iR]),')\n/ selectionrate = always\n/ selectionmode = values.era_index;\n</list>')
        dh0 <- paste0(dh0,as.character(era_val[iR]),')\n/ selectionrate = always\n/ selectionmode = values.era_index; \n</list>')
      }
    }
    write.table(dq0,paste0('era_loc-',as.character(rob_grid$iteration),'.txt'),row.names=F,col.names=F,quote=F)
    write.table(dh0,paste0('era_val-',as.character(rob_grid$iteration),'.txt'),row.names=F,col.names=F,quote=F)
    options("encoding" = "native.enc") # change encoding back to native
    
    options("encoding" = "UTF-8")
    dq0 <- NULL;
    nR <- length(att_loc);
    for (iR in 1:nR){
      if (iR==1){
        dq0 <- paste0('<list att_loc>\n/ items = (',as.integer(att_loc[iR]),',');
      } else if (iR > 1 && iR < nR){
        dq0 <- paste0(dq0,as.integer(att_loc[iR]),',');
      } else if (iR==nR){
        dq0 <- paste0(dq0,as.integer(att_loc[iR]),')\n/ selectionrate = always\n/ selectionmode = values.att_index;\n</list>')
      }
    }
    write.table(dq0,paste0('att_loc-',as.character(rob_grid$iteration),'.txt'),row.names=F,col.names=F,quote=F)
    options("encoding" = "native.enc") # change encoding back to native
    
    
  }
  
  
  
  

  
  if (animate==TRUE){
    loc <- round(tt$get_pvec(), 2)
    setwd('~/clock2/animation/')
    for (ii in 1:nrow(cc)) {
      if (bb$erasure_segments$trial_type[ii] == "erasure") {
        ss <- paste(", er:", round(bb$erasure_segments[ii,"segment_min"], 2))
      } else {
        ss <- ""
      }
      pdf(paste0('image-',ii,'.pdf'),height=12,width=12)
      gg1 <- plot(cc[ii,], type="l", main=paste("Trial", ii, "epoch", tt$epoch[ii], ss), ylim = range(cc), xaxt='n')
      axis(side = 1, at = seq(1, 360, by=30), labels= loc[seq(1, 360, by=30)])
      print(gg1)
      dev.off()
    }
  }
  
  
  
} else {
  
  #for (f in 1:niterations) {
  for (f in 1:2) {
    
    if (!test) {
      system(
        paste0(
          "cd ", sbatch_dir, "; ",
          "sbatch --time=23:00:00 --mem=1g",
          " --export=sourcefilestart=", f,
          " sbatch_clock2_sim.bash"
        )
      )
      #write compute level to temporary file
    }
    if (!silent) {
      cat(
        paste0(
          "cd ", sbatch_dir, "; ",
          "sbatch --time=23:00:00 --mem=1g",
          " --export=sourcefilestart=", f,
          " sbatch_clock2_sim.bash\n"
        )
      )
    }
  }
}
