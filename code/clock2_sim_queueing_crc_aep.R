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
generate_inquisit_lists <- F

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
idf_list <- rob_grid %>% group_split(iteration)
niterations <- length(idf_list)
for (f in 1:niterations) {data.table::fwrite(idf_list[[f]], file = paste0("grid_", f, ".csv"))}

if (sum(stringr::str_detect(Sys.info(), "andypapale"))>1)  {
  source('~/clock2/code/clock2_sim_crc.R')
  setwd(output_dir)
  i = 1
  j = 1
  files <- list.files(pattern = "grid_")
  df <- data.table::fread(file = files[j])
  set.seed(df$iteration[i])
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
  tt <- iterate_sim(df, bump_prominence, ncenters, centers, values, width_sd, i, j)
  bb <- tt
  cc <- round(bb$get_values_matrix(),0) 
  rm(tt)
  
  set.seed(df$iteration[i])
  ncenters <- 9 # how many gaussians there are
  mean_val <- 10 # mean reward rate
  sd_val <- 2 # standard deviation of reward / range of rewards
  centers <- sample(seq(0, 2*pi, by = pi/20), ncenters, replace = FALSE) # line up gaussians here
  values <- sample(truncnorm::rtruncnorm(ncenters, a = 0, mean = mean_val, sd = sd_val))
  width_sd <- 0.349 # fixed, how wide are the underlying Gaussians
  sanity_checks = T # diagnostic plots inside simulation loop
  ntrials = 300
  i = 1
  j = 1
  cat(sprintf("In loop i: %d, j: %d\n", i, j), file = "run_log.txt", append=T)
  # set up contingency
  bump_prominence <- 10
  bump_value <- mean_val * bump_prominence
  bump_center <- sample(seq(0, 2*pi, by = pi/20), 1, replace = FALSE)
  setwd(base_dir)
  contingency <- vm_circle_contingency(centers = c(centers, bump_center), weights = c(values, bump_value), widths = rep(width_sd, ncenters + 1), units = "radians")
  qq <- troll_world$new(n_trials=ntrials, values=contingency$get_wfunc(), drift_sd=1)
  qq$apply_flex(high_avg = 1, high_spread = 0, low_avg = df$low_avg[i], spread_max = 100, jump_high = T)
  values <- data.frame(round(t(qq$get_values_matrix())),0)
  #values <- values %>% mutate(timepoint = row_number()) %>% rowwise() %>% pivot_longer(cols = starts_with("X"), names_to = "trial")
  values <- values %>% mutate(timepoint = row_number()) %>% rowwise() %>% pivot_longer(cols = starts_with("X"), names_to = "trial") %>%
    mutate(trial = extract_numeric(trial)) %>% group_by(trial) %>% summarise(vmax = max(value),
                                                                             vmax_location = timepoint[which.max(value)])
  plot(values$vmax_location)
  
  inq_tri <- round(qq$get_values_matrix(),0) # original value matrix
  aa <- inq_tri
  
  if (generate_inquisit_lists){
    inq_tri <- data.frame(inq_tri)
    inq_tri <- inq_tri %>% mutate(trial = row_number()) %>% rowwise() %>% pivot_longer(cols = starts_with("X"), names_to = "RT") %>% mutate(RT = extract_numeric(RT))
    inq_tri <- inq_tri %>% arrange(trial,RT)
    write.csv(inq_tri,'Design-Matrix-6820.csv')
    epoch <- data.frame(qq$epoch)
    write.csv(epoch,'epoch-6820.csv')
    
    setwd('~/clock2')
    # generate value, RT and trial lists as 1 x (nT x nRT) inquisit lists
    options("encoding" = "UTF-8") # encode in UTF-8 as suggested here https://forums.millisecond.com/Topic15777.aspx#15778
    df0 <- NULL;
    dq0 <- NULL;
    dz0 <- NULL;
    nR <- nrow(inq_tri);
    for (iR in 1:nR){
      if (iR==1){
        df0 <- paste0('<list values>\n/ items = (',as.character(inq_tri$value[iR]),',');
        dq0 <- paste0('<list RT>\n/ items = (',as.character(inq_tri$RT[iR]),',');
        dz0 <- paste0('<list trial>\n/ items = (',as.character(inq_tri$trial[iR]),',');
      } else if (iR > 1 && iR < nR){
        df0 <- paste0(df0,as.character(inq_tri$value[iR]),',');
        dq0 <- paste0(dq0,as.character(inq_tri$RT[iR]),',');
        dz0 <- paste0(dz0,as.character(inq_tri$trial[iR]),',');
      } else if (iR==nR){
        df0 <- paste0(df0,as.character(inq_tri$value[iR]),')\n/ selectionrate = always\n/ selectionmode = values.master_idx;\n</list>')
        dq0 <- paste0(dq0,as.character(inq_tri$RT[iR]),')\n/ selectionrate = always\n/ selectionmode = values.master_idx;\n</list>')
        dz0 <- paste0(dz0,as.character(inq_tri$trial[iR]),')\n/ selectionrate = always\n/ selectionmode = values.master_idx;\n</list>')
      }
      if ((iR %% 1000)==0){
        print(iR/nR);
      }
    }
    write.table(df0,'values-6820.txt',row.names=F,col.names=F,quote=F)
    write.table(dq0,'RTs-6820.txt',row.names=F,col.names=F,quote=F)
    write.table(dz0,'trials-6820.txt',row.names=F,col.names=F,quote=F)
    options("encoding" = "native.enc") # change encoding back to native
    
    
    # write erasure schedule
    era_loc <- zero_to_2pi((bb$erasure_segments$segment_max+bb$erasure_segments$segment_min)/2)*180/pi
    trial_type <- bb$erasure_segments$trial_type
    drift <- bb$get_drift_vec()
    spread <- bb$spread
    options("encoding" = "UTF-8")
    df0 <- NULL;
    dq0 <- NULL;
    dz0 <- NULL;
    dh0 <- NULL;
    nR <- length(era_loc);
    for (iR in 1:nR){
      
      if (is.na(era_loc[iR])){
        temp <- 'NULL';
      } else {
        temp <- era_loc[iR];
      }
      
      if (iR==1){
        df0 <- paste0('<list erasure_locations>\n/ items = (',as.character(temp),',');
        dq0 <- paste0('<list trial_type>\n/ items = ("',as.character(trial_type[iR]),'",');
        dz0 <- paste0('<list drift>\n/ items = (',as.character(drift[iR]),',');
        dh0 <- paste0('<list spread>\n/ items = (',as.character(spread[iR]),',');
      } else if (iR > 1 && iR < nR){
        df0 <- paste0(df0,as.character(temp),',');
        dq0 <- paste0(dq0,'"',as.character(trial_type[iR]),'",');
        dz0 <- paste0(dz0,as.character(drift[iR]),',');
        dh0 <- paste0(dh0,as.character(spread[iR]),',');
      } else if (iR==nR){
        df0 <- paste0(df0,as.character(temp),')\n/ selectionrate = always\n/ selectionmode = values.era_loc_index;\n</list>')
        dq0 <- paste0(dq0,'"',as.character(trial_type[iR]),'")\n/ selectionrate = always\n/ selectionmode = values.trial;\n</list>')
        dz0 <- paste0(dz0,as.character(drift[iR]),')\n/ selectionrate = always\n/ selectionmode = values.trial; \n</list>')
        dh0 <- paste0(dh0,as.character(spread[iR]),')\n/ selectionrate = always\n/ selectionmode = values.trial; \n</list>')
      }
    }
    write.table(df0,'era_loc-6820.txt',row.names=F,col.names=F,quote=F)
    write.table(dq0,'trial_type-6820.txt',row.names=F,col.names=F,quote=F)
    write.table(dz0,'drift-vector-6820.txt',row.names=F,col.names=F,quote=F)
    write.table(dh0,'spread-6820.txt',row.names=F,col.names=F,quote=F)
    options("encoding" = "native.enc") # change encoding back to native
    
    era_val = NULL
    era_loc_last = 0
    iC = 1;
    for (iT in 1:300){
      if (bb$erasure_segments$trial_type[iT] == "erasure" && era_loc_last != era_loc[iT]){
        era_loc_last = era_loc[iT]
        print(cc[iT,era_loc[iT]])
        era_val[iC] = cc[iT,era_loc[iT]];
        iC = iC + 1;
      }
    }  
  }

  
  if (animate==TRUE){
    loc <- round(bb$get_pvec(), 2)
    setwd('~/clock2/animation/')
    for (ii in 1:nrow(cc)) {
      if (bb$erasure_segments$trial_type[ii] == "erasure") {
        ss <- paste(", er:", round(bb$erasure_segments[ii,"segment_min"], 2))
      } else {
        ss <- ""
      }
      pdf(paste0('image-',ii,'.pdf'),height=12,width=12)
      gg1 <- plot(cc[ii,], type="l", main=paste("Trial", ii, "epoch", bb$epoch[ii], ss), ylim = range(cc), xaxt='n')
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
