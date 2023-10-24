# R script for handling cluster queuing for clock2 simulations
library(tidyverse)

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
# if (test && test_on_mac) {
#   basedir <- "~/OneDrive - University of Pittsburgh/Momentum_EMA/eeg_data_t_split/"
#   output_dir <- "~/OneDrive - University of Pittsburgh/Momentum_EMA/"
# }
setwd(output_dir)
silent <- F

# set up simulation grid, write files
if (sum(stringr::str_detect(Sys.info(), "andypapale"))>1)  {
  rob_grid <- expand.grid(alpha = c(0.2), gamma = c(0.1),                 # model params
                          beta = c(1), # at very high betas, h and u are decorrelated, no need to test
                          epsilon_u = c(0.3), # 0.0833 is at chance, low correlation -- not worth testing
                          block_length = c(10), # block length > 15 had higher correlations, not worth testing
                          low_avg = c(10),
                          iteration = c(8014),
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
  files <- list.files(pattern = "grid_")
  df <- data.table::fread(file = files[j])
  set.seed(df$iteration[i])
  ncenters <- 9 # how many gaussians there are
  mean_val <- 10 # mean reward rate
  sd_val <- 2 # standard deviation of reward / range of rewards
  centers <- sample(seq(0, 2*pi, by = pi/20), ncenters, replace = FALSE) # line up gaussians here
  values <- sample(truncnorm::rtruncnorm(ncenters, a = 0, mean = mean_val, sd = sd_val))
  width_sd <- 20 # fixed, how wide are the underlying Gaussians
  sanity_checks = T # diagnostic plots inside simulation loop
  ntrials = 300
  i = 1
  j = 1
  tt <- iterate_sim(df, bump_prominence, ncenters, centers, values, width_sd, i, j)
  
  inq_val <- round(t(tt$get_values_matrix()),0)
  inq_tri <- data.frame(t(inq_val))
  inq_tri <- inq_tri %>% mutate(trial = row_number()) %>% rowwise() %>% pivot_longer(cols = starts_with("X"), names_to = "RT") %>% mutate(RT = extract_numeric(RT))
  inq_tri <- inq_tri %>% arrange(trial,RT)
  
  
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
      df0 <- paste0(df0,as.character(inq_tri$value[iR]),')\n/ selectionrate = always\n</list>')
      dq0 <- paste0(dq0,as.character(inq_tri$RT[iR]),')\n/ selectionrate = always\n</list>')
      dz0 <- paste0(dz0,as.character(inq_tri$trial[iR]),')\n/ selectionrate = always\n</list>')
    }
    if ((iR %% 1000)==0){
      print(iR/nR);
    }
  }
  write.table(df0,'values.txt',row.names=F,col.names=F,quote=F)
  write.table(dq0,'RTs.txt',row.names=F,col.names=F,quote=F)
  write.table(dz0,'trials.txt',row.names=F,col.names=F,quote=F)
  options("encoding" = "native.enc") # change encoding back to native
  
  aa <- tt$get_values_matrix("objective", quiet=F) # all manipulations, matrix of expected values, p_reward =0.7 fixed for now
  loc <- round(tt$get_pvec(), 2)
  for (ii in 1:nrow(aa)) {
    if (tt$erasure_segments$trial_type[ii] == "erasure") {
      ss <- paste(", er:", round(tt$erasure_segments[ii,"segment_min"], 2))
    } else {
      ss <- ""
    }
    plot(aa[ii,], type="l", main=paste("Trial", ii, "epoch", tt$epoch[ii], ss), ylim = range(aa), xaxt='n')
    axis(side = 1, at = seq(1, 360, by=30), labels= loc[seq(1, 360, by=30)])
    Sys.sleep(.2)
    
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