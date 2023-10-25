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
animate <- FALSE
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
                          iteration = c(5888),
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
  
  values <- data.frame(round(t(tt$get_values_matrix())),0)
  #values <- values %>% mutate(timepoint = row_number()) %>% rowwise() %>% pivot_longer(cols = starts_with("X"), names_to = "trial")
  values <- values %>% mutate(timepoint = row_number()) %>% rowwise() %>% pivot_longer(cols = starts_with("X"), names_to = "trial") %>%
    mutate(trial = extract_numeric(trial)) %>% group_by(trial) %>% summarise(vmax = max(value),
                                                                             vmax_location = timepoint[which.max(value)])
  plot(values$vmax_location)
  
  inq_tri <- round(tt$get_values_matrix("objective", quiet=F),0) # all manipulations, matrix of expected values
  aa <- inq_tri
  inq_tri <- data.frame(inq_val)
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
  write.table(df0,'values-5888.txt',row.names=F,col.names=F,quote=F)
  write.table(dq0,'RTs-5888.txt',row.names=F,col.names=F,quote=F)
  write.table(dz0,'trials-5888.txt',row.names=F,col.names=F,quote=F)
  options("encoding" = "native.enc") # change encoding back to native
  
  
  # write erasure schedule
  era_loc <- zero_to_2pi((tt$erasure_segments$segment_max+tt$erasure_segments$segment_min)/2)*180/pi
  trial_type <- tt$erasure_segments$trial_type
  options("encoding" = "UTF-8")
  df0 <- NULL;
  dq0 <- NULL;
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
    } else if (iR > 1 && iR < nR){
      df0 <- paste0(df0,as.character(temp),',');
      dq0 <- paste0(dq0,'"',as.character(trial_type[iR]),'",');
    } else if (iR==nR-9){
      df0 <- paste0(df0,as.character(temp),')\n/ selectionrate = always\n/ selectionmode = values.trial;\n</list>')
      dq0 <- paste0(dq0,'"',as.character(trial_type[iR]),'")\n/ selectionrate = always\n/ selectionmode = values.trial;\n</list>')
    }
  }
  write.table(df0,'era_loc-5888.txt',row.names=F,col.names=F,quote=F)
  write.table(dq0,'trial_type-5888.txt',row.names=F,col.names=F,quote=F)
  options("encoding" = "native.enc") # change encoding back to native
  
  
  
  
  
  if (animate==TRUE){
    loc <- round(tt$get_pvec(), 2)
    setwd('~/clock2/animation/')
    for (ii in 1:nrow(aa)) {
      if (tt$erasure_segments$trial_type[ii] == "erasure") {
        ss <- paste(", er:", round(tt$erasure_segments[ii,"segment_min"], 2))
      } else {
        ss <- ""
      }
      pdf(paste0('image-',ii,'.pdf'),height=12,width=12)
      gg1 <- plot(aa[ii,], type="l", main=paste("Trial", ii, "epoch", tt$epoch[ii], ss), ylim = range(aa), xaxt='n')
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