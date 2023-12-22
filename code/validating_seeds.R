rob_grid <- expand.grid(alpha = c(0.2), gamma = c(0.1),                 # model params
                        beta = c(1), # at very high betas, h and u are decorrelated, no need to test
                        epsilon_u = c(0.9999), # 0.0833 is at chance, low correlation -- not worth testing
                        block_length = c(10), # block length > 15 had higher correlations, not worth testing
                        low_avg = c(10),
                        iteration = c(6520),
                        #drift = c(1, 2, 4), bump_prom = c(8, 10, 15),
                        seed = 1)

source('~/clock2/code/clock2_sim_crc_aep.R')
setwd(output_dir)
i = 1
j = 1
set.seed(rob_grid$iteration[i])
ncenters <- 9 # how many gaussians there are
mean_val <- 10 # mean reward rate
sd_val <- 2 # standard deviation of reward / range of rewards
centers <- sample(seq(0, 2*pi, by = pi/20), ncenters, replace = FALSE) # line up gaussians here
values <- sample(truncnorm::rtruncnorm(ncenters, a = 0, mean = mean_val, sd = sd_val))
width_sd <- 0.349 # fixed, how wide are the underlying Gaussians
sanity_checks = F # diagnostic plots inside simulation loop
ntrials = 300
cat(sprintf("In loop i: %d, j: %d\n", i, j), file = "run_log.txt", append=T)
bump_prominence <- 10
bump_value <- mean_val * bump_prominence
bump_center <- sample(seq(0, 2*pi, by = pi/20), 1, replace = FALSE)
setwd(base_dir)
tt <- iterate_sim(rob_grid, bump_prominence, ncenters, centers, values, width_sd, i, j)
cc <- round(tt$get_values_matrix(type = 'objective', quiet = F),0)
 # for (r in 1:nrow(cc)) {
 #    plot(cc[r,])
 #    print(r)
 #    Sys.sleep(.2)
 # }
rm(tt)
gg <- iterate_sim(rob_grid, bump_prominence, ncenters, centers, values, width_sd, i, j)
cc <- round(gg$get_values_matrix(type = 'objective', quiet = F),0)
# for (r in 1:nrow(cc)) {
#    plot(cc[r,])
#    print(r)
#    Sys.sleep(.2)
# }

inq_tri <- data.frame(cc)
inq_tri <- inq_tri %>% mutate(trial = row_number()) %>% rowwise() %>% pivot_longer(cols = starts_with("X"), names_to = "RT") %>% mutate(RT = readr::parse_number(RT))
inq_tri <- inq_tri %>% arrange(trial,RT)
setwd('/Users/andypapale/clock2/Inquisit_design_files/')
write.csv(inq_tri,paste0('Design-Matrix-with-Erasures-',as.character(rob_grid$iteration),'.csv'))

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
dd <- iterate_sim(rob_grid, bump_prominence, ncenters, centers, values, width_sd, i, j)
rm(dd)
contingency <- vm_circle_contingency(centers = c(centers, bump_center), weights = c(values, bump_value), widths = rep(width_sd, ncenters + 1), units = "radians")
qq <- troll_world$new(n_trials=ntrials, values=contingency$get_wfunc(), drift_sd=1)
qq$apply_flex(high_avg = 1, high_spread = 0, low_avg = rob_grid$low_avg[i], spread_max = 100, jump_high = T)
dd <- round(qq$get_values_matrix(type = 'objective', quiet = F),0)

# for (r in 1:nrow(dd)) {
#   plot(dd[r,])
#   print(r)
#   Sys.sleep(.2)
# }


df <- read_csv('/Users/andypapale/Inquisit Code/EEG_clock/Clock_v2/data/clock_2_1_1_seed_6520_testing_raw_1_2023-12-22-18-49-18-373.csv')
#df <- read_csv('/Users/andypapale/Inquisit Code/EEG_clock/Clock_v2/data/clock_2_1_1_seed_5815_testing_raw_1_2023-12-21-22-36-42-177.csv')
df0 <- df %>% filter(trialcode=="getValueVector")
df0$trial <- df0$trial...18 # can make this step not necessary in final program
d9 <- read_csv(paste0('/Users/andypapale/clock2/Inquisit_design_files/Design-Matrix-',as.character(rob_grid$iteration),'.csv'))

d9_w_era <- read_csv(paste0('/Users/andypapale/clock2/Inquisit_design_files/Design-Matrix-with-Erasures-',as.character(rob_grid$iteration),'.csv'))

if (!file.exists(paste0('/Users/andypapale/Desktop/seed-',as.character(rob_grid$iteration),'-validation'))){
  dir.create(paste0('/Users/andypapale/Desktop/seed-',as.character(rob_grid$iteration),'-validation'))
}
setwd(paste0('/Users/andypapale/Desktop/seed-',as.character(rob_grid$iteration),'-validation'))

columns = c('e_pos','e_val');

for (iT in 1:300){
  e_val <- d9_w_era %>% filter(trial==iT);
  e_val <- e_val$value[round(e_val$RT)==round(era_loc[iT])];
  df_era <- data.frame(e_pos = era_loc[iT],e_val = e_val)
  inq_e_val = df0 %>% filter(trial==iT); 
  df_era <- data.frame(matrix(nrow = 0, ncol = length(columns)))
  colnames(df_era) = columns
  if (gg$erasure_segments$trial_type[iT]=='erasure'){
    era_loc2 <- inq_e_val$erasure_RT[1];
    e_val_inq <- inq_e_val$erasure_value[1];
    df_era <- data.frame(e_pos = era_loc[iT],e_val = e_val)
  } else if (gg$erasure_segments$trial_type[iT]=='attention'){
    era_loc2 <- inq_e_val$stim_center_deg[1];
    e_val_inq <- NULL;
  } else if (gg$erasure_segments$trial_type[iT]=='no erasure'){
    era_loc2 <- NULL;
    e_val_inq <- NULL;
  }
  pdf(paste0('seed-',as.character(rob_grid$iteration),'-trial-',iT,'.pdf'),height=12,width=12); 
  if (nrow(df_era)==1){
    gg1 <- ggplot(d9 %>% filter(trial==iT), aes(x=RT,y=value), color='black') + 
     geom_line() + geom_point(size=5,color='black') + 
     geom_line(data = df0 %>% filter(trial==iT),aes(x=x_final,y=y_final),color='red') + 
     geom_point(data = df0 %>% filter(trial==iT),aes(x=x_final,y=y_final),color='red') + 
     ylim(0,200) + geom_point(data = df_era, aes(x=e_pos,y=e_val),color='black',shape=2,size=10) +
     ggtitle(paste0('trial ',iT,' trial_type = ', gg$erasure_segments$trial_type[iT], ' trial_type_inq = ',inq_e_val$trial_type[1],' e_val = ',e_val, ' inq_e_val =', e_val_inq, ' e_max = ',round(gg$erasure_segments$segment_max[iT]*180/pi), ' inq_e_loc = ', era_loc2)); 
  } else {
    gg1 <-  ggplot(d9 %>% filter(trial==iT), aes(x=RT,y=value), color='black') + 
      geom_line() + geom_point(size=5,color='black') + 
      geom_line(data = df0 %>% filter(trial==iT),aes(x=x_final,y=y_final),color='red') + 
      geom_point(data = df0 %>% filter(trial==iT),aes(x=x_final,y=y_final),color='red') + 
      ylim(0,200) +
      ggtitle(paste0('trial ',iT,' trial_type = ', gg$erasure_segments$trial_type[iT], ' trial_type_inq = ',inq_e_val$trial_type[1])); 
  }
      # gg1 <- ggplot(d9 %>% filter(trial==iT), aes(x=RT,y=value)) + geom_line() + geom_point(size=5,color='black') +  ylim(0,200) +
  #   ggtitle(paste0('trial ',iT,' trial_type = ', bb$erasure_segments$trial_type[iT], ' e_val = ',e_val, ' e_max = ',round(tt$erasure_segments$segment_max[iT]*180/pi)))
  print(gg1); 
  dev.off()
}
