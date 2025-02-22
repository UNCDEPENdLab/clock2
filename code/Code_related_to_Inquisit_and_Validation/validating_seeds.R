df <- data.frame(alpha = c(0.2), gamma = c(0.1),                 # model params
                 beta = c(1), # at very high betas, h and u are decorrelated, no need to test
                 epsilon_u = c(0.9999), # 0.0833 is at chance, low correlation -- not worth testing
                 block_length = c(10), # block length > 15 had higher correlations, not worth testing
                 low_avg = c(10),
                 iteration = c(152),
                 #drift = c(1, 2, 4), bump_prom = c(8, 10, 15),
                 seed = 1)

source('~/clock2/code/clock2_sim_crc_aep.R')
setwd(output_dir)
i = 1
j = 1

## get values with erasures
set.seed(df$iteration[i])
ncenters <- 9 # how many gaussians there are
mean_val <- 10 # mean reward rate
sd_val <- 2 # standard deviation of reward / range of rewards
centers <- sample(seq(0, 2*pi, by = pi/20), ncenters, replace = FALSE) # line up gaussians here
values <- sample(truncnorm::rtruncnorm(ncenters, a = 0, mean = mean_val, sd = sd_val))
width_sd <- 0.349 # fixed, how wide are the underlying Gaussians
ntrials = 300
bump_prominence <- 10
bump_value <- mean_val * bump_prominence
bump_center <- sample(seq(0, 2*pi, by = pi/20), 1, replace = FALSE)
contingency <- vm_circle_contingency(centers = c(centers, bump_center), weights = c(values, bump_value), widths = rep(width_sd, ncenters + 1), units = "radians")
tt <- troll_world$new(n_trials=ntrials, values=contingency$get_wfunc(), drift_sd=1)
tt$apply_flex(high_avg = 1, high_spread = 0, low_avg = df$low_avg[i], spread_max = 100, jump_high = T)
tt$setup_erasure_blocks(disappear_clicks = 2, timeout_trials = 2, block_length = df$block_length[i])
sceptic_agent <- scepticc$new(n_basis=12, n_points=200, contingency=tt)
sceptic_agent$alpha <- alpha <- df$alpha[i]
sceptic_agent$beta <- beta <- df$beta[i]
sceptic_agent$gamma <- gamma <- df$gamma[i]
sceptic_agent$epsilon_u <- epsilon_u <- df$epsilon_u[i]
for (j in 1:length(df$seed)){
  set.seed(df$seed[j])
  learning_history <- sceptic_agent$run_contingency(optimize = FALSE)
}
cc <- round(tt$get_values_matrix(type = 'objective', quiet = F),0)
inq_tri <- data.frame(cc)
inq_tri <- inq_tri %>% mutate(trial = row_number()) %>% rowwise() %>% pivot_longer(cols = starts_with("X"), names_to = "RT") %>% mutate(RT = readr::parse_number(RT))
inq_tri <- inq_tri %>% arrange(trial,RT)
setwd('/Users/andypapale/clock2/Inquisit_design_files/')
write.csv(inq_tri,paste0('testing-Design-Matrix-with-Erasures-',as.character(df$iteration[i]),'.csv'))

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
d9_w_era <- read_csv(paste0('/Users/andypapale/clock2/Inquisit_design_files/testing-Design-Matrix-with-Erasures-',as.character(df$iteration[i]),'.csv'))


## get values without erasures
set.seed(df$iteration[i])
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
# set up contingency
bump_prominence <- 10
bump_value <- mean_val * bump_prominence
bump_center <- sample(seq(0, 2*pi, by = pi/20), 1, replace = FALSE)
contingency <- vm_circle_contingency(centers = c(centers, bump_center), weights = c(values, bump_value), widths = rep(width_sd, ncenters + 1), units = "radians")
qq <- troll_world$new(n_trials=ntrials, values=contingency$get_wfunc(), drift_sd=1)
qq$apply_flex(high_avg = 1, high_spread = 0, low_avg = df$low_avg[i], spread_max = 100, jump_high = T)
aa <- round(qq$get_values_matrix(type = 'objective', quiet = F),0)
inq_tri <- data.frame(aa)
inq_tri <- inq_tri %>% mutate(trial = row_number()) %>% rowwise() %>% pivot_longer(cols = starts_with("X"), names_to = "RT") %>% mutate(RT = readr::parse_number(RT))
inq_tri <- inq_tri %>% arrange(trial,RT)
setwd('/Users/andypapale/clock2/Inquisit_design_files/')
write.csv(inq_tri,paste0('testing-Design-Matrix-',as.character(df$iteration[i]),'.csv'))

d9 <- read_csv(paste0('/Users/andypapale/clock2/Inquisit_design_files/testing-Design-Matrix-',as.character(df$iteration[i]),'.csv'))

#setwd('/Users/andypapale/Desktop/2023-12-22-Testing/data')
# no clicks testing

df_inq <- read_csv('/Users/andypapale/Inquisit Code/EEG_clock/Clock_v2/data/clock_2_1_1_seed_152_testing_raw_1_2024-03-22-21-28-30-908.csv')
#df_inq <- read_csv('/Users/andypapale/Inquisit Code/EEG_clock/Clock_v2/data/clock_2_1_1_seed_152_testing_raw_1_2024-01-24-18-15-16-107.csv')
#df_inq <- read_csv('/Users/andypapale/Inquisit Code/EEG_clock/Clock_v2/data/clock_2_1_1_seed_152_testing_raw_1_2024-01-08-18-34-01-988.csv');
#df_inq <- read_csv('/Users/andypapale/Inquisit Code/EEG_clock/Clock_v2/data/clock_2_1_1_seed_6520_fMRI_raw_1_2024-01-10-16-55-04-979.csv')
#df_inq <- read_csv('/Users/andypapale/Inquisit Code/EEG_clock/Clock_v2/data/clock_2_1_1_seed_152_testing_raw_1_2024-01-23-21-25-55-768.csv')
#df_inq <- read_csv('/Users/andypapale/Inquisit Code/EEG_clock/Clock_v2/data/clock_2_1_1_seed_152_testing_raw_1_2024-01-24-18-03-24-721.csv')
#df <- read_csv('/Users/andypapale/Inquisit Code/EEG_clock/Clock_v2/data/clock_2_1_1_seed_152_testing_raw_1_2024-01-04-17-23-27-494.csv')
#df <- read_csv('/Users/andypapale/Inquisit Code/EEG_clock/Clock_v2/data/clock_2_1_1_seed_2534_testing_raw_1_2024-01-04-17-32-42-019.csv')
#df <- read_csv('/Users/andypapale/Inquisit Code/EEG_clock/Clock_v2/data/clock_2_1_1_seed_152_testing_raw_1_2024-01-02-23-57-16-667.csv')
#df <- read_csv('clock_2_1_1_seed_868_testing_raw_1_2023-12-22-20-58-38-919.csv')
#df <- read_csv('clock_2_1_1_seed_1464_testing_raw_1_2023-12-22-20-46-44-340.csv')
#df<- read_csv('clock_2_1_1_seed_1752_testing_raw_1_2023-12-22-20-26-53-565.csv')
#df <- read_csv('clock_2_1_1_seed_2534_testing_raw_1_2023-12-22-20-07-56-051.csv')
#df <- read_csv('clock_2_1_1_seed_4938_testing_raw_1_2023-12-22-19-48-02-863.csv')
#df <- read_csv('clock_2_1_1_seed_5094_testing_raw_1_2023-12-22-19-36-38-316.csv')
#df <- read_csv('clock_2_1_1_seed_5173_testing_raw_1_2023-12-22-19-20-17-939.csv')
#df <- read_csv('clock_2_1_1_seed_5815_testing_raw_1_2023-12-22-15-44-04-427.csv')
#df <- read_csv('clock_2_1_1_seed_6520_testing_raw_1_2023-12-22-18-49-18-373.csv')

df0 <- df_inq %>% filter(trialcode=="getValueVector")
#df0$trial <- df0$trial...18 # can make this step not necessary in final program
d9 <- read_csv(paste0('/Users/andypapale/clock2/Inquisit_design_files/Design-Matrix-',as.character(df$iteration),'.csv'))

d9_w_era <- read_csv(paste0('/Users/andypapale/clock2/Inquisit_design_files/Design-Matrix-with-Erasures-',as.character(df$iteration),'.csv'))

if (!file.exists(paste0('/Users/andypapale/Desktop/seed-',as.character(df$iteration),'-validation'))){
  dir.create(paste0('/Users/andypapale/Desktop/seed-',as.character(df$iteration),'-validation'))
}
setwd(paste0('/Users/andypapale/Desktop/seed-',as.character(df$iteration),'-validation'))

columns = c('e_pos','e_val');

val_cat <- NULL;
oval_cat <- NULL;
loc_cat <- NULL;
e_trial <- NULL;
for (iT in 1:300){
  e_val <- d9_w_era %>% filter(trial==iT);
  e_val <- e_val$value[round(e_val$RT)==round(era_loc[iT])];
  df_era <- data.frame(e_pos = era_loc[iT],e_val = e_val)
  inq_e_val = df0 %>% filter(trial==iT); 
  oeval <- d9 %>% filter(trial==iT);
  oeval <- oeval$value[round(oeval$RT)==round(era_loc[iT])];
  df_era <- data.frame(matrix(nrow = 0, ncol = length(columns)))
  colnames(df_era) = columns
  if (tt$erasure_segments$trial_type[iT]=='erasure'){
    era_loc2 <- inq_e_val$erasure_RT[1];
    e_val_inq <- inq_e_val$erasure_value[1];
    df_era <- data.frame(e_pos = era_loc[iT],e_val = e_val)
    val_cat <- rbind(val_cat,df_era$e_val)
    oval_cat <- rbind(oval_cat, oeval)
    loc_cat <- rbind(loc_cat,era_loc[iT])
    e_trial <- rbind(e_trial,iT)
  } else if (tt$erasure_segments$trial_type[iT]=='attention'){
    era_loc2 <- inq_e_val$stim_center_deg[1];
    e_val_inq <- NULL;
  } else if (tt$erasure_segments$trial_type[iT]=='no erasure'){
    era_loc2 <- NULL;
    e_val_inq <- NULL;
  }
  pdf(paste0('seed-',as.character(df$iteration),'-trial-',iT,'.pdf'),height=12,width=12);
  if (nrow(df_era)==1){
    gg1 <- ggplot(d9 %>% filter(trial==iT), aes(x=RT,y=value), color='black') +
     geom_line() + geom_point(size=5,color='black') +
     geom_line(data = df0 %>% filter(trial==iT),aes(x=x_final,y=y_final),color='red') +
     geom_point(data = df0 %>% filter(trial==iT),aes(x=x_final,y=y_final),color='red') +
     ylim(0,200) + geom_point(data = df_era, aes(x=e_pos,y=e_val),color='black',shape=2,size=10) +
     ggtitle(paste0('trial ',iT,' trial_type = ', tt$erasure_segments$trial_type[iT], ' trial_type_inq = ',inq_e_val$trial_type[1],' e_val = ',e_val, ' inq_e_val =', e_val_inq, ' e_max = ',round(tt$erasure_segments$segment_max[iT]*180/pi), ' inq_e_loc = ', era_loc2));
  } else {
    gg1 <-  ggplot(d9 %>% filter(trial==iT), aes(x=RT,y=value), color='black') +
      geom_line() + geom_point(size=5,color='black') +
      geom_line(data = df0 %>% filter(trial==iT),aes(x=x_final,y=y_final),color='red') +
      geom_point(data = df0 %>% filter(trial==iT),aes(x=x_final,y=y_final),color='red') +
      ylim(0,200) +
      ggtitle(paste0('trial ',iT,' trial_type = ', tt$erasure_segments$trial_type[iT], ' trial_type_inq = ',inq_e_val$trial_type[1]));
  }
      # gg1 <- ggplot(d9 %>% filter(trial==iT), aes(x=RT,y=value)) + geom_line() + geom_point(size=5,color='black') +  ylim(0,200) +
  #   ggtitle(paste0('trial ',iT,' trial_type = ', bb$erasure_segments$trial_type[iT], ' e_val = ',e_val, ' e_max = ',round(tt$erasure_segments$segment_max[iT]*180/pi)))
  print(gg1);
  dev.off()
}

setwd('~/clock2/Inquisit_design_files/')
write.table(val_cat,paste0('erasure-values-',as.character(df$iteration),'.csv'),row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(loc_cat,paste0('erasure-locations-',as.character(df$iteration),'.csv'),row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(oval_cat,paste0('original-values-',as.character(df$iteration),'.csv'),row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(e_trial,paste0('erasure-trial-',as.character(df$iteration),'.csv'),row.names = FALSE, col.names = FALSE, quote = FALSE)