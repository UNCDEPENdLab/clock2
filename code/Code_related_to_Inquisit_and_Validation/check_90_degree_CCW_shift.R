# 2025-02-03 AndyP
# Check 90 degree CCW shift

library(tidyverse)
library(data.table)
library(lme4)
library(BAMBI)
library(pracma)

design_file <- '/Users/andypapale/clock2/Inquisit_design_files/Design-Matrix-152.csv'
d1 <- read_csv(design_file)
#dx <- NULL; for (iT in 1:300){ dx <- rbind(dx, as.matrix(circshift(d1$value[d1$trial==iT],90)))}
#d1$value = dx;
d1 <- d1 %>% group_by(trial) %>% summarise(vmax = max(value),
                                           vmax_location = RT[which.max(value)]) %>% mutate(seed = 152)

design_file <- '/Users/andypapale/clock2/Inquisit_design_files/Design-Matrix-868.csv'
d2 <- read_csv(design_file)
#dx <- NULL; for (iT in 1:300){ dx <- rbind(dx, as.matrix(circshift(d2$value[d2$trial==iT],90)))}
#d2$value = dx;
d2 <- d2 %>% group_by(trial) %>% summarise(vmax = max(value),
                                           vmax_location = RT[which.max(value)]) %>% mutate(seed = 868)
design_file <- '/Users/andypapale/clock2/Inquisit_design_files/Design-Matrix-1464.csv'
d3 <- read_csv(design_file)
#dx <- NULL; for (iT in 1:300){ dx <- rbind(dx, as.matrix(circshift(d3$value[d3$trial==iT],90)))}
#d3$value = dx;
d3 <- d3 %>% group_by(trial) %>% summarise(vmax = max(value),
                                           vmax_location = RT[which.max(value)]) %>% mutate(seed = 1464)
design_file <- '/Users/andypapale/clock2/Inquisit_design_files/Design-Matrix-1752.csv'
d4 <- read_csv(design_file)
#dx <- NULL; for (iT in 1:300){ dx <- rbind(dx, as.matrix(circshift(d4$value[d4$trial==iT],90)))}
#d4$value = dx;
d4 <- d4 %>% group_by(trial) %>% summarise(vmax = max(value),
                                           vmax_location = RT[which.max(value)]) %>% mutate(seed = 1752)
design_file <- '/Users/andypapale/clock2/Inquisit_design_files/Design-Matrix-2534.csv'
d5 <- read_csv(design_file)
#dx <- NULL; for (iT in 1:300){ dx <- rbind(dx, as.matrix(circshift(d5$value[d5$trial==iT],90)))}
#d5$value = dx;
d5 <- d5 %>% group_by(trial) %>% summarise(vmax = max(value),
                                           vmax_location = RT[which.max(value)]) %>% mutate(seed = 2534)
design_file <- '/Users/andypapale/clock2/Inquisit_design_files/Design-Matrix-4938.csv'
d6 <- read_csv(design_file)
#dx <- NULL; for (iT in 1:300){ dx <- rbind(dx, as.matrix(circshift(d6$value[d6$trial==iT],90)))}
#d6$value = dx;
d6 <- d6 %>% group_by(trial) %>% summarise(vmax = max(value),
                                           vmax_location = RT[which.max(value)]) %>% mutate(seed = 4938)
design_file <- '/Users/andypapale/clock2/Inquisit_design_files/Design-Matrix-5094.csv'
d7 <- read_csv(design_file)
#dx <- NULL; for (iT in 1:300){ dx <- rbind(dx, as.matrix(circshift(d7$value[d7$trial==iT],90)))}
#d7$value = dx;
d7 <- d7 %>% group_by(trial) %>% summarise(vmax = max(value),
                                           vmax_location = RT[which.max(value)]) %>% mutate(seed = 5094)
design_file <- '/Users/andypapale/clock2/Inquisit_design_files/Design-Matrix-5173.csv'
d8 <- read_csv(design_file)
#dx <- NULL; for (iT in 1:300){ dx <- rbind(dx, as.matrix(circshift(d8$value[d8$trial==iT],90)))}
#d8$value = dx;
d8 <- d8 %>% group_by(trial) %>% summarise(vmax = max(value),
                                           vmax_location = RT[which.max(value)]) %>% mutate(seed = 5173)
design_file <- '/Users/andypapale/clock2/Inquisit_design_files/Design-Matrix-5815.csv'
d9 <- read_csv(design_file)
#dx <- NULL; for (iT in 1:300){ dx <- rbind(dx, as.matrix(circshift(d9$value[d9$trial==iT],90)))}
#d9$value = dx;
d9 <- d9 %>% group_by(trial) %>% summarise(vmax = max(value),
                                           vmax_location = RT[which.max(value)]) %>% mutate(seed = 5815)
design_file <- '/Users/andypapale/clock2/Inquisit_design_files/Design-Matrix-6520.csv'
d10 <- read_csv(design_file)
#dx <- NULL; for (iT in 1:300){ dx <- rbind(dx, as.matrix(circshift(d10$value[d10$trial==iT],90)))}
#d10$value = dx;
d10 <- d10 %>% group_by(trial) %>% summarise(vmax = max(value),
                                             vmax_location = RT[which.max(value)]) %>% mutate(seed = 6520)
design <- rbind(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10)

epoch_file <- '/Users/andypapale/clock2/Inquisit_design_files/epoch-152.csv'
e1 <- read_csv(epoch_file)
e1 <- e1 %>% as_tibble() %>% select(-`...1`) %>% mutate(trial = row_number()) %>% rename(epoch = tt.epoch)
e1 <- e1 %>% mutate(epoch_bin = rep(1:60, each=5),seed = 152) # just do every 10 trials for now

epoch_file <- '/Users/andypapale/clock2/Inquisit_design_files/epoch-868.csv'
e2 <- read_csv(epoch_file)
e2 <- e2 %>% as_tibble() %>% select(-`...1`) %>% mutate(trial = row_number()) %>% rename(epoch = tt.epoch)
e2 <- e2 %>% mutate(epoch_bin = rep(1:60, each=5),seed = 868) # just do every 10 trials for now

epoch_file <- '/Users/andypapale/clock2/Inquisit_design_files/epoch-1464.csv'
e3 <- read_csv(epoch_file)
e3 <- e3 %>% as_tibble() %>% select(-`...1`) %>% mutate(trial = row_number()) %>% rename(epoch = tt.epoch)
e3 <- e3 %>% mutate(epoch_bin = rep(1:60, each=5),seed = 1464) # just do every 10 trials for now

epoch_file <- '/Users/andypapale/clock2/Inquisit_design_files/epoch-1752.csv'
e4 <- read_csv(epoch_file)
e4 <- e4 %>% as_tibble() %>% select(-`...1`) %>% mutate(trial = row_number()) %>% rename(epoch = tt.epoch)
e4 <- e4 %>% mutate(epoch_bin = rep(1:60, each=5),seed = 1752) # just do every 10 trials for now

epoch_file <- '/Users/andypapale/clock2/Inquisit_design_files/epoch-2534.csv'
e5 <- read_csv(epoch_file)
e5 <- e5 %>% as_tibble() %>% select(-`...1`) %>% mutate(trial = row_number()) %>% rename(epoch = tt.epoch)
e5 <- e5 %>% mutate(epoch_bin = rep(1:60, each=5),seed = 2534) # just do every 10 trials for now

epoch_file <- '/Users/andypapale/clock2/Inquisit_design_files/epoch-4938.csv'
e6 <- read_csv(epoch_file)
e6 <- e6 %>% as_tibble() %>% select(-`...1`) %>% mutate(trial = row_number()) %>% rename(epoch = tt.epoch)
e6 <- e6 %>% mutate(epoch_bin = rep(1:60, each=5),seed = 4938) # just do every 10 trials for now

epoch_file <- '/Users/andypapale/clock2/Inquisit_design_files/epoch-5094.csv'
e7 <- read_csv(epoch_file)
e7 <- e7 %>% as_tibble() %>% select(-`...1`) %>% mutate(trial = row_number()) %>% rename(epoch = tt.epoch)
e7 <- e7 %>% mutate(epoch_bin = rep(1:60, each=5),seed = 5094) # just do every 10 trials for now

epoch_file <- '/Users/andypapale/clock2/Inquisit_design_files/epoch-5173.csv'
e8 <- read_csv(epoch_file)
e8 <- e8 %>% as_tibble() %>% select(-`...1`) %>% mutate(trial = row_number()) %>% rename(epoch = tt.epoch)
e8 <- e8 %>% mutate(epoch_bin = rep(1:60, each=5),seed = 5173) # just do every 10 trials for now

epoch_file <- '/Users/andypapale/clock2/Inquisit_design_files/epoch-5815.csv'
e9 <- read_csv(epoch_file)
e9 <- e9 %>% as_tibble() %>% select(-`...1`) %>% mutate(trial = row_number()) %>% rename(epoch = tt.epoch)
e9 <- e9 %>% mutate(epoch_bin = rep(1:60, each=5),seed = 5815) # just do every 10 trials for now

epoch_file <- '/Users/andypapale/clock2/Inquisit_design_files/epoch-6520.csv'
e10 <- read_csv(epoch_file)
e10 <- e10 %>% as_tibble() %>% select(-`...1`) %>% mutate(trial = row_number()) %>% rename(epoch = tt.epoch)
e10 <- e10 %>% mutate(epoch_bin = rep(1:60, each=5),seed = 6520) # just do every 10 trials for now

epoch <- rbind(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10)

design <- inner_join(design, epoch,by = c('trial','seed'))

df0 <- read_csv('/Users/andypapale/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Documents - DNPLskinner/skinner/data/prolific/clock_v2.1_pilot/prolific_01-09/raw/2024-01-09-Raw-noValueVector.csv')

# these subjects have a data formatting problem
df0 <- df0 %>% filter(subject != '6449707327ff66156c264c6f' & subject!='6464b45cf8e8a0f06fe011b9')


df0 <- df0 %>% filter(trialcode == 'feedback')
df0 <- df0 %>% mutate(u_present = case_when(trial_type == 'erasure' ~ TRUE,
                                            trial_type != 'erasure' ~ FALSE),
                      att_present = case_when(trial_type == 'attention' ~ TRUE,
                                              trial_type != 'attention' ~ FALSE),
                      u_location = case_when(trial_type == 'erasure' ~ stim_center_deg),
                      att_location = case_when(trial_type == 'attention' ~ stim_center_deg)
)

#df0$trial <- df0$trial...18 # can make this step not necessary in final program

df0 <- inner_join(df0,design,by=c('trial','seed'))

#df0 <- df0 %>% mutate(rt_index1 = case_when(rt_index+90 <= 360 ~ rt_index+90, rt_index+90 > 360 ~ abs(360-(rt_index+90))))
df0 <- df0 %>% mutate(resp_theta = case_when(pos_shifted < 90 ~ (pos_shifted + 270) * pi/180,
                                             pos_shifted >= 90 ~ (pos_shifted - 90) * pi/180),
                      resp_theta_c = zero_to_2pi(resp_theta),
                      vmax_theta = vmax_location * pi/180,
                      vmax_theta_c = zero_to_2pi(vmax_theta),
                      u_theta = case_when(u_location < 90 ~ (u_location + 270) * pi/180,
                                          u_location >= 90 ~ (u_location - 90) * pi/180),
                      u_theta_c = zero_to_2pi(u_theta),
                      att_theta = case_when(att_location < 90 ~ (att_location + 270) * pi/180,
                                            att_location >= 90 ~ (att_location - 90) * pi/180),
                      att_theta_c = zero_to_2pi(att_theta)) %>% 
  group_by(subject) %>% 
  arrange(trial) %>%
  mutate(resp_theta_c_lag = lag(resp_theta_c),
         vmax_theta_c_lag = lag(vmax_theta_c),
         vmax_loc_scaled = scale(vmax_location),
         outcome_lag = lag(inc_rg),
         omission = 0.7 <= rng,
         omission_lag = lag(omission),
         reward_lag = !omission_lag) %>% ungroup()
df0$pre_clock_freeze = df0$list.preClockFreeze.currentvalue

df0 <- df0 %>% group_by(epoch_bin,seed) %>% mutate(vmax_loc_mean = mean(vmax_location), vmax_mean = mean(vmax), mean_u_loc = mean(u_location, na.rm=TRUE)) %>% ungroup()
df0 <- df0 %>% add_count(epoch_bin, epoch,seed) %>% group_by(epoch_bin,seed) %>% mutate(epoch_mode = epoch[n==max(n)][1]) %>% select(-n) %>% ungroup()

df0 <- df0 %>% arrange(subject,block,trial)

df0 <- df0 %>% mutate(rt_index_theta_c = zero_to_2pi(rt_index*pi/180))  %>% 
  group_by(subject) %>% 
  arrange(trial) %>%
  mutate(rt_index_theta_lag = lag(rt_index_theta_c)) %>% ungroup()


df2 <- read_csv('/Users/andypapale/Downloads/df2.csv')

design_file <- '/Users/andypapale/clock2/Inquisit_design_files/Design-Matrix-1752.csv'
d152 <- read_csv(design_file)
df_152 <- df2 %>% filter(seed=='1752')
df_152 <- df_152 %>% mutate(resp_theta_deg = case_when(pos_shifted < 90 ~ (pos_shifted + 270),
                                                 pos_shifted >= 90 ~ (pos_shifted - 90)))
df_152_correct <- df_152 %>% mutate(RT = round(resp_theta_deg))
df_152_correct <- inner_join(d152,df_152_correct,by=c('trial','RT'))

hist(df_152_correct$value-df_152_correct$mag,breaks=100)

# these mismatches are due to the decay back to baseline after an erasure
df_mismatch_152_correct = df_152_correct %>% filter(df_152_correct$value - df_152_correct$mag != 0)

hist(df_mismatch_152_correct$value - df_mismatch_152_correct$mag, breaks=100)

ggplot(df_mismatch_152_correct, aes(x=value,y=mag,color=as.factor(trial_since_erasure_gen<=10))) + geom_point()

df_152_incorrect <- df_152 %>% mutate(RT = round(pos_shifted))
df_152_incorrect <- inner_join(d152,df_152_incorrect,by=c('trial','RT'))

hist(df_152_incorrect$value-df_152_incorrect$mag,breaks=100)

df_mismatch_152_incorrect <- df_152_incorrect %>% filter(df_152_incorrect$value - df_152_incorrect$mag != 0)

hist(df_mismatch_152_incorrect$value - df_mismatch_152_incorrect$mag, breaks=100)

                  