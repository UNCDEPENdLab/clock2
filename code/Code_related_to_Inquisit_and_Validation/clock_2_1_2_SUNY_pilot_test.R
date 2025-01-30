
library(tidyverse)
library(data.table)
library(lme4)
library(BAMBI)
library(pracma)

#######################
## read design files ##
#######################

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

#####################
### read raw data ###
#####################

dfa <- read.table('/Users/andypapale/Downloads/Spaceship Behavioral Files/clock_2_1_2_seed_1752_EEG_raw_3883_2025-01-22-17-00-36-698.iqdat', sep = '\t',header = TRUE, fileEncoding = "UTF-8")
dfb <- read.table('/Users/andypapale/Downloads/Spaceship Behavioral Files/clock_2_1_2_seed_2534_EEG_raw_0085_2025-01-21-20-21-58-920.iqdat', sep = '\t',header = TRUE, fileEncoding = "UTF-8")
dfc <- read.table('/Users/andypapale/Downloads/Spaceship Behavioral Files/clock_2_1_2_seed_5094_EEG_raw_0019_2025-01-13-17-01-17-514.iqdat', sep = '\t',header = TRUE, fileEncoding = "UTF-8")
dfd <- read.table('/Users/andypapale/Downloads/Spaceship Behavioral Files/clock_2_1_2_seed_5094_EEG_raw_0050_2025-01-15-16-33-33-620.iqdat', sep = '\t',header = TRUE, fileEncoding = "UTF-8")
dfcomb <- rbind(dfa,dfb,dfc,dfd)


dfcomb <- inner_join(dfcomb,design,by=c('trial','seed'))

dfcomb <- dfcomb %>% filter(trialcode == 'feedback')
dfcomb <- dfcomb %>% mutate(u_present = case_when(trial_type == 'erasure' ~ TRUE,
                                            trial_type != 'erasure' ~ FALSE),
                      att_present = case_when(trial_type == 'attention' ~ TRUE,
                                              trial_type != 'attention' ~ FALSE),
                      u_location = case_when(trial_type == 'erasure' ~ stim_center_deg),
                      att_location = case_when(trial_type == 'attention' ~ stim_center_deg)
)
#df0 <- df0 %>% mutate(rt_index1 = case_when(rt_index+90 <= 360 ~ rt_index+90, rt_index+90 > 360 ~ abs(360-(rt_index+90))))
dfcomb <- dfcomb %>% mutate(resp_theta = case_when(pos_shifted < 90 ~ (pos_shifted + 270) * pi/180,
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
dfcomb$pre_clock_freeze = dfcomb$list.preClockFreeze.currentvalue

dfcomb <- dfcomb %>% arrange(subject,block,trial)

dfcomb <- dfcomb %>% mutate(rt_index_theta_c = zero_to_2pi(rt_index*pi/180))  %>% 
  group_by(subject) %>% 
  arrange(trial) %>%
  mutate(rt_index_theta_lag = lag(rt_index_theta_c)) %>% ungroup()

pdf('resp-RTvmax-SUNY-pilot.pdf',height=12,width=12)
dfcomb1 <- dfcomb %>% filter(mag > 100)
gg1 <- ggplot(dfcomb1, aes(x=trial,y=minuspi_to_pi(rt_index_theta_c - vmax_theta_c),color=mag)) + geom_point() + facet_wrap(~subject) +ylim(-3,3) + ggtitle('RT minus vmax_location')
print(gg1)
dev.off()

m1 <- lm(resp_theta_c ~ scale(vmax_location)*scale(vmax) + scale(resp_theta_c_lag)*outcome_lag, dfcomb)
summary(m1)
