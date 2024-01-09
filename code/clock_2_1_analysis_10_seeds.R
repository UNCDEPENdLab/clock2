# 2023-12-11 AndyP
# Clock 2.1.1 analysis of 10 seed Prolific pilot

# 2023-11-29 AndyP
# clock 2_1 analysis

library(tidyverse)
library(data.table)
library(lme4)
library(BAMBI)
library(pracma)
design_file <- '/Users/andypapale/clock2/Inquisit_design_files/Design-Matrix-152.csv'
d1 <- read_csv(design_file)
d1$value <- circshift(d1$value, 90)
d1 <- d1 %>% group_by(trial) %>% summarise(vmax = max(value),
                                                   vmax_location = RT[which.max(value)]) %>% mutate(seed = 152)

design_file <- '/Users/andypapale/clock2/Inquisit_design_files/Design-Matrix-868.csv'
d2 <- read_csv(design_file)
d2$value <- circshift(d2$value, 90)
d2 <- d2 %>% group_by(trial) %>% summarise(vmax = max(value),
                                           vmax_location = RT[which.max(value)]) %>% mutate(seed = 868)
design_file <- '/Users/andypapale/clock2/Inquisit_design_files/Design-Matrix-1464.csv'
d3 <- read_csv(design_file)
d3$value <- circshift(d3$value, 90)
d3 <- d3 %>% group_by(trial) %>% summarise(vmax = max(value),
                                           vmax_location = RT[which.max(value)]) %>% mutate(seed = 1464)
design_file <- '/Users/andypapale/clock2/Inquisit_design_files/Design-Matrix-1752.csv'
d4 <- read_csv(design_file)
d4$value <- circshift(d4$value, 90)
d4 <- d4 %>% group_by(trial) %>% summarise(vmax = max(value),
                                           vmax_location = RT[which.max(value)]) %>% mutate(seed = 1752)
design_file <- '/Users/andypapale/clock2/Inquisit_design_files/Design-Matrix-2534.csv'
d5 <- read_csv(design_file)
d5$value <- circshift(d5$value, 90)
d5 <- d5 %>% group_by(trial) %>% summarise(vmax = max(value),
                                           vmax_location = RT[which.max(value)]) %>% mutate(seed = 2534)
design_file <- '/Users/andypapale/clock2/Inquisit_design_files/Design-Matrix-4938.csv'
d6 <- read_csv(design_file)
d6$value <- circshift(d6$value, 90)
d6 <- d6 %>% group_by(trial) %>% summarise(vmax = max(value),
                                           vmax_location = RT[which.max(value)]) %>% mutate(seed = 4938)
design_file <- '/Users/andypapale/clock2/Inquisit_design_files/Design-Matrix-5094.csv'
d7 <- read_csv(design_file)
d7$value <- circshift(d7$value, 90)
d7 <- d7 %>% group_by(trial) %>% summarise(vmax = max(value),
                                           vmax_location = RT[which.max(value)]) %>% mutate(seed = 5094)
design_file <- '/Users/andypapale/clock2/Inquisit_design_files/Design-Matrix-5173.csv'
d8 <- read_csv(design_file)
d8$value <- circshift(d8$value, 90)
d8 <- d8 %>% group_by(trial) %>% summarise(vmax = max(value),
                                           vmax_location = RT[which.max(value)]) %>% mutate(seed = 5173)
design_file <- '/Users/andypapale/clock2/Inquisit_design_files/Design-Matrix-5815.csv'
d9 <- read_csv(design_file)
d9$value <- circshift(d9$value, 90)
d9 <- d9 %>% group_by(trial) %>% summarise(vmax = max(value),
                                           vmax_location = RT[which.max(value)]) %>% mutate(seed = 5815)
design_file <- '/Users/andypapale/clock2/Inquisit_design_files/Design-Matrix-6520.csv'
d10 <- read_csv(design_file)
d10$value <- circshift(d10$value, 90)
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

#df1 <- read_csv('/Users/andypapale/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Documents - DNPLskinner/skinner/data/prolific/clock_v2.1_pilot/prolific_12-15/raw/papalea_prosper_clock_2_1_1_seed_152_raw_2312151732.csv')
#df2 <- read_csv('/Users/andypapale/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Documents - DNPLskinner/skinner/data/prolific/clock_v2.1_pilot/prolific_12-15/raw/papalea_prosper_clock_2_1_1_seed_868_raw_2312151815.csv')
#df3 <- read_csv('/Users/andypapale/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Documents - DNPLskinner/skinner/data/prolific/clock_v2.1_pilot/prolific_12-15/raw/papalea_prosper_clock_2_1_1_seed_1464_raw_2312151733.csv')
#df4 <- read_csv('/Users/andypapale/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Documents - DNPLskinner/skinner/data/prolific/clock_v2.1_pilot/prolific_12-15/raw/papalea_prosper_clock_2_1_1_seed_1752_raw_2312151735.csv')
#df5 <- read_csv('/Users/andypapale/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Documents - DNPLskinner/skinner/data/prolific/clock_v2.1_pilot/prolific_12-15/raw/papalea_prosper_clock_2_1_1_seed_2534_raw_2312151919.csv')
#df6 <- read_csv('/Users/andypapale/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Documents - DNPLskinner/skinner/data/prolific/clock_v2.1_pilot/prolific_12-15/raw/papalea_prosper_clock_2_1_1_seed_2534_raw_2312151920.csv')
#df7 <- read_csv('/Users/andypapale/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Documents - DNPLskinner/skinner/data/prolific/clock_v2.1_pilot/prolific_12-15/raw/papalea_prosper_clock_2_1_1_seed_4938_raw_2312151904.csv')
#df8 <- read_csv('/Users/andypapale/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Documents - DNPLskinner/skinner/data/prolific/clock_v2.1_pilot/prolific_12-15/raw/papalea_prosper_clock_2_1_1_seed_5094_raw_2312151734.csv')
#df9 <- read_csv('/Users/andypapale/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Documents - DNPLskinner/skinner/data/prolific/clock_v2.1_pilot/prolific_12-15/raw/papalea_prosper_clock_2_1_1_seed_6520_raw_2312151734.csv')
#df10 <- read_csv('/Users/andypapale/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Documents - DNPLskinner/skinner/data/prolific/clock_v2.1_pilot/prolific_12-15/raw/papalea_prosper_clock_2_1_1_seed_5173_raw_2312151732.csv')
#df11 <- read_csv('/Users/andypapale/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Documents - DNPLskinner/skinner/data/prolific/clock_v2.1_pilot/prolific_12-15/raw/papalea_prosper_clock_2_1_1_seed_5815_raw_2312151729.csv')
#df12 <- read_csv('/Users/andypapale/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Documents - DNPLskinner/skinner/data/prolific/clock_v2.1_pilot/prolific_12-08/raw/papalea_prosper_clock_2_1_1_seed_5173_raw_2312081911.csv')
#df13 <- read_csv('/Users/andypapale/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Documents - DNPLskinner/skinner/data/prolific/clock_v2.1_pilot/prolific_12-08/raw/papalea_prosper_clock_2_1_1_seed_5173_raw_2312090301.csv')
#df14 <- read_csv('/Users/andypapale/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Documents - DNPLskinner/skinner/data/prolific/clock_v2.1_pilot/prolific_12-08/raw/papalea_prosper_clock_2_1_1_seed_5815_raw_2312081912.csv')
#df15 <- read_csv('/Users/andypapale/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Documents - DNPLskinner/skinner/data/prolific/clock_v2.1_pilot/prolific_12-08/raw/papalea_prosper_clock_2_1_1_seed_6520_raw_2312081906.csv')

#df0 <- rbind(df1,df2,df3,df4,df5,df6,df7,df8,df9,df10,df11,df12,df13,df14,df15)
# rm(df0)
# common_cols = intersect(colnames(df1),colnames(df2))
# df0 <- rbind(subset(df1,select = common_cols),subset(df2,select = common_cols))
# common_cols = intersect(colnames(df0),colnames(df3))
# df0 <- rbind(subset(df0,select = common_cols),subset(df3,select = common_cols))
# common_cols = intersect(colnames(df0),colnames(df4))
# df0 <- rbind(subset(df0,select = common_cols),subset(df4,select = common_cols))
# common_cols = intersect(colnames(df0),colnames(df5))
# df0 <- rbind(subset(df0,select = common_cols),subset(df5,select = common_cols))
# common_cols = intersect(colnames(df0),colnames(df6))
# df0 <- rbind(subset(df0,select = common_cols),subset(df6,select = common_cols))
# common_cols = intersect(colnames(df0),colnames(df7))
# df0 <- rbind(subset(df0,select = common_cols),subset(df7,select = common_cols))
# common_cols = intersect(colnames(df0),colnames(df8))
# df0 <- rbind(subset(df0,select = common_cols),subset(df8,select = common_cols))
# common_cols = intersect(colnames(df0),colnames(df9))
# df0 <- rbind(subset(df0,select = common_cols),subset(df9,select = common_cols))
# common_cols = intersect(colnames(df0),colnames(df10))
# df0 <- rbind(subset(df0,select = common_cols),subset(df10,select = common_cols))
# common_cols = intersect(colnames(df0),colnames(df11))
# df0 <- rbind(subset(df0,select = common_cols),subset(df11,select = common_cols))

#df0 <- read_csv('/Users/andypapale/Inquisit Code/EEG_clock/Clock_v2/data/clock_2_1_1_seed_5815_testing_raw_1_2023-12-19-18-37-19-595.csv')

df0 <- read_csv('/Users/andypapale/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Documents - DNPLskinner/skinner/data/prolific/clock_v2.1_pilot/prolific_01-09/raw/2023-01-09-Raw-noValueVector.csv')

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


df0 <- df0 %>% mutate(resp_theta = pos_shifted * pi/180,
                      resp_theta_c = zero_to_2pi(resp_theta),
                      vmax_theta = vmax_location * pi/180,
                      vmax_theta_c = zero_to_2pi(vmax_theta),
                      u_theta = u_location * pi/180,
                      u_theta_c = zero_to_2pi(u_theta),
                      att_theta = att_location * pi/180,
                      att_theta_c = zero_to_2pi(att_theta)) %>% group_by(subject) %>% arrange(trial) %>%
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

# these subjects have a data formatting problem
df0 <- df0 %>% filter(subject != '6449707327ff66156c264c6f' & subject!='6464b45cf8e8a0f06fe011b9')

pdf('resp-RTvmax.pdf',height=12,width=12)
gg1 <- ggplot(df0, aes(x=trial,y=minuspi_to_pi(resp_theta_c - vmax_theta_c),color=inc_rg)) + geom_point() + facet_wrap(~subject) + ggtitle('RT minus vmax_location')
print(gg1)
dev.off()
ggplot(df0, aes(x=trial,y=minuspi_to_pi(resp_theta_c - u_theta_c),color=subject,size=inc_rg)) + geom_point() + facet_wrap(~subject) + ggtitle('RT minus u_location')
ggplot(df0, aes(x=trial,y=minuspi_to_pi(resp_theta_c - att_theta_c),color=subject)) + geom_point() + facet_wrap(~subject) + ggtitle('RT minus att_location')



# by seed
pdf('resp-RTvmax_seed1752.pdf',height=12,width=12)
gg1 <- ggplot(df0 %>% filter(seed==1752), aes(x=pos_shifted,y=inc_rg)) + facet_wrap(~epoch_bin) + geom_line() + geom_point() + geom_point(aes(x=vmax_loc_mean,y= 150), shape = 2,color='green') + geom_point(aes(x=u_location,y=erasure_value),shape=3,color='blue')
print(gg1)
dev.off()

pdf('resp-RTvmax_seed5815.pdf',height=12,width=12)
gg1 <- ggplot(df0 %>% filter(seed==5815), aes(x=pos_shifted,y=inc_rg)) + facet_wrap(~epoch_bin) + geom_line() + geom_point() + geom_point(aes(x=vmax_loc_mean,y= 150), shape = 2,color='green') + geom_point(aes(x=u_location,y=erasure_value),shape=3,color='blue')
print(gg1)
dev.off()

pdf('resp-RTvmax_seed152.pdf',height=12,width=12)
gg1 <- ggplot(df0 %>% filter(seed==152), aes(x=pos_shifted,y=inc_rg)) + facet_wrap(~epoch_bin) + geom_line() + geom_point() + geom_point(aes(x=vmax_loc_mean,y= 150), shape = 2,color='green')  + geom_point(aes(x=u_location,y=erasure_value),shape=3,color='blue')
print(gg1)
dev.off()

pdf('resp-RTvmax_seed868.pdf',height=12,width=12)
gg1 <- ggplot(df0 %>% filter(seed==868), aes(x=pos_shifted,y=inc_rg)) + facet_wrap(~epoch_bin) + geom_line() + geom_point() + geom_point(aes(x=vmax_loc_mean,y= 150), shape = 2,color='green') + geom_point(aes(x=u_location,y=erasure_value),shape=3,color='blue')
print(gg1)
dev.off()

pdf('resp-RTvmax_seed5094.pdf',height=12,width=12)
gg1 <- ggplot(df0 %>% filter(seed==5094), aes(x=pos_shifted,y=inc_rg)) + facet_wrap(~epoch_bin) + geom_line() + geom_point() + geom_point(aes(x=vmax_loc_mean,y= 150), shape = 2,color='green')  + geom_point(aes(x=u_location,y=erasure_value),shape=3,color='blue')
print(gg1)
dev.off()

pdf('resp-RTvmax_seed5173.pdf',height=12,width=12)
gg1 <- ggplot(df0 %>% filter(seed==5173), aes(x=pos_shifted,y=inc_rg)) + facet_wrap(~epoch_bin) + geom_line() + geom_point() + geom_point(aes(x=vmax_loc_mean,y= 150), shape = 2,color='green')  + geom_point(aes(x=u_location,y=erasure_value),shape=3,color='blue')
print(gg1)
dev.off()

pdf('resp-RTvmax_seed2534.pdf',height=12,width=12)
gg1 <- ggplot(df0 %>% filter(seed==2534), aes(x=pos_shifted,y=inc_rg)) + facet_wrap(~epoch_bin) + geom_line() + geom_point() + geom_point(aes(x=vmax_loc_mean,y= 150), shape = 2,color='green')  + geom_point(aes(x=u_location,y=erasure_value),shape=3,color='blue')
print(gg1)
dev.off()

pdf('resp-RTvmax_seed1464.pdf',height=12,width=12)
gg1 <- ggplot(df0 %>% filter(seed==1464), aes(x=pos_shifted,y=inc_rg)) + facet_wrap(~epoch_bin) + geom_line() + geom_point() + geom_point(aes(x=vmax_loc_mean,y= 150), shape = 2,color='green')  + geom_point(aes(x=u_location,y=erasure_value),shape=3,color='blue')
print(gg1)
dev.off()

pdf('resp-RTvmax_seed6520.pdf',height=12,width=12)
gg1 <- ggplot(df0 %>% filter(seed==6520), aes(x=pos_shifted,y=inc_rg)) + facet_wrap(~epoch_bin) + geom_line() + geom_point() + geom_point(aes(x=vmax_loc_mean,y= 150), shape = 2,color='green')  + geom_point(aes(x=u_location,y=erasure_value),shape=3,color='blue')
print(gg1)
dev.off()

pdf('resp-RTvmax_seed4938.pdf',height=12,width=12)
gg1 <- ggplot(df0 %>% filter(seed==4938), aes(x=pos_shifted,y=inc_rg)) + facet_wrap(~epoch_bin) + geom_line() + geom_point() + geom_point(aes(x=vmax_loc_mean,y= 150), shape = 2,color='green')  + geom_point(aes(x=u_location,y=erasure_value),shape=3,color='blue')
print(gg1)
dev.off()


ggplot(df0, aes(trial, vmax_theta_c, color = vmax)) + geom_line() + scale_color_viridis_c() + facet_wrap(~seed)
ggplot(design, aes(trial, vmax)) + geom_line() + scale_color_viridis_c()
ggplot(df0,aes(x=trial,y=minuspi_to_pi(resp_theta_c - vmax_theta_c))) + geom_point() + facet_grid(epoch~trial_type) + geom_hline(yintercept = 0) + ggtitle('RT minus vmax_location')
ggplot(df0,aes(x=trial,y=minuspi_to_pi(resp_theta_c - u_theta_c))) + geom_point() + facet_grid(epoch~trial_type) + geom_hline(yintercept = 0) + ggtitle('RT - u_location')

m1 <- lmerTest::lmer(pos_shifted ~ scale(vmax_location)*scale(vmax) + scale(resp_theta_c_lag)*outcome_lag + (1|subject), df0)
summary(m1)

df <- df0
mlist <- list()
for (i in 1:1000) {
  df$u_location[!df$u_present] <- runif(length(df$u_location[!df$u_present]), 0, 360)
  df$att_location[!df$att_present] <- runif(length(df$att_location[!df$att_present]), 0, 360)
  mi <- lmerTest::lmer(pos_shifted ~ scale(vmax_location)*scale(vmax) + 
                         scale(u_location)*u_present + 
                         scale(att_location)*att_present +
                         resp_theta_c_lag*outcome_lag +
                         (1|subject), 
                       df)
  mdf <- broom.mixed::tidy(mi)
  mdf$i <- i
  mlist[[i]] <- mdf
  print(i)
}
ddf <- rbindlist(mlist)
str(ddf)
mean_ddf <- ddf %>% filter(effect == "fixed") %>% select(-i, -effect, -group) %>% group_by(term) %>% summarise_all(mean)

