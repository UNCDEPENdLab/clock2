library(tidyverse)
library(BAMBI)

loc_1 <- read.csv('/Users/andypapale/clock2/Inquisit_design_files/erasure-locations-152.csv',header=FALSE) %>% mutate(seed=152,type = 'location')
loc_2 <- read.csv('/Users/andypapale/clock2/Inquisit_design_files/erasure-locations-868.csv',header=FALSE) %>% mutate(seed=868,type = 'location')
loc_3 <- read.csv('/Users/andypapale/clock2/Inquisit_design_files/erasure-locations-1464.csv',header=FALSE) %>% mutate(seed=1464,type = 'location')
loc_4 <- read.csv('/Users/andypapale/clock2/Inquisit_design_files/erasure-locations-1752.csv',header=FALSE) %>% mutate(seed=1752,type = 'location')
loc_5 <- read.csv('/Users/andypapale/clock2/Inquisit_design_files/erasure-locations-2534.csv',header=FALSE) %>% mutate(seed=2534,type = 'location')
loc_6 <- read.csv('/Users/andypapale/clock2/Inquisit_design_files/erasure-locations-4938.csv',header=FALSE) %>% mutate(seed=4938,type = 'location')
loc_7 <- read.csv('/Users/andypapale/clock2/Inquisit_design_files/erasure-locations-5094.csv',header=FALSE) %>% mutate(seed=5094,type = 'location')
loc_8 <- read.csv('/Users/andypapale/clock2/Inquisit_design_files/erasure-locations-5173.csv',header=FALSE) %>% mutate(seed=1573,type = 'location')
loc_9 <- read.csv('/Users/andypapale/clock2/Inquisit_design_files/erasure-locations-5815.csv',header=FALSE) %>% mutate(seed=1585,type = 'location')
loc_10 <- read.csv('/Users/andypapale/clock2/Inquisit_design_files/erasure-locations-6520.csv',header=FALSE) %>% mutate(seed=6520,type = 'location')

loc <- rbind(loc_1,loc_2,loc_3,loc_4,loc_5,loc_6,loc_7,loc_8,loc_9,loc_10)

ev_1 <- read.csv('/Users/andypapale/clock2/Inquisit_design_files/erasure-values-152.csv',header=FALSE) %>% mutate(seed=152)
ev_2 <- read.csv('/Users/andypapale/clock2/Inquisit_design_files/erasure-values-868.csv',header=FALSE) %>% mutate(seed=868)
ev_3 <- read.csv('/Users/andypapale/clock2/Inquisit_design_files/erasure-values-1464.csv',header=FALSE) %>% mutate(seed=1464)
ev_4 <- read.csv('/Users/andypapale/clock2/Inquisit_design_files/erasure-values-1752.csv',header=FALSE) %>% mutate(seed=1752)
ev_5 <- read.csv('/Users/andypapale/clock2/Inquisit_design_files/erasure-values-2534.csv',header=FALSE) %>% mutate(seed=2534)
ev_6 <- read.csv('/Users/andypapale/clock2/Inquisit_design_files/erasure-values-4938.csv',header=FALSE) %>% mutate(seed=4938)
ev_7 <- read.csv('/Users/andypapale/clock2/Inquisit_design_files/erasure-values-5094.csv',header=FALSE) %>% mutate(seed=5094)
ev_8 <- read.csv('/Users/andypapale/clock2/Inquisit_design_files/erasure-values-5173.csv',header=FALSE) %>% mutate(seed=1573)
ev_9 <- read.csv('/Users/andypapale/clock2/Inquisit_design_files/erasure-values-5815.csv',header=FALSE) %>% mutate(seed=1585)
ev_10 <- read.csv('/Users/andypapale/clock2/Inquisit_design_files/erasure-values-6520.csv',header=FALSE) %>% mutate(seed=6520)

ev <- rbind(ev_1,ev_2,ev_3,ev_4,ev_5,ev_6,ev_7,ev_8,ev_9,ev_10)
ev <- ev %>% rename(erasure_value = V1)

ov_1 <- read.csv('/Users/andypapale/clock2/Inquisit_design_files/original-values-152.csv',header=FALSE)
ov_2 <- read.csv('/Users/andypapale/clock2/Inquisit_design_files/original-values-868.csv',header=FALSE)
ov_3 <- read.csv('/Users/andypapale/clock2/Inquisit_design_files/original-values-1464.csv',header=FALSE)
ov_4 <- read.csv('/Users/andypapale/clock2/Inquisit_design_files/original-values-1752.csv',header=FALSE)
ov_5 <- read.csv('/Users/andypapale/clock2/Inquisit_design_files/original-values-2534.csv',header=FALSE)
ov_6 <- read.csv('/Users/andypapale/clock2/Inquisit_design_files/original-values-4938.csv',header=FALSE)
ov_7 <- read.csv('/Users/andypapale/clock2/Inquisit_design_files/original-values-5094.csv',header=FALSE)
ov_8 <- read.csv('/Users/andypapale/clock2/Inquisit_design_files/original-values-5173.csv',header=FALSE)
ov_9 <- read.csv('/Users/andypapale/clock2/Inquisit_design_files/original-values-5815.csv',header=FALSE)
ov_10 <- read.csv('/Users/andypapale/clock2/Inquisit_design_files/original-values-6520.csv',header=FALSE) 

ov <- rbind(ov_1,ov_2,ov_3,ov_4,ov_5,ov_6,ov_7,ov_8,ov_9,ov_10)
ov <- ov %>% rename(original_value = V1)

design_file <- '/Users/andypapale/clock2/Inquisit_design_files/Design-Matrix-152.csv'
d1 <- read_csv(design_file)
d1 <- d1 %>% group_by(trial) %>% summarise(vmax = max(value),
                                           vmax_location = RT[which.max(value)]) %>% mutate(seed = 152)
design_file <- '/Users/andypapale/clock2/Inquisit_design_files/Design-Matrix-868.csv'
d2 <- read_csv(design_file)
d2 <- d2 %>% group_by(trial) %>% summarise(vmax = max(value),
                                           vmax_location = RT[which.max(value)]) %>% mutate(seed = 868)
design_file <- '/Users/andypapale/clock2/Inquisit_design_files/Design-Matrix-1464.csv'
d3 <- read_csv(design_file)
d3 <- d3 %>% group_by(trial) %>% summarise(vmax = max(value),
                                           vmax_location = RT[which.max(value)]) %>% mutate(seed = 1464)
design_file <- '/Users/andypapale/clock2/Inquisit_design_files/Design-Matrix-1752.csv'
d4 <- read_csv(design_file)
d4 <- d4 %>% group_by(trial) %>% summarise(vmax = max(value),
                                           vmax_location = RT[which.max(value)]) %>% mutate(seed = 1752)
design_file <- '/Users/andypapale/clock2/Inquisit_design_files/Design-Matrix-2534.csv'
d5 <- read_csv(design_file)
d5 <- d5 %>% group_by(trial) %>% summarise(vmax = max(value),
                                           vmax_location = RT[which.max(value)]) %>% mutate(seed = 2534)
design_file <- '/Users/andypapale/clock2/Inquisit_design_files/Design-Matrix-4938.csv'
d6 <- read_csv(design_file)
d6 <- d6 %>% group_by(trial) %>% summarise(vmax = max(value),
                                           vmax_location = RT[which.max(value)]) %>% mutate(seed = 4938)
design_file <- '/Users/andypapale/clock2/Inquisit_design_files/Design-Matrix-5094.csv'
d7 <- read_csv(design_file)
d7 <- d7 %>% group_by(trial) %>% summarise(vmax = max(value),
                                           vmax_location = RT[which.max(value)]) %>% mutate(seed = 5094)
design_file <- '/Users/andypapale/clock2/Inquisit_design_files/Design-Matrix-5173.csv'
d8 <- read_csv(design_file)
d8 <- d8 %>% group_by(trial) %>% summarise(vmax = max(value),
                                           vmax_location = RT[which.max(value)]) %>% mutate(seed = 5173)
design_file <- '/Users/andypapale/clock2/Inquisit_design_files/Design-Matrix-5815.csv'
d9 <- read_csv(design_file)
d9 <- d9 %>% group_by(trial) %>% summarise(vmax = max(value),
                                           vmax_location = RT[which.max(value)]) %>% mutate(seed = 5815)
design_file <- '/Users/andypapale/clock2/Inquisit_design_files/Design-Matrix-6520.csv'
d10 <- read_csv(design_file)
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

values <- cbind(ev,ov) %>% mutate(erasure_minus_original = erasure_value - original_value)


ggplot(values, aes(x=erasure_minus_original)) + geom_histogram(bins=50) + facet_wrap(~seed)




erasure_RT_file <- '/Users/andypapale/clock2/Inquisit_design_files/erasure-locations-152.csv'
e1 <- read.csv(erasure_RT_file)
e1 <- e1 %>% rename(era_loc = colnames(e1)[1]) %>% mutate(seed=152)

erasure_RT_file <- '/Users/andypapale/clock2/Inquisit_design_files/erasure-locations-868.csv'
e2 <- read.csv(erasure_RT_file)
e2 <- e2 %>% rename(era_loc = colnames(e2)[1]) %>% mutate(seed=868)

erasure_RT_file <- '/Users/andypapale/clock2/Inquisit_design_files/erasure-locations-1464.csv'
e3 <- read.csv(erasure_RT_file)
e3 <- e3 %>% rename(era_loc = colnames(e3)[1]) %>% mutate(seed=1464)

erasure_RT_file <- '/Users/andypapale/clock2/Inquisit_design_files/erasure-locations-1752.csv'
e4 <- read.csv(erasure_RT_file)
e4 <- e4 %>% rename(era_loc = colnames(e4)[1]) %>% mutate(seed=1752)

erasure_RT_file <- '/Users/andypapale/clock2/Inquisit_design_files/erasure-locations-2534.csv'
e5 <- read.csv(erasure_RT_file)
e5 <- e5 %>% rename(era_loc = colnames(e5)[1]) %>% mutate(seed=2534)

erasure_RT_file <- '/Users/andypapale/clock2/Inquisit_design_files/erasure-locations-4938.csv'
e6 <- read.csv(erasure_RT_file)
e6 <- e6 %>% rename(era_loc = colnames(e6)[1]) %>% mutate(seed=4938)

erasure_RT_file <- '/Users/andypapale/clock2/Inquisit_design_files/erasure-locations-5094.csv'
e7 <- read.csv(erasure_RT_file)
e7 <- e7 %>% rename(era_loc = colnames(e7)[1]) %>% mutate(seed=5094)

erasure_RT_file <- '/Users/andypapale/clock2/Inquisit_design_files/erasure-locations-5173.csv'
e8 <- read.csv(erasure_RT_file)
e8 <- e8 %>% rename(era_loc = colnames(e8)[1]) %>% mutate(seed=5173)

erasure_RT_file <- '/Users/andypapale/clock2/Inquisit_design_files/erasure-locations-5815.csv'
e9 <- read.csv(erasure_RT_file)
e9 <- e9 %>% rename(era_loc = colnames(e9)[1]) %>% mutate(seed=5815)

erasure_RT_file <- '/Users/andypapale/clock2/Inquisit_design_files/erasure-locations-6520.csv'
e10 <- read.csv(erasure_RT_file)
e10 <- e10 %>% rename(era_loc = colnames(e10)[1]) %>% mutate(seed=6520)

era_loc <- rbind(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10)

erasure_trial_file <- '/Users/andypapale/clock2/Inquisit_design_files/erasure-trial-152.csv'
e1 <- read.csv(erasure_trial_file)
e1 <- e1 %>% rename(era_trial = colnames(e1)[1]) 

erasure_trial_file <- '/Users/andypapale/clock2/Inquisit_design_files/erasure-trial-868.csv'
e2 <- read.csv(erasure_trial_file)
e2 <- e2 %>% rename(era_trial = colnames(e2)[1]) 

erasure_trial_file <- '/Users/andypapale/clock2/Inquisit_design_files/erasure-trial-1464.csv'
e3 <- read.csv(erasure_trial_file)
e3 <- e3 %>% rename(era_trial = colnames(e3)[1]) 

erasure_trial_file <- '/Users/andypapale/clock2/Inquisit_design_files/erasure-trial-1752.csv'
e4 <- read.csv(erasure_trial_file)
e4 <- e4 %>% rename(era_trial = colnames(e4)[1]) 

erasure_trial_file <- '/Users/andypapale/clock2/Inquisit_design_files/erasure-trial-2534.csv'
e5 <- read.csv(erasure_trial_file)
e5 <- e5 %>% rename(era_trial = colnames(e5)[1]) 

erasure_trial_file <- '/Users/andypapale/clock2/Inquisit_design_files/erasure-trial-4938.csv'
e6 <- read.csv(erasure_trial_file)
e6 <- e6 %>% rename(era_trial = colnames(e6)[1]) 

erasure_trial_file <- '/Users/andypapale/clock2/Inquisit_design_files/erasure-trial-5094.csv'
e7 <- read.csv(erasure_trial_file)
e7 <- e7 %>% rename(era_trial = colnames(e7)[1]) 

erasure_trial_file <- '/Users/andypapale/clock2/Inquisit_design_files/erasure-trial-5173.csv'
e8 <- read.csv(erasure_trial_file)
e8 <- e8 %>% rename(era_trial = colnames(e8)[1]) 

erasure_trial_file <- '/Users/andypapale/clock2/Inquisit_design_files/erasure-trial-5815.csv'
e9 <- read.csv(erasure_trial_file)
e9 <- e9 %>% rename(era_trial = colnames(e9)[1]) 

erasure_trial_file <- '/Users/andypapale/clock2/Inquisit_design_files/erasure-trial-6520.csv'
e10 <- read.csv(erasure_trial_file)
e10 <- e10 %>% rename(era_trial = colnames(e10)[1]) 

era_trial <- rbind(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10)

era <- cbind(era_trial,era_loc)
era <- era %>% rename(trial = era_trial)

design1 <- inner_join(design,era,by=c('seed','trial'))

design1 <- design1 %>% mutate(erasure_RT_minus_vmax_RT = zero_to_2pi(era_loc - vmax_location)*180/pi)
ggplot(design1, aes(x=erasure_RT_minus_vmax_RT)) + geom_histogram(bins=20) + facet_wrap(~seed)

sum <- design1 %>% group_by(seed) %>% summarize(eRTp15vmaxRT = sum(erasure_RT_minus_vmax_RT <= 15, na.rm=TRUE), eRTm15vmaxRT = sum(erasure_RT_minus_vmax_RT >=345, na.rm=TRUE)) %>% ungroup()
sum <- sum %>% mutate(eRTpm15vmaxRT = eRTp15vmaxRT + eRTm15vmaxRT) %>% select(seed,eRTpm15vmaxRT)
