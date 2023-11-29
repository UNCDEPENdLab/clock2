# 2023-11-29 AndyP
# clock 2_1 analysis

library(tidyverse)
library(data.table)
library(lme4)
library(BAMBI)
design_file <- '/Users/andypapale/clock2/Design-Matrix-2455.csv'
design <- read_csv(design_file)
design <- design %>% group_by(trial) %>% summarise(vmax = max(value),
                                            vmax_location = RT[which.max(value)])


epoch_file <- '/Users/andypapale/clock2/epoch-2455.csv'
epoch <- read_csv(epoch_file)
epoch <- epoch %>% as_tibble() %>% select(-`...1`) %>% mutate(trial = row_number()) %>% rename(epoch = qq.epoch)
epoch <- epoch %>% mutate(epoch_bin = rep(1:50, each=6)) # just do every 10 trials for now

design <- inner_join(design,epoch,by='trial')

df1 <- read_csv('/Users/andypapale/clock2/pilot_data/papalea_prosper_clock_2_1_1_raw_2311291846.csv')
df1 <- df1 %>% filter(trialcode == 'feedback')
df1 <- df1 %>% mutate(u_present = case_when(trial_type == 'erasure' ~ TRUE,
                                            trial_type != 'erasure' ~ FALSE),
                      att_present = case_when(trial_type == 'attention' ~ TRUE,
                                              trial_type != 'attention' ~ FALSE),
                      u_location = case_when(trial_type == 'erasure' ~ stim_center_deg),
                      att_location = case_when(trial_type == 'attention' ~ stim_center_deg)
)
                      
df1$trial <- df1$trial...18 # can make this step not necessary in final program

df1 <- inner_join(df1,design,by='trial')
df1 <- df1 %>% mutate(resp_theta = pos_shifted * pi/180,
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
df1$pre_clock_freeze = df1$list.preClockFreeze.currentvalue

df1 <- df1 %>% group_by(epoch_bin) %>% mutate(vmax_loc_mean = mean(vmax_location), vmax_mean = mean(vmax)) %>% ungroup()
df1 <- df1 %>% add_count(epoch_bin, epoch) %>% group_by(epoch_bin) %>% mutate(epoch_mode = epoch[n==max(n)][1]) %>% select(-n) %>% ungroup()


df1 <- df1 %>% arrange(subject,block,trial)

ggplot(df1, aes(x=trial,y=minuspi_to_pi(resp_theta_c - vmax_theta_c),color=subject)) + geom_point() + ggtitle('RT minus vmax_location')
ggplot(df1, aes(x=trial,y=minuspi_to_pi(resp_theta_c - u_theta_c),color=subject)) + geom_point() + ggtitle('RT minus u_location')
ggplot(df1, aes(x=trial,y=minuspi_to_pi(resp_theta_c - att_theta_c),color=subject)) + geom_point() + ggtitle('RT minus att_location')
ggplot(df1, aes(x=pos_shifted,y=inc_rg)) + facet_wrap(~epoch_bin) + geom_line() + geom_point(aes(color=omission)) + geom_point(aes(x=vmax_loc_mean,y= 150), shape = 2,color='green')
ggplot(df1, aes(trial, vmax_theta_c, color = vmax)) + geom_line() + scale_color_viridis_c()
ggplot(design, aes(trial, vmax)) + geom_line() + scale_color_viridis_c()
ggplot(df1,aes(x=trial,y=minuspi_to_pi(resp_theta_c - vmax_theta_c), color=subject)) + geom_point() + facet_grid(epoch~trial_type) + geom_hline(yintercept = 0) + ggtitle('RT minus vmax_location')
ggplot(df1,aes(x=trial,y=minuspi_to_pi(resp_theta_c - u_theta_c), color=subject)) + geom_point() + facet_grid(epoch~trial_type) + geom_hline(yintercept = 0) + ggtitle('RT - u_location')

m1 <- lmerTest::lmer(pos_shifted ~ scale(vmax_location)*scale(vmax) + scale(resp_theta_c_lag)*outcome_lag + (1|subject), df1)
summary(m1)

df <- df1
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
