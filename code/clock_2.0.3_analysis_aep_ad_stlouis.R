library(tidyverse)
library(reshape2)
library(lme4)
library(brms)
library(data.table)
# 2023-05-03 AndyP
# clock 2.0.1 analysis


# set paths (will be good to make this more user-general)
sys <- Sys.info()
if (sum(str_detect(Sys.info(), "Alex"))>1) { base_dir <- "~/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Documents/"
} else if (sum(str_detect(Sys.info(), "alexdombrovski"))>1) {
  repo_dir = "~/code/clock2/"
  base_dir <- '~/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Documents - DNPLskinner' }# on skinner}
output_dir <- "~/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Documents/skinner/data/prolific/clock_v2_pilot"
# design_file <- '~/code/clock2/2022-04-25-DesignFile.csv' # in Michael's clock2 repo
# design_file <- '/Volumes/Users/Andrew/2022-05-08-DesignFile.csv'

data_wave <- "Sept" # "Jul" before balanced U/A/0 + fixed magnitude, "Sept" with balancing + jittered magnitude

design_file <- '~/code/clock2/2022-05-08-DesignFile.csv'
design <- as.matrix(read_csv(design_file)) %>% 
  as_tibble() %>% select(-`...1`) %>% mutate(timepoint = row_number()) %>% rowwise() %>% pivot_longer(cols = starts_with("V"), names_to = "trial") %>%
  mutate(trial = extract_numeric(trial)) %>% group_by(trial) %>% summarise(vmax = max(value),
                                                                           vmax_location = timepoint[which.max(value)])
design <- design %>% mutate(epoch = case_when(vmax > 74 ~ 'High Value', vmax <= 74 ~ 'Low Value'))
design <- design %>% mutate(epoch_bin = as.factor(case_when(epoch=='High Value' & trial < 33 ~ 1,
                                                            epoch=='High Value' & trial > 33 & trial < 74 ~ 2,
                                                            epoch=='High Value' & trial > 74 & trial < 116 ~ 3,
                                                            epoch=='High Value' & trial > 116 & trial < 161 ~ 4,
                                                            epoch=='High Value' & trial > 161 & trial < 205 ~ 5,
                                                            epoch=='High Value' & trial > 205 & trial < 250 ~ 6,
                                                            epoch=='High Value' & trial > 250 & trial < 291 ~ 7,
                                                            epoch=='High Value' & trial > 291 ~ 8,
                                                            epoch=='Low Value' & trial < 74 ~ -1,
                                                            epoch=='Low Value' & trial >= 74 & trial < 116 ~ -2,
                                                            epoch=='Low Value' & trial >= 116 & trial < 161 ~ -3,
                                                            epoch=='Low Value' & trial >=161 & trial < 205 ~ -4,
                                                            epoch=='Low Value' & trial >=205 & trial < 250 ~ -5,
                                                            epoch=='Low Value' & trial >=250 & trial < 291 ~ -6)))
design <- design %>% group_by(epoch_bin) %>% mutate(vmax_loc_mean = mean(vmax_location), vmax_mean = mean(vmax)) %>% ungroup()
design <- design %>% mutate(vmax_loc_lag1 = case_when(
  epoch_bin==2 ~ mean(design$vmax_loc_mean[design$epoch_bin==1],na.rm=TRUE),
  epoch_bin==3 ~ mean(design$vmax_loc_mean[design$epoch_bin==2],na.rm=TRUE),
  epoch_bin==4 ~ mean(design$vmax_loc_mean[design$epoch_bin==3],na.rm=TRUE),
  epoch_bin==5 ~ mean(design$vmax_loc_mean[design$epoch_bin==4],na.rm=TRUE),
  epoch_bin==6 ~ mean(design$vmax_loc_mean[design$epoch_bin==5],na.rm=TRUE),
  epoch_bin==7 ~ mean(design$vmax_loc_mean[design$epoch_bin==6],na.rm=TRUE),
  epoch_bin==-1 ~ mean(design$vmax_loc_mean[design$epoch_bin==1],na.rm=TRUE),
  epoch_bin==-2 ~ mean(design$vmax_loc_mean[design$epoch_bin==2],na.rm=TRUE),
  epoch_bin==-3 ~ mean(design$vmax_loc_mean[design$epoch_bin==3],na.rm=TRUE),
  epoch_bin==-4 ~ mean(design$vmax_loc_mean[design$epoch_bin==4],na.rm=TRUE),
  epoch_bin==-5 ~ mean(design$vmax_loc_mean[design$epoch_bin==5],na.rm=TRUE),
  epoch_bin==-6 ~ mean(design$vmax_loc_mean[design$epoch_bin==6],na.rm=TRUE))
  ,
  vmax_loc_lag2 = case_when(
    epoch_bin==3 ~ mean(design$vmax_loc_mean[design$epoch_bin==1],na.rm=TRUE),
    epoch_bin==4 ~ mean(design$vmax_loc_mean[design$epoch_bin==2],na.rm=TRUE),
    epoch_bin==5 ~ mean(design$vmax_loc_mean[design$epoch_bin==3],na.rm=TRUE),
    epoch_bin==6 ~ mean(design$vmax_loc_mean[design$epoch_bin==4],na.rm=TRUE),
    epoch_bin==7 ~ mean(design$vmax_loc_mean[design$epoch_bin==5],na.rm=TRUE),
    epoch_bin==-2 ~ mean(design$vmax_loc_mean[design$epoch_bin==1],na.rm=TRUE),
    epoch_bin==-3 ~ mean(design$vmax_loc_mean[design$epoch_bin==2],na.rm=TRUE),
    epoch_bin==-4 ~ mean(design$vmax_loc_mean[design$epoch_bin==3],na.rm=TRUE),
    epoch_bin==-5 ~ mean(design$vmax_loc_mean[design$epoch_bin==4],na.rm=TRUE),
    epoch_bin==-6 ~ mean(design$vmax_loc_mean[design$epoch_bin==5],na.rm=TRUE))
)
# initial sanity checks on the disign file
#str(design)
#ggplot(design, aes(trial, vmax_location, color = vmax)) + geom_line() + scale_color_viridis_c()
#ggplot(design, aes(trial, vmax)) + geom_line() + scale_color_viridis_c()
#median_vmax = median(design$vmax)

if (data_wave == "Jul") {
  # combine behavioral data from two batches of 25 Prolific subjects
  df1_m <- read_csv('~/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Documents/skinner/data/prolific/clock_v2_pilot/pilot_v3_03_07-31-23/raw/papalea_prosper_eeg_clock_v2_0_3_raw_2308031715.csv')
  df2_m <- read_csv('~/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Documents/skinner/data/prolific/clock_v2_pilot/pilot_v3_03_07-28-23/raw/papalea_prosper_eeg_clock_v2_0_3_raw_2307282152.csv')
  df3_m <- read_csv('~/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Documents/skinner/data/prolific/clock_v2_pilot/pilot_v3_03_07-03-23/raw/papalea_prosper_eeg_clock_v2_0_3_raw_2307031911.csv')
  df1_m <- df1_m %>% select(!date)
  df2_m <- df2_m %>% select(!date)
  df3_m <- df3_m %>% select(!date)
  dfm <- rbind(df1_m,df2_m,df3_m) %>% filter(trialcode=='dispFeedback_noU' | trialcode=='dispFeedback_U') # the score is calculated here, the variable of interest is rt_shifted
  dfm <- dfm %>% mutate(experiment ='mushrooms')
  df1_t <- read_csv('~/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Documents/skinner/data/prolific/clock_v2_pilot/pilot_v3_Trolls_08-10-23/raw/papalea_prosper_eeg_clock_v2_0_3_trolls_raw_2308150248.csv')
  df2_t <- read_csv('~/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Documents/skinner/data/prolific/clock_v2_pilot/pilot_v3_Trolls_08-10-23/raw/papalea_prosper_eeg_clock_v2_0_3_trolls_raw_2308151454.csv')
  df3_t <- read_csv('~/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Documents/skinner/data/prolific/clock_v2_pilot/pilot_v3_Trolls_08-10-23/raw/papalea_prosper_eeg_clock_v2_0_3_trolls_raw_2308132214.csv')
  df4_t <- read_csv('~/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Documents/skinner/data/prolific/clock_v2_pilot/pilot_v3_Trolls_08-10-23/raw/papalea_prosper_eeg_clock_v2_0_3_trolls_raw_2308132213.csv')
  df5_t <- read_csv('~/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Documents/skinner/data/prolific/clock_v2_pilot/pilot_v3_Trolls_08-10-23/raw/papalea_prosper_eeg_clock_v2_0_3_trolls_raw_2308132210.csv')
  df1_t <- df1_t %>% select(!date)
  df2_t <- df2_t %>% select(!date)
  df3_t <- df3_t %>% select(!date)
  df4_t <- df4_t %>% select(!date)
  df5_t <- df5_t %>% select(!date)
  df1 <- rbind(df1_t,df2_t,df3_t,df4_t,df5_t) %>% filter(trialcode=='dispFeedback_noU' | trialcode=='dispFeedback_U') # the score is calculated here, the variable of interest is rt_shifted
  df1 <- df1 %>% mutate(experiment = 'trolls')
  df1 <- rbind(df1, dfm) } else if (data_wave == "Sept") {
    df1 <- read_csv(file.path(paste0(repo_dir, "/data/raw_data_good_9-22.csv"))) %>% filter(trialcode=='dispFeedback_noU' | trialcode=='dispFeedback_U')
  }
df1 <- df1 %>%  mutate(u_present = case_when(!is.na(windPos_out) & local_uncertainty == "wind" ~ TRUE,
                                             !is.na(fogPos_out) & local_uncertainty == "cloud" ~ TRUE,
                                             is.na(windPos_out) ~ FALSE,
                                             is.na(fogPos_out) ~ FALSE
),
att_present = case_when(!is.na(windPos_out) & local_uncertainty == "cloud" ~ TRUE,
                        !is.na(fogPos_out) & local_uncertainty == "wind" ~ TRUE,
                        is.na(windPos_out) ~ FALSE,
                        is.na(fogPos_out) ~ FALSE
),
u_location = case_when(
  !is.na(windPos_out) & local_uncertainty == "wind" ~ windPos_out,
  !is.na(fogPos_out) & local_uncertainty == "cloud" ~ fogPos_out
),
att_location = case_when(
  !is.na(windPos_out) & attentional_control == "wind" ~ windPos_out,
  !is.na(fogPos_out) & attentional_control == "cloud" ~ fogPos_out
)
)

df1 <- df1 %>% mutate(trial_to_merge = list.rt_0.currentindex)
design <- design %>% mutate(trial_to_merge = trial)

df1 <- full_join(df1,design,by='trial_to_merge')
df1 <- df1 %>% mutate(resp_theta = pos_shifted * pi/180,
                      resp_theta_c = resp_theta - pi,
                      vmax_theta = vmax_location * pi/180,
                      vmax_theta_c = vmax_theta - pi,
                      u_theta = u_location * pi/180,
                      u_theta_c = u_theta - pi,
                      att_theta = att_location * pi/180,
                      att_theta_c = att_theta - pi,
                      outcome_sc = as.vector(scale(as.numeric(Earnings)))) %>% group_by(subject) %>% arrange(meta_trialCount) %>%
  mutate(resp_theta_c_lag = lag(resp_theta_c),
         vmax_theta_c_lag = lag(vmax_theta_c),
         vmax_loc_scaled = as.vector(scale(vmax_location)),
         outcome_sc_lag = lag(outcome_sc),
         omission = Earnings< .01,
         omission_lag = lag(omission),
         reward_lag = !omission_lag) %>% ungroup()
df1$pre_clock_freeze = df1$list.preClockFreeze.currentvalue

df1 <- df1 %>% arrange(subject,blocknum,meta_trialCount)
# sanity checks on behavioral data: ensure we know where erasures/control stimuli were
#ggplot(df1, aes(trial, u_present, color = as.factor(chooseUncertainty))) + geom_point() + facet_wrap(~subject)
#pdf('test.pdf',height=18,width=16)
#gg1 <- ggplot(df1, aes(trial, local_uncertainty, color = as.factor(chooseFog1))) + geom_point(size=1) + facet_wrap(~subject) + geom_text(aes(x = 50, y = .5, label = local_uncertainty))
#print(gg1)
#dev.off()
# ggplot(df, aes(trial, uncertainty_block, color = as.factor(chooseWind2))) + geom_point() + facet_wrap(~subject) + geom_text(aes(x = 50, y = .5, label = local_uncertainty))
# ggplot(df, aes(trial, uncertainty_block, color = as.factor(uncertainty_present))) + geom_point() + facet_wrap(~subject) + geom_text(aes(x = 50, y = .5, label = local_uncertainty))
# 
# ggplot(df, aes(trial, uncertainty_block, color = as.factor(chooseAttentionalControl))) + geom_point() + facet_wrap(~subject)
# 
# ggplot(df, aes(trial, uncertainty_block, color = as.factor(chooseUncertainty))) + geom_jitter()
# 
# 
# ggplot(df, aes(vmax_location, chosen_location, alpha = vmax, color = subject)) + geom_point() + geom_smooth(method = "loess")
# ggplot(df, aes(uncertainty_location, chosen_location,  color = subject)) + geom_point() + geom_smooth(method = "loess")
df <- df1
df$u_location[!df$u_present] <- runif(length(df$u_location[!df$u_present]), 0, 360)
df$att_location[!df$att_present] <- runif(length(df$att_location[!df$att_present]), 0, 360)
df2 <- df %>% group_by(subject) %>% mutate(vmax_scaled = scale(vmax), u_loc_scaled = scale(u_location), att_loc_scaled = scale(att_location)) %>% ungroup()
df2 <- df2 %>% filter(!is.na(pos_shifted) & !is.na(vmax_loc_scaled) & !is.na(vmax_scaled) & !is.na(u_loc_scaled) & !is.na(u_present) & !is.na(resp_theta_c_lag) & !is.na(omission_lag) & !is.na(att_present) & !is.na(att_loc_scaled) & !is.na(subject))
setwd("~/code/clock2/data")
if (data_wave == "Jul") {df2 <- df2 %>% 
  select(experiment,pos_shifted,vmax_loc_scaled,vmax_scaled,u_loc_scaled,u_present,att_present,att_loc_scaled,resp_theta_c_lag,omission_lag,subject)
df2 <- df2 %>% arrange(subject)
df2 <- transform(df2, id=match(subject, unique(subject)))
df2$id <- as.numeric(df2$id)
df2$pos_shifted <- df2$pos_shifted * pi/180
# library(bpnreg)
# m5 <- bpnreg::bpnme(pred.I = pos_shifted ~ vmax_loc_scaled*vmax_scaled*experiment + u_loc_scaled*u_present*vmax_scaled + att_loc_scaled*att_present*vmax_scaled + resp_theta_c_lag*omission_lag + (1|id), data = df2, its = 2000, burn = 100, n.lag=3, seed=121)
df2 <- df2 %>% select(subject, experiment, outcome_sc, outcome_sc_lag, meta_trialCount, resp_theta_c, att_theta_c, u_theta_c,
                      u_present, u_location, att_location, att_present, pos_shifted,
                      vmax_theta, vmax_location, ev)
write_csv2(df1, file = "trolls_mushrooms_n202_Jul2023.csv")
# initial Gaussian glms: good for manipulation checks
# only the value bump
m1 <- lmer(pos_shifted ~ scale(vmax_location)*scale(vmax) + experiment + (1|subject), df1 )
summary(m1)
car::Anova(m1, '3')
# add the previous response/outcome (win-stay/lose-shift)
m1a_1 <- lmer(pos_shifted ~ scale(vmax_location)*scale(vmax)*experiment + resp_theta_c_lag * outcome_lag*experiment + (1|subject), df1 )
summary(m1a_1)
car::Anova(m1a_1, '3')
car::vif(m1a_1)
# add uncertainty
m2_1 <- lmerTest::lmer(pos_shifted ~ scale(vmax_location)*scale(vmax)*experiment + scale(u_location)*scale(vmax)*experiment + (1|subject), 
                       df1)
summary(m2_1)
car::Anova(m2_1, '3')
# full model with uncertainty and attentional control

m2_2 <- lmerTest::lmer(pos_shifted ~ scale(vmax_location)*scale(vmax) + scale(u_location)*scale(vmax) + 
                         resp_theta_c_lag*omission_lag + (1|subject), 
                       df1)
summary(m2_2)
car::Anova(m2_2, '3')

m2_3 <- lmerTest::lmer(pos_shifted ~ scale(vmax_location)*scale(vmax) + scale(att_location)*scale(vmax) + 
                         resp_theta_c_lag*omission_lag+ (1|subject), 
                       df1)
summary(m2_3)
car::Anova(m2_3, '3')

m5_1 <- lmerTest::lmer(pos_shifted ~ scale(vmax_location):pre_clock_freeze + scale(vmax_location):scale(vmax):pre_clock_freeze + 
                         scale(u_location):u_present:pre_clock_freeze + pre_clock_freeze*scale(att_location):att_present + 
                         resp_theta_c_lag*omission_lag*pre_clock_freeze + (1|subject), 
                       df)
summary(m5_1)
car::Anova(m5_1, '3')
} else if (
  data_wave == "Sept") {
  df2 <- df2 %>% arrange(subject)
  df2 <- transform(df2, id=match(subject, unique(subject)))
  df2$id <- as.numeric(df2$id)
  df2$pos_shifted <- df2$pos_shifted * pi/180
  df2 <- df2 %>% 
    select(id, pos_shifted,vmax_loc_scaled,vmax_scaled,u_loc_scaled,u_present,att_present,att_loc_scaled, resp_theta_c, att_theta_c, u_theta_c,
           resp_theta_c_lag,omission_lag,subject, meta_trialCount, chooseUncertainty, chooseAttentionalControl,
           Earnings, vmax_location) %>%
    rename(trial = meta_trialCount, score = Earnings)
  write_csv2(df1, file = "mushrooms_n180_Sept2023.csv")}
library(data.table)
df2 <- transform(df2, Counter = ave(chooseUncertainty, rleid(u_present, chooseUncertainty), FUN = seq_along)) %>% 
  group_by(subject) %>% arrange(subject, trial) %>% mutate(
    ntrial_u_sampled = case_when(
      chooseUncertainty == 0 ~ NA,
      chooseUncertainty == 1 & Counter > 1 ~ NA,
      chooseUncertainty == 1 & !lag(u_present) ~ 1,
      chooseUncertainty == 1 & lag(u_present) & Counter == 1 ~ lag(Counter) + 1
    )
  )
df2 <- df2 %>% arrange(subject, trial)

df2 <- transform(df2, clickCounter = ave(u_present, rleid(u_present), FUN = seq_along)) %>% 
  group_by(subject) %>% arrange(subject, trial) %>% mutate(
    clickCounter = case_when(
      !u_present ~ NA,
      T ~ clickCounter
    )
  )

sdf <- df2 %>% group_by(id) %>% summarize(u_choices = sum(chooseUncertainty),
                                          u_ntrials = sum(u_present),
                                          att_choices = sum(chooseAttentionalControl),
                                          att_ntrials = sum(att_present),
                                          u_percent_sampled = u_choices*100/u_ntrials,
                                          att_percent_sampled = att_choices*100/att_ntrials,
                                          u_att_choices_diff = u_choices - att_choices,
                                          u_att_percent_sampled_diff = u_percent_sampled - att_percent_sampled
) %>% ungroup() %>% mutate(
  u_rank = rank(u_att_percent_sampled_diff)
)
df2 <- df2 %>% inner_join(sdf, by = "id")
# note, 1% of pos_shifted have astronomical values, censor
df2$pos_shifted[df2$pos_shifted>360] <- NA

ggplot(sdf, aes(att_choices, u_choices)) + geom_point() + geom_abline(slope = 1)
ggplot(sdf, aes(att_percent_sampled, u_percent_sampled, color = u_rank)) + geom_point() + geom_abline(slope = 1)
ggplot(sdf, aes(att_ntrials, u_ntrials, color = u_rank)) + geom_point() + geom_abline(slope = 1)

ggplot(sdf, aes(u_percent_sampled)) + geom_histogram()
ggplot(sdf, aes(att_percent_sampled)) + geom_histogram()
ggplot(sdf, aes(u_att_percent_sampled_diff)) + geom_histogram()
# diagnostic plots to understand how subjects experienced conditions
ggplot(sdf, aes(u_percent_sampled, u_ntrials)) + geom_point() 
ggplot(sdf, aes(att_percent_sampled, att_ntrials)) + geom_point() 


ggplot(df2 %>% filter(clickCounter < 7), aes(clickCounter, chooseUncertainty)) + geom_smooth(method = "loess")
pdf("u_experience_by_sub.pdf", height = 30, width = 30)
ggplot(df2, aes(trial, as.numeric(u_present), color = u_rank)) + geom_line() + 
  facet_wrap(~id)
dev.off()
ggplot(df2, aes(ntrial_u_sampled)) + geom_bar()
ggplot(df2, aes(trial, outcome_sc)) + geom_smooth(method = "loess", se = F)
# after how many trials do they click on erasures?
pdf("u_counter_by_sub.pdf", height = 30, width = 30)
ggplot(df2, aes(trial, clickCounter, color = u_rank)) + geom_line() + 
  facet_wrap(~id)
dev.off()



m1 <- lmer(pos_shifted ~ scale(vmax_location)*vmax_scaled + (1|subject), df2)
summary(m1)
car::Anova(m1, '3')
# add the previous response/outcome (win-stay/lose-shift)
m1a_1 <- lmer(pos_shifted ~ scale(vmax_location)*scale(vmax)*experiment + resp_theta_c_lag * outcome_lag*experiment + (1|subject), df1 )
summary(m1a_1)
car::Anova(m1a_1, '3')
car::vif(m1a_1)
# add uncertainty
m2_1 <- lmerTest::lmer(pos_shifted ~ scale(vmax_location)*scale(vmax)*experiment + scale(u_location)*scale(vmax)*experiment + (1|subject), 
                       df1)
summary(m2_1)
car::Anova(m2_1, '3')
# full model with uncertainty and attentional control

m2_2 <- lmerTest::lmer(pos_shifted ~ scale(vmax_location)*scale(vmax) + scale(u_location)*scale(vmax) + 
                         resp_theta_c_lag*omission_lag + (1|subject), 
                       df1)
summary(m2_2)
car::Anova(m2_2, '3')

m2_3 <- lmerTest::lmer(pos_shifted ~ scale(vmax_location)*scale(vmax) + scale(att_location)*scale(vmax) + 
                         resp_theta_c_lag*omission_lag+ (1|subject), 
                       df1)
summary(m2_3)
car::Anova(m2_3, '3')

m5_1 <- lmerTest::lmer(pos_shifted ~ scale(vmax_location):pre_clock_freeze + scale(vmax_location):scale(vmax):pre_clock_freeze + 
                         scale(u_location):u_present:pre_clock_freeze + pre_clock_freeze*scale(att_location):att_present + 
                         resp_theta_c_lag*omission_lag*pre_clock_freeze + (1|subject), 
                       df)
summary(m5_1)
car::Anova(m5_1, '3')


# simulate various phantom erasure/control locations
# simpler alternative: use the last location
mlist <- list()
for (i in 1:1000) {
  df$u_location[!df$u_present] <- runif(length(df$u_location[!df$u_present]), 0, 360)
  df$att_location[!df$att_present] <- runif(length(df$att_location[!df$att_present]), 0, 360)
  mi <- lmer(pos_shifted ~ scale(vmax_location):epoch + 
               scale(u_location)*u_present + 
               scale(att_location)*att_present +
               resp_theta_c_lag*outcome_lag:epoch + outcome_lag +
               vmax_loc_lag1:epoch + vmax_loc_lag2:epoch +
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

m4 <- lmer(pos_shifted ~ scale(vmax_location)*scale(vmax) + 
             scale(u_location)*u_present + 
             scale(att_location) * att_present +
             resp_theta_c_lag * outcome_lag +
             (1|subject), 
           df)
summary(m4)
car::Anova(m4, '3')


# draft brms analyses with von Mises distribution

library(rstan)
rstan_options(auto_write = TRUE) 
options(mc.cores = parallel::detectCores()-1)
library(brms)
library(cmdstanr)
#set home and data directories
#SettingsØ

setwd(output_dir)

br.formula <- brmsformula(resp_theta_c ~ scale(vmax_theta_c) + scale(vmax_theta_c):scale(vmax) +
                            resp_theta_c_lag*omission_lag +
                            (1|subject), kappa ~ reward_lag) + von_mises()

br.formula <- brmsformula(resp_theta_c ~ scale(vmax_theta_c) + scale(vmax_theta_c):scale(vmax) + 
                            scale(u_theta_c):u_present + scale(u_theta_c):u_present:scale(vmax) + scale(att_theta_c):att_present + scale(att_theta_c):att_present:scale(vmax) + 
                            resp_theta_c_lag*omission_lag +
                            (1|subject), kappa ~ 1 + (1|subject)) + von_mises()
# noting that entering linear RT and sin/cos transformations of theta as predictors makes things only worse
# sincos.formula <- brmsformula(resp_theta_c ~ sin(vmax_theta_c) + cos(vmax_theta_c) + (1|subject), kappa ~ 1) + von_mises()
# *stan prior choice wiki*
brmsfamily = von_mises(link = "tan_half", link_kappa = "log");
prior_fit_von <- get_prior(br.formula, df1, family = brmsfamily,nl=TRUE)
# %>% filter(!class == "b" & !dpar == "kappa"),
# set_prior('normal(0, 1)', class='b', lb = -pi, ub = pi),
# # set_prior('normal(0, 1)', class='Intercept', dpar = "", lb = -pi, ub = pi),
# set_prior("gamma(2, .2)", class = "Intercept", dpar = "kappa"),
# set_prior('normal(0, .5)', class='b', dpar = "kappa")) 
# prior_fit_von$prior[prior_fit_von$class == "b"] <- "normal(0,1)"
# # with RT swings, kappa should be quite low, so should be the scale of gamma, but 0.1 seems too restrictive and 0.5 not enough
# prior_fit_von[prior_fit_von$class == "b"] <- set_prior('normal(0, 1)', class='b', lb = -pi, ub = pi)

# check Intercept prior!
# Prior_Theta = c(set_prior('normal(0, 1)', class='b', lb = -pi, ub = pi),
#                 set_prior("lkj(2)", class = "cor"),
#                 # set_prior("gamma(2, .2)", class = "Intercept", dpar = "kappa"),
#                 set_prior("normal(5, .8)", class = "Intercept", dpar = "kappa"),
#                 # set_prior("normal(0, 1) ", class = "intercept", dpar = "kappa", lb = 0),
#                 # set_prior("gamma(2, .2)", class = "b", dpar = "kappa"),
#                 set_prior('normal(0, .5)', class='b', dpar = "kappa"),
#                 # set_prior("normal(1,1)", class = "Intercept", dpar = "kappa"),
#                 set_prior("student_t(3, 0, 2.5)", class = "sd", lb = 0))
# # prior_fit_von$prior[prior_fit_von$dpar == "kappa"] <- "cauchy(0, 1)"
# bm6 <- brm(sincos.formula, data = df, prior = prior_fit_von, chains = 4, iter = 2000, backend = "cmdstanr", threads = threading(4))
bm_base <- brm(formula = br.formula, data = df1, prior = prior_fit_von, chains = 4, backend = "cmdstanr", iter=5000,cores=26,warmup = 1000,
               algorithm = "sampling")
summary(bm13)
bm13$prior

save(bm11, file = "clock2_brm_von_mises_m11.RData")
save(bm12, file = "clock2_brm_von_mises_m12.RData")

# bm2 -- good simple model, 

save(bm2, file = "clock2_brm_von_mises_m2.RData")

# simulate distribution from a prior

sim_formula <- brmsformula(resp_theta_c ~ 1 + 
                             #resp_theta_c_lag + resp_theta_c_lag:outcome_lag +   
                             (1|subject), kappa ~ 1) + von_mises()
prior_fit_von <- get_prior(sim_formula, df, family = "von_mises")

Prior_Theta = c(
  # set_prior("gamma(2, 0.01)", class = "Intercept", dpar = "kappa"),
  
  set_prior("student_t(3, 0, 10)", class = "sd"))
sm2= brm(sim_formula, data=df, prior = Prior_Theta,  sample_prior="only")
summary(sm2)
newY=predict(sm2)
str(newY)
hist(newY[,1])
hist(newY[,2])
hist(newY[,3])


# do we need to enter sin and cos for circular predictors? No.
sincos <- function(theta) {
  sine = sin(theta)
  cosine = cos(theta)
  print(c(sine, cosine), digits = 2)
}

# sanity check
ggplot(df, aes(trial, uncertainty_block, color = as.factor(chooseUncertainty))) + geom_jitter() + facet_wrap(~subject)
ggplot(df, aes(trial, uncertainty_block, color = as.factor(chooseFog1))) + geom_jitter() + facet_wrap(~subject)

ggplot(df, aes(meta_trialCount, uncertainty_block, color = as.factor(chooseAttentionalControl))) + geom_jitter() + facet_wrap(~subject)




# latency is the latency of the trial (for feedback always 1000, to check latency of choice/RT referenced to starting position filter by trialcode=='experiment_U' | trialcode=='experiment_noU')
# rt_shifted is the position around the clock chosen
# startPos is the starting position range [0 100] for one complete revolution, 12 o'clock is 75=100, 9 o'clock is 50, 6 o'clock is 25, 3 o'clock is 0
# meta_trialCount is trial [1... total number of trials ~ 300]
# trialCount is trial [1 ... trials within block].
# blockCount is block [1-8] blocks are 32-43 trials and contain 1 attentional control block and one uncertainty manipulation block
# ntrials is the number of trials per block (32-43)
# attentional_control is either "wind" or "cloud" this is the attentional control that subjects learn has no effect
# local_uncertainty is either "wind" or "cloud" this is the local uncertainty manipulation that subjects learn changes the number of mushrooms
# uncertaintyBlock is either 0,1,2 within each block, corresponding to the attentionalControl or localUncertainty manipulation
# chooseFog1 is 1 if "cloud" was chosen, something seems to be wrong with this counter 2022-05-03 AndyP, possibly not resetting to 0 after the block is complete?
# chooseWind2 is 1 if "wind" was chosen
# fogPos1 determines the position of the "cloud" stimulus, either [0,1,2,5,10,15,20,25,30,45,50] in 30 degree rotations counter clockwise starting at 12 o'clock for 0
# windPos1 determines the position of the "wind" stimulus, either [0,1,2,5,10,15,20,25,30,45,50] in 30 degree rotations counter clockwise starting at 12 o'clock for 0
# totalEarnings is cumulative earnings in units of 0.5c / mushroom, multiply by 2 to get # mushrooms
# Earnings is earnings on the current trial
# unc_att_start1 is the starting trial of the "cloud" stimulus on each block
# unc_att_start2 is the starting trial fo the "wind" stimulus on each block
# choose_unc_att1 counts the number of times the "cloud" stimulus is chosen (should not exceed 2)
# choose_unc_att2 counts the number of times the "wind" stimulus is chosen (should not exceed 2)
# inc is the number of points on that trial for the given rt_shifted
# rng is the random number generated to determine if the trial is an 'omission' (1 mushroom given).  We should have a fixed probability (freq) of 0.7 for all RT's and trials
# 