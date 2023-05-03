library(tidyverse)
library(reshape2)
library(lme4)
# 2023-05-03 AndyP
# clock 2.0.1 analysis

sys <- Sys.info()
if (sum(str_detect(Sys.info(), "Alex"))>1) { base_dir <- "~/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Documents/"} else {
  base_dir <- '~/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Documents - DNPLskinner' }# on skinner}

design_file <- '~/code/clock2/2022-04-25-DesignFile.csv' # in Michael's clock2 repo

design <- as.matrix(read_csv(design_file)) %>% 
  as_tibble() %>% select(-`...1`) %>% mutate(trial = row_number()) %>% rowwise() %>% pivot_longer(cols = starts_with("V"), names_to = "timepoint") %>%
  mutate(timepoint = extract_numeric(timepoint)) %>% group_by(trial) %>% summarise(vmax = max(value),
                                                                                   vmax_location = timepoint[which.max(value)])
str(design)
ggplot(design, aes(trial, loc_vmax, color = vmax)) + geom_line() + scale_color_viridis_c()

ggplot(design, aes(timepoint, value, color = trial)) + geom_line()


df <- read_csv(paste0(base_dir,'/skinner/data/prolific/clock_v2.1_pilot/In_lab_pilot_05-2023/papalea_prosper_eeg_clock_v3_raw_2305011503.csv')) %>% 
  filter(trialcode=='dispFeedback_noU' | trialcode=='dispFeedback_U') %>%# the score is calculated here, the variable of interest is rt_shifted
  mutate(chooseFog1 = chooseFog1*(uncertainty_block>0),
         chooseWind2 = chooseWind2*(uncertainty_block>0),
         chooseUncertainty = chooseUncertainty*(uncertainty_block>0),
         chooseAttentionalControl = chooseAttentionalControl*(uncertainty_block>0),
         chosen_location = round(rt_shifted*72/1000),
         uncertainty_location = case_when(
           attentional_control == "wind" ~ dense_rank(fogPos1)*30 - 30,
           attentional_control == "cloud" ~ dense_rank(windPos2)*30 - 30),
         control_location = case_when(
           attentional_control == "cloud" ~ dense_rank(fogPos1)*30 - 30,
           attentional_control == "wind" ~ dense_rank(windPos2)*30 - 30),
         uncertainty_present = case_when(
           local_uncertainty == "cloud" ~ uncertainty_block == 1,
           local_uncertainty == "wind" ~ uncertainty_block == 2),
         control_present = case_when(
           local_uncertainty == "cloud" ~ uncertainty_block == 2,
           local_uncertainty == "wind" ~ uncertainty_block == 1)
         
  ) %>% rename(trial = meta_trialCount) %>% inner_join(design, by = "trial")

ggplot(df, aes(trial, uncertainty_block, color = as.factor(chooseUncertainty))) + geom_point() + facet_wrap(~subject)
ggplot(df, aes(trial, uncertainty_block, color = as.factor(chooseFog1))) + geom_point() + facet_wrap(~subject) + geom_text(aes(x = 50, y = .5, label = local_uncertainty))
ggplot(df, aes(trial, uncertainty_block, color = as.factor(chooseWind2))) + geom_point() + facet_wrap(~subject) + geom_text(aes(x = 50, y = .5, label = local_uncertainty))
ggplot(df, aes(trial, uncertainty_block, color = as.factor(uncertainty_present))) + geom_point() + facet_wrap(~subject) + geom_text(aes(x = 50, y = .5, label = local_uncertainty))

ggplot(df, aes(trial, uncertainty_block, color = as.factor(chooseAttentionalControl))) + geom_point() + facet_wrap(~subject)

ggplot(df, aes(trial, uncertainty_block, color = as.factor(chooseUncertainty))) + geom_jitter()


ggplot(df, aes(vmax_location, chosen_location, alpha = vmax, color = subject)) + geom_point() + geom_smooth(method = "loess")
ggplot(df, aes(uncertainty_location, chosen_location,  color = subject)) + geom_point() + geom_smooth(method = "loess")

# need to switch to bpnreg::fit.bpnme; it does not handle missing values
# library(bpnreg)


m1 <- lmer(chosen_location ~ scale(vmax_location)*scale(vmax) + (1|subject), df )
summary(m1)
car::Anova(m1, '3')

m2 <- lmer(chosen_location ~ scale(vmax_location)*scale(vmax) + scale(uncertainty_location)*uncertainty_present + scale(control_location)*control_present + (1|subject), 
           df)
summary(m2)
car::Anova(m2, '3')


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