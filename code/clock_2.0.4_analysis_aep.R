library(tidyverse)
library(reshape2)
library(lme4)
#library(brms)
library(data.table)
# 2023-05-03 AndyP
# clock 2.0.1 analysis


design_file <- '/Users/andypapale/clock2/2023-09-26-Design-File-asMatrix.csv'
design <- as.matrix(read_csv(design_file)) %>% 
  as_tibble() %>% select(-`...1`) %>% mutate(timepoint = row_number()) %>% rowwise() %>% pivot_longer(cols = starts_with("V"), names_to = "trial") %>%
  mutate(trial = extract_numeric(trial)) %>% group_by(trial) %>% summarise(vmax = max(value),
                                                                                   vmax_location = timepoint[which.max(value)])
design <- design %>% mutate(epoch = case_when(vmax > 120 ~ 'High Value', vmax <= 120 ~ 'Low Value'))
design <- design %>% mutate(epoch_bin = as.factor(case_when(epoch=='High Value' & trial < 36 ~ 1,
                                                  epoch=='High Value' & trial > 36 & trial < 80 ~ 2,
                                                  epoch=='High Value' & trial > 80 & trial < 116 ~ 3,
                                                  epoch=='High Value' & trial > 116 & trial < 154 ~ 4,
                                                  epoch=='High Value' & trial > 154 & trial < 195 ~ 5,
                                                  epoch=='High Value' & trial > 195 & trial < 233 ~ 6,
                                                  epoch=='High Value' & trial > 233 & trial < 270 ~ 7,
                                                  epoch=='High Value' & trial > 270 & trial < 306 ~ 8,
                                                  epoch=='High Value' & trial > 306 ~ 9,
                                                  epoch=='Low Value' & trial < 84 ~ -1,
                                                  epoch=='Low Value' & trial > 84 & trial < 120 ~ -2,
                                                  epoch=='Low Value' & trial > 120 & trial < 158 ~ -3,
                                                  epoch=='Low Value' & trial > 158 & trial < 199 ~ -4,
                                                  epoch=='Low Value' & trial > 199 & trial < 237 ~ -5,
                                                  epoch=='Low Value' & trial > 237 & trial < 274 ~ -6,
                                                  epoch=='Low Value' & trial > 274 & trial < 310 ~ -7,
                                                  epoch=='Low Value' & trial > 310 ~ -8)))
design <- design %>% group_by(epoch_bin) %>% mutate(vmax_loc_mean = mean(vmax_location), vmax_mean = mean(vmax)) %>% ungroup()
design <- design %>% mutate(vmax_loc_lag1 = case_when(
                                                      epoch_bin==2 ~ mean(design$vmax_loc_mean[design$epoch_bin==1],na.rm=TRUE),
                                                      epoch_bin==3 ~ mean(design$vmax_loc_mean[design$epoch_bin==2],na.rm=TRUE),
                                                      epoch_bin==4 ~ mean(design$vmax_loc_mean[design$epoch_bin==3],na.rm=TRUE),
                                                      epoch_bin==5 ~ mean(design$vmax_loc_mean[design$epoch_bin==4],na.rm=TRUE),
                                                      epoch_bin==6 ~ mean(design$vmax_loc_mean[design$epoch_bin==5],na.rm=TRUE),
                                                      epoch_bin==7 ~ mean(design$vmax_loc_mean[design$epoch_bin==6],na.rm=TRUE),
                                                      epoch_bin==8 ~ mean(design$vmax_loc_mean[design$epoch_bin==7],na.rm=TRUE),
                                                      epoch_bin==9 ~ mean(design$vmax_loc_mean[design$epoch_bin==8],na.rm=TRUE),
                                                      epoch_bin==-1 ~ mean(design$vmax_loc_mean[design$epoch_bin==1],na.rm=TRUE),
                                                      epoch_bin==-2 ~ mean(design$vmax_loc_mean[design$epoch_bin==2],na.rm=TRUE),
                                                      epoch_bin==-3 ~ mean(design$vmax_loc_mean[design$epoch_bin==3],na.rm=TRUE),
                                                      epoch_bin==-4 ~ mean(design$vmax_loc_mean[design$epoch_bin==4],na.rm=TRUE),
                                                      epoch_bin==-5 ~ mean(design$vmax_loc_mean[design$epoch_bin==5],na.rm=TRUE),
                                                      epoch_bin==-6 ~ mean(design$vmax_loc_mean[design$epoch_bin==6],na.rm=TRUE),
                                                      epoch_bin==-7 ~ mean(design$vmax_loc_mean[design$epoch_bin==7],na.rm=TRUE),
                                                      epoch_bin==-8 ~ mean(design$vmax_loc_mean[design$epoch_bin==8],na.rm=TRUE))
                            ,
                            vmax_loc_lag2 = case_when(
                              epoch_bin==3 ~ mean(design$vmax_loc_mean[design$epoch_bin==1],na.rm=TRUE),
                              epoch_bin==4 ~ mean(design$vmax_loc_mean[design$epoch_bin==2],na.rm=TRUE),
                              epoch_bin==5 ~ mean(design$vmax_loc_mean[design$epoch_bin==3],na.rm=TRUE),
                              epoch_bin==6 ~ mean(design$vmax_loc_mean[design$epoch_bin==4],na.rm=TRUE),
                              epoch_bin==7 ~ mean(design$vmax_loc_mean[design$epoch_bin==5],na.rm=TRUE),
                              epoch_bin==8 ~ mean(design$vmax_loc_mean[design$epoch_bin==6],na.rm=TRUE),
                              epoch_bin==9 ~ mean(design$vmax_loc_mean[design$epoch_bin==7],na.rm=TRUE),
                              epoch_bin==-2 ~ mean(design$vmax_loc_mean[design$epoch_bin==1],na.rm=TRUE),
                              epoch_bin==-3 ~ mean(design$vmax_loc_mean[design$epoch_bin==2],na.rm=TRUE),
                              epoch_bin==-4 ~ mean(design$vmax_loc_mean[design$epoch_bin==3],na.rm=TRUE),
                              epoch_bin==-5 ~ mean(design$vmax_loc_mean[design$epoch_bin==4],na.rm=TRUE),
                              epoch_bin==-6 ~ mean(design$vmax_loc_mean[design$epoch_bin==5],na.rm=TRUE),
                              epoch_bin==-7 ~ mean(design$vmax_loc_mean[design$epoch_bin==6],na.rm=TRUE),
                              epoch_bin==-8 ~ mean(design$vmax_loc_mean[design$epoch_bin==7],na.rm=TRUE))
)
# initial sanity checks on the disign file
str(design)
ggplot(design, aes(trial, vmax_location, color = vmax)) + geom_line() + scale_color_viridis_c()
ggplot(design, aes(trial, vmax)) + geom_line() + scale_color_viridis_c()
#median_vmax = median(design$vmax)

# combine behavioral data from two batches of 25 Prolific subjects
df1 <- read_csv('/Users/andypapale/Downloads/clock_raw_data_all_9-24.csv') %>%
  filter(trialcode=='dispFeedback_noU' | trialcode=='dispFeedback_U') # the score is calculated here, the variable of interest is rt_shifted
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
  
  df1 <- inner_join(df1,design,by='trial_to_merge')
  df1 <- df1 %>% mutate(resp_theta = pos_shifted * pi/180,
         resp_theta_c = resp_theta - pi,
         vmax_theta = vmax_location * pi/180,
         vmax_theta_c = vmax_theta - pi,
         u_theta = u_location * pi/180,
         u_theta_c = u_theta - pi,
         att_theta = att_location * pi/180,
         att_theta_c = att_theta - pi) %>% group_by(subject) %>% arrange(trial_to_merge) %>%
  mutate(resp_theta_c_lag = lag(resp_theta_c),
         vmax_theta_c_lag = lag(vmax_theta_c),
         vmax_loc_scaled = scale(vmax_location),
         outcome_lag = lag(Earnings),
         omission = 0.7 <= rng,
         omission_lag = lag(omission),
         reward_lag = !omission_lag) %>% ungroup()
  df1$pre_clock_freeze = df1$list.preClockFreeze.currentvalue
  
  df1 <- df1 %>% arrange(subject,blocknum,meta_trialCount)

df2 <- df1 %>% filter(epoch=='High Value')
ggplot(df2, aes(x=pos_shifted,y=inc_rg)) + facet_wrap(~epoch_bin) + geom_line() + geom_point(aes(color=omission)) + geom_point(aes(x=vmax_location,y= 150), shape = 2,color='green')
    
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
df2 <- df2 %>% select(pos_shifted,vmax_loc_scaled,vmax_scaled,resp_theta_c_lag,omission_lag,subject)
df2 <- df2 %>% arrange(subject)
df2 <- transform(df2, id=match(subject, unique(subject)))
df2$id <- as.numeric(df2$id)
df2$pos_shifted <- df2$pos_shifted * pi/180
library(bpnreg)
m5 <- bpnreg::bpnme(pred.I = pos_shifted ~ vmax_loc_scaled*vmax_scaled + resp_theta_c_lag*omission_lag + (1|id), data = df2, its = 200, burn = 100, n.lag=3, seed=121)


# initial Gaussian glms: good for manipulation checks
# only the value bump
m1 <- lmerTest::lmer(pos_shifted ~ scale(vmax_location)*scale(vmax) + (1|subject), df1 )
summary(m1)
car::Anova(m1, '3')
# add the previous response/outcome (win-stay/lose-shift)
m1a_1 <- lmerTest::lmer(pos_shifted ~ scale(vmax_location)*scale(vmax)*experiment + resp_theta_c_lag * outcome_lag*experiment + (1|subject), df1 )
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
  mi <- lmerTest::lmer(pos_shifted ~ scale(vmax_location)*scale(vmax)*gauss_sd + 
               scale(u_location)*u_present*gauss_sd + 
               scale(att_location)*att_present*gauss_sd +
               resp_theta_c_lag*outcome_lag*gauss_sd +
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
set_cmdstan_path(path='~/.cmdstan/cmdstan-2.33.1')
setwd(getwd())

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
# ') %>% filter(trialcode=='dispFeedback_noU' | trialcode=='dispFeedback_U') # the score is calculated here, the variable of interest is rt_shifted
df1 <- df1 %>% mutate(experiment = 'trolls')
df1 <- rbind(df1, dfm)
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
                      att_theta_c = att_theta - pi) %>% group_by(subject) %>% arrange(trial_to_merge) %>%
  mutate(resp_theta_c_lag = lag(resp_theta_c),
         vmax_theta_c_lag = lag(vmax_theta_c),
         vmax_loc_scaled = scale(vmax_location),
         outcome_lag = lag(Earnings),
         omission = 0.7 >= rng,
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
df2 <- df2 %>% select(experiment,pos_shifted,vmax_loc_scaled,vmax_scaled,u_loc_scaled,u_present,att_present,att_loc_scaled,resp_theta_c_lag,omission_lag,subject)
df2 <- df2 %>% arrange(subject)
df2 <- transform(df2, id=match(subject, unique(subject)))
df2$id <- as.numeric(df2$id)
df2$pos_shifted <- df2$pos_shifted * pi/180
library(bpnreg)
m5 <- bpnreg::bpnme(pred.I = pos_shifted ~ vmax_loc_scaled*vmax_scaled*experiment + u_loc_scaled*u_present*vmax_scaled + att_loc_scaled*att_present*vmax_scaled + resp_theta_c_lag*omission_lag + (1|id), data = df2, its = 2000, burn = 100, n.lag=3, seed=121)


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








