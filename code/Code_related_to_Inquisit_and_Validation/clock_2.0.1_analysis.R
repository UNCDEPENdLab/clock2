library(tidyverse)
library(reshape2)
library(lme4)
library(brms)
# 2023-05-03 AndyP
# clock 2.0.1 analysis


# set paths (will be good to make this more user-general)
sys <- Sys.info()
if (sum(str_detect(Sys.info(), "Alex"))>1) { base_dir <- "~/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Documents/"} else {
  base_dir <- '~/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Documents - DNPLskinner' }# on skinner}
output_dir <- "~/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Documents/skinner/data/prolific/clock_v2_pilot"
# design_file <- '~/code/clock2/2022-04-25-DesignFile.csv' # in Michael's clock2 repo
design_file <- '/Volumes/Users/Andrew/2022-05-08-DesignFile.csv'
design <- as.matrix(read_csv(design_file)) %>% 
  as_tibble() %>% select(-`...1`) %>% mutate(timepoint = row_number()) %>% rowwise() %>% pivot_longer(cols = starts_with("V"), names_to = "trial") %>%
  mutate(trial = extract_numeric(trial)) %>% group_by(trial) %>% summarise(vmax = max(value),
                                                                                   vmax_location = timepoint[which.max(value)])
# initial sanity checks on the disign file
str(design)
ggplot(design, aes(trial, vmax_location, color = vmax)) + geom_line() + scale_color_viridis_c()
ggplot(design, aes(trial, vmax)) + geom_line() + scale_color_viridis_c()
median_vmax = median(design$vmax)

# combine behavioral data from two batches of 25 Prolific subjects
df1 <- read_csv('/Volumes/Users/Andrew/papalea_prosper_eeg_clock_v2_0_3_raw_2307031911.csv') %>% filter(trialcode=='dispFeedback_noU' | trialcode=='dispFeedback_U') # the score is calculated here, the variable of interest is rt_shifted
df1 <- df1 %>%  mutate(u_present = (!is.na(windPos_out) & local_uncertainty == "wind") | (!is.na(fogPos_out) & local_uncertainty == "cloud"), 
         att_present = (!is.na(windPos_out) & local_uncertainty == "cloud") | (!is.na(fogPos_out) & local_uncertainty == "wind"),
         u_location = case_when(
           !is.na(windPos_out) & local_uncertainty == "wind" ~ windPos_out,
           !is.na(fogPos_out) & local_uncertainty == "fog" ~ fogPos_out,
           T ~ runif(1, 0, 360)
         ),
         att_location = case_when(
           !is.na(windPos_out) & local_uncertainty == "fog" ~ windPos_out,
           !is.na(fogPos_out) & local_uncertainty == "wind" ~ fogPos_out,
           T ~ runif(1, 0, 360))
         ) %>% rename(trial = meta_trialCount)
  
  df1 <- df1 %>% mutate(trial_to_merge = list.rt_0.currentindex)
  design <- design %>% mutate(trial_to_merge = trial)
  
  df1 <- inner_join(df1,design,by='trial')
  df1 <- df1 %>%mutate(resp_theta = pos_shifted * pi/180,
         resp_theta_c = resp_theta - pi,
         vmax_theta = vmax_location * pi/180,
         vmax_theta_c = vmax_theta - pi,
         u_theta = u_location * pi/180,
         u_theta_c = u_location - pi,
         att_theta = att_location * pi/180,
         att_theta_c = att_theta - pi) %>% group_by(subject) %>% arrange(trial) %>%
  mutate(resp_theta_c_lag = lag(resp_theta_c),
         vmax_theta_c_lag = lag(vmax_theta_c),
         outcome_lag = scale(lag(inc)),
         omission = inc == 1,
         omission_lag = lag(omission),
         reward_lag = !omission_lag)

# sanity checks on behavioral data: ensure we know where erasures/control stimuli were
ggplot(df1, aes(trial, u_present, color = as.factor(chooseUncertainty))) + geom_point() + facet_wrap(~subject)
pdf('test.pdf',height=18,width=16)
gg1 <- ggplot(df1, aes(trial, local_uncertainty, color = as.factor(chooseFog1))) + geom_point(size=1) + facet_wrap(~subject) + geom_text(aes(x = 50, y = .5, label = local_uncertainty))
print(gg1)
dev.off()
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


# initial Gaussian glms: good for manipulation checks
# only the value bump
m1 <- lmer(pos_shifted ~ scale(vmax_location)*scale(vmax) + (1|subject), df1 )
summary(m1)
car::Anova(m1, '3')
# add the previous response/outcome (win-stay/lose-shift)
m1a_1 <- lmer(pos_shifted ~ scale(vmax_location)*scale(vmax) + resp_theta_c_lag * outcome_lag + (1|subject), df1 )
summary(m1a_1)
car::Anova(m1a_1, '3')
car::vif(m1a_1)
# add uncertainty
m2_1 <- lmerTest::lmer(pos_shifted ~ scale(vmax_location)*scale(vmax) + scale(u_location)*scale(vmax) + (1|subject), 
           df1)
summary(m2_1)
car::Anova(m2_1, '3')
# full model with uncertainty and attentional control
m5_1 <- lmerTest::lmer(pos_shifted ~ scale(vmax_location) + scale(vmax_location):scale(vmax) + 
               scale(u_location):u_present + scale(u_location):u_present:scale(vmax) + scale(att_location):att_present + scale(att_location):att_present:scale(vmax) + 
             resp_theta_c_lag*omission_lag +
             (1|subject), 
           df1)
summary(m5_1)
car::Anova(m5_1, '3')

# simulate various phantom erasure/control locations
# simpler alternative: use the last location
mlist <- list()
for (i in 1:1000) {
  df$u_location[!df$u_present] <- runif(length(df$u_location[!df$u_present]), 0, 360)
  df$att_location[!df$att_present] <- runif(length(df$att_location[!df$att_present]), 0, 360)
  mi <- lmer(pos_shifted ~ scale(vmax_location)*scale(vmax) + 
               scale(u_location)*u_present*scale(vmax) + 
               scale(att_location) * att_present*scale(vmax) +
               resp_theta_c_lag * outcome_lag +
               (1|subject), 
             df1)
  mdf <- broom.mixed::tidy(mi)
  mdf$i <- i
  mlist[[i]] <- mdf
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

br.formula <- brmsformula(resp_theta_c ~ vmax_theta_c + vmax_theta_c + 
                            resp_theta_c_lag +
                            # resp_theta_c_lag + resp_theta_c_lag:reward_lag +
                            (1|subject), kappa ~ reward_lag) + von_mises()

# noting that entering linear RT and sin/cos transformations of theta as predictors makes things only worse
# sincos.formula <- brmsformula(resp_theta_c ~ sin(vmax_theta_c) + cos(vmax_theta_c) + (1|subject), kappa ~ 1) + von_mises()
# *stan prior choice wiki*
prior_fit_von <- c(get_prior(br.formula, df1, family = "von_mises")) 
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
bm12 <- update(bm12, formula = br.formula, newdata = df1, prior = prior_fit_von, chains = 4, iter = 8000, backend = "cmdstanr", threads = threading(5),
               control = list(adapt_delta = 0.8))
summary(bm12)
bm12$prior
bm13 <- brm(formula = br.formula, data = df1, prior = prior_fit_von, chains = 4, backend = "cmdstanr", 
               algorithm = "meanfield")
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