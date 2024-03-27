# 2024-03-27 AndyP
# Spline in Troll World using extracted value vector
# without erasures, and pulling out points from erasures 
# to then reimplement spline one erasure at a time the way it is done in Inquisit
plot = T

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