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

if (grepl("Nidhi", getwd())){
  thispath <- "C:/Users/Nidhi/OneDrive - University of North Carolina at Chapel Hill/Documents/GitHub/clock2/code"
  output_dir <- 'C:/Users/Nidhi/OneDrive - University of North Carolina at Chapel Hill/Documents/GitHub/clock2/Inquisit_design_files'
  source(file.path(thispath, "von_mises_basis.R"))
  source(file.path(thispath, "clock2_troll_world.R"))
  source(file.path(thispath, "scepticc_Nov_2023_commit.R"))
} else {
  output_dir <- '/Users/andypapale/clock2/Inquisit_design_files/'
  source('~/clock2/code/clock2_sim_crc_aep.R')
  source("~/clock2/code/von_mises_basis.R")
  source("~/clock2/code/clock2_troll_world.R")
  source("~/clock2/code/scepticc_Nov_2023_commit.R")
}
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
library(tidyverse)
inq_tri <- inq_tri %>% mutate(trial = row_number()) %>% rowwise() %>% pivot_longer(cols = starts_with("X"), names_to = "RT") %>% mutate(RT = readr::parse_number(RT))
inq_tri <- inq_tri %>% arrange(trial,RT)
if (!grepl("Nidhi", getwd())){
  setwd('/Users/andypapale/clock2/Inquisit_design_files/')
  write.csv(inq_tri,paste0('testing-Design-Matrix-with-Erasures-',as.character(df$iteration[i]),'.csv'))
} else{
  write.csv(inq_tri,paste0(output_dir, 'testing-Design-Matrix-with-Erasures-',as.character(df$iteration[i]),'.csv'))
}

# write erasure schedule
seg_min <- tt$erasure_segments$segment_min*180/pi;
seg_max <- tt$erasure_segments$segment_max*180/pi;
era_loc <- NULL;
era_val <- NULL;
val_min <- NULL;
val_max <- NULL;
for (iT in 1:300){
  if (!is.na(seg_max[iT]) && tt$erasure_segments$trial_type=="erasure"){
    if (seg_max[iT] < 31){
      # then era_loc >= 360
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
    temp <- inq_tri %>% filter(trial==iT);
    temp1 <- temp %>% filter(RT==round(era_loc[iT]));
    era_val[iT] <- temp1$value;
    temp2 <- temp %>% filter(RT==round(seg_min[iT]));
    val_min[iT] <- temp2$value;
    temp3 <- temp %>% filter(RT==round(seg_max[iT]));
    val_max[iT] <- temp3$value;
  } else {
    era_loc[iT] <- NA;
    era_val[iT] <- NA;
    val_min[iT] <- NA;
    val_max[iT] <- NA;
  }
}


## get values without erasures
set.seed(df$iteration[i])
ncenters <- 9 # how many gaussians there are
mean_val <- 10 # mean reward rate
sd_val <- 2 # standard deviation of reward / range of rewards
centers <- sample(seq(0, 2*pi, by = pi/20), ncenters, replace = FALSE) # line up gaussians here
values <- sample(truncnorm::rtruncnorm(ncenters, a = 0, mean = mean_val, sd = sd_val))
width_sd <- 0.349 # fixed, how wide are the underlying Gaussians
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
if (!grepl("Nidhi", getwd())){
  setwd('/Users/andypapale/clock2/Inquisit_design_files/')
  write.csv(inq_tri,paste0('testing-Design-Matrix-',as.character(df$iteration[i]),'.csv'))
  d9 <- read_csv(paste0('/Users/andypapale/clock2/Inquisit_design_files/testing-Design-Matrix-',as.character(df$iteration[i]),'.csv'))
} else {
  write.csv(inq_tri,paste0(output_dir,'/testing-Design-Matrix-',as.character(df$iteration[i]),'.csv'))
  d9 <- read_csv(paste0(output_dir,'/testing-Design-Matrix-',as.character(df$iteration[i]),'.csv'))
}

############################################
# Create splines from Inquisit data
############################################
df_wo_era <- d9 # read Inquisit design matrix without erasures
pacman::p_load(officer, ggplot2)

df_splines <- data.frame(trial = as.numeric(),
                         RTs = as.numeric(),
                         value = as.numeric(),
                         era_loc = as.numeric(),
                         era_val = as.numeric())
# ppt <- read_pptx()  # Create a new PowerPoint document
for (t in 1:300){
  if (!is.na(era_loc[t])){
    v <- df_wo_era$value[df_wo_era$trial==t] # 360 values for this trial before removal of erasure locations
    v_w_era <- v
    if (seg_max[t] < 31){ # erasure segment contains the end points
      erasure_segment <- c(seg_min[t]:360, 1:seg_max[t])
      v_w_era[erasure_segment] <- NA
      v_w_era[era_loc[t]] <- era_val[t]
      v_temp <- v_w_era
      v_temp <- c(v_temp[181:360],v_temp[1:180]) # bringing the erasure to the middle so that the spline works
      v_temp <- zoo::na.spline(v_temp)
      v_w_era <- c(v_temp[181:360],v_temp[1:180]) # bringing the erasure back to the original location
    } else {
      erasure_segment <- seg_min[t]:seg_max[t]
      v_w_era[erasure_segment] <- NA
      v_w_era[era_loc[t]] <- era_val[t]
      v_w_era <- zoo::na.spline(v_w_era)
    }
    
    df_temp <- data.frame(trial = t, RTs = 1:360, v=v, v_w_era=v_w_era)
    df_splines <- rbind(df_splines, data.frame(trial = t, RTs = 1:360, value = v_w_era, 
                                               era_loc = era_loc[t], era_val = era_val[t]))
    print(t)
    # gg1 <- ggplot(df_temp, aes(x=RTs, y=v_w_era)) +
      # geom_line() +
      # geom_point(size=1) +
      # geom_line(aes(y=v_w_era), color='red', size=1) + 
      # points(era_loc[t], era_val[t], col='blue', pch=19) +
      # ylim(0,max(df_wo_era$value)) +
      # ggtitle(paste0('Trial ',t, ' with erasure at ', era_loc[t], ' and value ', era_val[t]))
    # slide <- add_slide(ppt, layout='Title and Content')
    # slide <- ph_with(slide, value = gg1, location = ph_location_fullsize())
    # print(ppt, paste0("v_with_era_testing_seed-",df$iteration[1],".pptx"))
  }
}

############################################
# Compare splines in R to Inquisit data
############################################

# Pull in Inquisit output files
if (grepl("Nidhi", getwd())){
  df_inq <- read_csv("C:/Users/Nidhi/OneDrive - University of North Carolina at Chapel Hill/Documents/GitHub/Clock_v2/data/clock_2_1_1_seed_152_testing_raw_1_2024-01-23-21-25-55-768.csv")
} else {
  df_inq <- read_csv('/Users/andypapale/Desktop/2024-01-Testing/data/clock_2_1_1_seed_152_testing_raw_1_2024-01-24-18-03-24-721.csv')
}
df_inq <- df_inq %>% filter(trialcode=="getValueVector")
ppt <- read_pptx()  # Create a new PowerPoint document
for (t in 1:300){
  if (!is.na(era_loc[t])){
    print(t)
    gg1 <- ggplot(df_splines %>% filter(trial == t), aes(x=RTs, y=value)) +
      geom_point(col='black', pch=19, size=2) +
      geom_point(aes(x=era_loc[1],y=era_val[1]),color='blue',shape=2,size=8) +
      geom_point(data=df_inq %>% filter(trial == t), aes(x=x_final, y=y_final), color = 'red', size = 0.8) +
      ylim(0,max(df_wo_era$value)) +
      ggtitle(paste0('Trial ',t, ' with erasure at ', era_loc[t], ' and value ', era_val[t]))
    slide <- add_slide(ppt, layout='Title and Content')
    slide <- ph_with(slide, value = gg1, location = ph_location_fullsize())
    print(ppt, paste0("inq-output-red_inq-input-black_testing_seed-",df$iteration[1],".pptx"))
  }
}



