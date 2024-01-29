# 20204-01-06 AndyP
# Testing Seeds
generate_inquisit_lists = T
plot = F

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
i = 1
j = 1
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


if (plot){
  if (!file.exists(paste0('/Users/andypapale/Desktop/seed-',df$seed[j],'-iteration-',as.character(df$iteration[i]),'-validation2'))){
    dir.create(paste0('/Users/andypapale/Desktop/seed-',df$seed[j],'-iteration-',as.character(df$iteration[i]),'-validation2'))
  }
  
  ## compare original value distribution to value distribution with erasures from the same seed
  setwd(paste0('/Users/andypapale/Desktop/seed-',df$seed[j],'-iteration-',as.character(df$iteration[i]),'-validation2'))
  val_cat <- NULL;
  oval_cat <- NULL;
  loc_cat <- NULL;
  e_trial <- NULL;
  for (iT in 1:300){
    e_val <- d9_w_era %>% filter(trial==iT);
    no_eval <- d9 %>% filter(trial==iT);
    no_eval <- no_eval$value[round(e_val$RT)==round(era_loc[iT])]
    e_val <- e_val$value[round(e_val$RT)==round(era_loc[iT])];
    df_era <- data.frame(e_pos = era_loc[iT],e_val = e_val)
    df_noera <- data.frame(e_pos = era_loc[iT],no_eval = no_eval)
    pdf(paste0('iteration-',as.character(df$iteration[i]),'-trial-',iT,'.pdf'),height=12,width=12);
    if (tt$erasure_segments$trial_type[iT] == 'erasure'){
      gg1 <- ggplot(d9 %>% filter(trial==iT), aes(x=RT,y=value), color='black') +
        geom_line() + geom_point(size=5,color='black') +
        geom_line(data = d9_w_era %>% filter(trial==iT),aes(x=RT,y=value),color='red') +
        geom_point(data = d9_w_era %>% filter(trial==iT),aes(x=RT,y=value),color='red') +
        ylim(0,200) + geom_point(data = df_era, aes(x=e_pos,y=e_val),color='black',shape=2,size=10) +
        ggtitle(paste0('trial ',iT,' trial_type = ', tt$erasure_segments$trial_type[iT],' erasure_val = ',e_val, ' erasure_RT = ',era_loc[iT]));
    } else {
      gg1 <-  ggplot(d9 %>% filter(trial==iT), aes(x=RT,y=value), color='black') +
        geom_line() + geom_point(size=5,color='black') +
        geom_line(data = d9_w_era %>% filter(trial==iT),aes(x=RT,y=value),color='red') +
        geom_point(data = d9_w_era %>% filter(trial==iT),aes(x=RT,y=value),color='red') +
        ylim(0,200) +
        ggtitle(paste0('trial ',iT,' trial_type = ', tt$erasure_segments$trial_type[iT]));
    }
    # gg1 <- ggplot(d9 %>% filter(trial==iT), aes(x=RT,y=value)) + geom_line() + geom_point(size=5,color='black') +  ylim(0,200) +
    #   ggtitle(paste0('trial ',iT,' trial_type = ', bb$erasure_segments$trial_type[iT], ' e_val = ',e_val, ' e_max = ',round(tt$erasure_segments$segment_max[iT]*180/pi)))
    print(gg1);
    dev.off()
  }
}
if (generate_inquisit_lists){
  
  setwd('/Users/andypapale/clock2/Inquisit_design_files')
  
  aa <- data.frame(aa)
  aa <- aa %>% mutate(trial = row_number()) %>% rowwise() %>% pivot_longer(cols = starts_with("X"), names_to = "RT") %>% mutate(RT = extract_numeric(RT))
  aa <- aa %>% arrange(trial,RT)
  write.csv(aa,paste0('Design-Matrix-',as.character(df$iteration),'.csv'))
  epoch <- data.frame(tt$epoch)
  write.csv(epoch,paste0('epoch-',as.character(df$iteration),'.csv'))
  
  #generate value, RT and trial lists as 1 x (nT x nRT) inquisit lists
  options("encoding" = "UTF-8") # encode in UTF-8 as suggested here https://forums.millisecond.com/Topic15777.aspx#15778
  df0 <- NULL;
  dq0 <- NULL;
  dz0 <- NULL;
  nR <- nrow(aa);
  for (iR in 1:nR){
    if (iR==1){
      df0 <- paste0('<list values>\n/ items = (',as.character(aa$value[iR]),',');
      dq0 <- paste0('<list RT>\n/ items = (',as.character(aa$RT[iR]),',');
      dz0 <- paste0('<list trial>\n/ items = (',as.character(aa$trial[iR]),',');
    } else if (iR > 1 && iR < nR){
      df0 <- paste0(df0,as.character(aa$value[iR]),',');
      dq0 <- paste0(dq0,as.character(aa$RT[iR]),',');
      dz0 <- paste0(dz0,as.character(aa$trial[iR]),',');
    } else if (iR==nR){
      df0 <- paste0(df0,as.character(aa$value[iR]),')\n/ selectionrate = always\n/ selectionmode = values.master_idx;\n</list>')
      dq0 <- paste0(dq0,as.character(aa$RT[iR]),')\n/ selectionrate = always\n/ selectionmode = values.master_idx;\n</list>')
      dz0 <- paste0(dz0,as.character(aa$trial[iR]),')\n/ selectionrate = always\n/ selectionmode = values.master_idx;\n</list>')
    }
    if ((iR %% 1000)==0){
      print(iR/nR);
    }
  }
  write.table(df0,paste0('values-',as.character(df$iteration),'.txt'),row.names=F,col.names=F,quote=F)
  write.table(dq0,paste0('RTs-',as.character(df$iteration),'.txt'),row.names=F,col.names=F,quote=F)
  write.table(dz0,paste0('trials-',as.character(df$iteration),'.txt'),row.names=F,col.names=F,quote=F)
  options("encoding" = "native.enc") # change encoding back to native
  
  trial_type <- tt$erasure_segments$trial_type
  drift <- tt$get_drift_vec()
  spread <- tt$spread
  options("encoding" = "UTF-8")
  df0 <- NULL;
  dq0 <- NULL;
  dz0 <- NULL;
  dh0 <- NULL;
  nR <- length(era_loc);
  for (iR in 1:nR){
    
    #if (is.na(era_loc[iR])){
    #  temp <- 'NULL';
    #} else {
    #  temp <- era_loc[iR];
    #}
    
    if (iR==1){
      #df0 <- paste0('<list erasure_locations>\n/ items = (',as.character(temp),',');
      dq0 <- paste0('<list block_type>\n/ items = ("',as.character(trial_type[iR]),'",');
      dz0 <- paste0('<list drift>\n/ items = (',as.character(drift[iR]),',');
      dh0 <- paste0('<list spread>\n/ items = (',as.character(spread[iR]),',');
    } else if (iR > 1 && iR < nR){
      #df0 <- paste0(df0,as.character(temp),',');
      dq0 <- paste0(dq0,'"',as.character(trial_type[iR]),'",');
      dz0 <- paste0(dz0,as.character(drift[iR]),',');
      dh0 <- paste0(dh0,as.character(spread[iR]),',');
    } else if (iR==nR){
      #df0 <- paste0(df0,as.character(temp),')\n/ selectionrate = always\n/ selectionmode = values.era_loc_index;\n</list>')
      dq0 <- paste0(dq0,'"',as.character(trial_type[iR]),'")\n/ selectionrate = always\n/ selectionmode = values.trial;\n</list>')
      dz0 <- paste0(dz0,as.character(drift[iR]),')\n/ selectionrate = always\n/ selectionmode = values.trial; \n</list>')
      dh0 <- paste0(dh0,as.character(spread[iR]),')\n/ selectionrate = always\n/ selectionmode = values.trial; \n</list>')
    }
  }
  #write.table(df0,paste0('era_loc-',as.character(df$iteration),'.txt'),row.names=F,col.names=F,quote=F)
  write.table(dq0,paste0('block_type-',as.character(df$iteration),'.txt'),row.names=F,col.names=F,quote=F)
  write.table(dz0,paste0('drift-vector-',as.character(df$iteration),'.txt'),row.names=F,col.names=F,quote=F)
  write.table(dh0,paste0('spread-',as.character(df$iteration),'.txt'),row.names=F,col.names=F,quote=F)
  options("encoding" = "native.enc") # change encoding back to native
  
  era_val = NULL
  era_loc1 = NULL
  att_loc = NULL
  iC = 1;
  iD = 1;
  for (iT in 1:300){
    if (tt$erasure_segments$trial_type[iT]=='erasure' && tt$erasure_segments$clicks_remain[iT]==2){
      era_val[iC] <- cc[iT,era_loc[iT]];
      era_loc1[iC] <- era_loc[iT];
      iC <- iC + 1;
    }  
    if (tt$erasure_segments$trial_type[iT]=='attention' && (tt$erasure_segments$trial[iT] %% 10) == 0){
      att_loc[iD] <- era_loc[iT];
      iD <- iD + 1;
    }
  }
  
  options("encoding" = "UTF-8")
  dq0 <- NULL;
  dh0 <- NULL;
  nR <- length(era_loc1);
  for (iR in 1:nR){
    if (iR==1){
      dq0 <- paste0('<list era_loc>\n/ items = (',as.integer(era_loc1[iR]),',');
      dh0 <- paste0('<list era_val>\n/ items = (',as.character(era_val[iR]),',');
    } else if (iR > 1 && iR < nR){
      dq0 <- paste0(dq0,as.integer(era_loc1[iR]),',');
      dh0 <- paste0(dh0,as.character(era_val[iR]),',');
    } else if (iR==nR){
      dq0 <- paste0(dq0,as.integer(era_loc1[iR]),')\n/ selectionrate = always\n/ selectionmode = values.era_index;\n</list>')
      dh0 <- paste0(dh0,as.character(era_val[iR]),')\n/ selectionrate = always\n/ selectionmode = values.era_index; \n</list>')
    }
  }
  write.table(dq0,paste0('era_loc-',as.character(df$iteration),'.txt'),row.names=F,col.names=F,quote=F)
  write.table(dh0,paste0('era_val-',as.character(df$iteration),'.txt'),row.names=F,col.names=F,quote=F)
  options("encoding" = "native.enc") # change encoding back to native
  
  options("encoding" = "UTF-8")
  dq0 <- NULL;
  nR <- length(att_loc);
  for (iR in 1:nR){
    if (iR==1){
      dq0 <- paste0('<list att_loc>\n/ items = (',as.integer(att_loc[iR]),',');
    } else if (iR > 1 && iR < nR){
      dq0 <- paste0(dq0,as.integer(att_loc[iR]),',');
    } else if (iR==nR){
      dq0 <- paste0(dq0,as.integer(att_loc[iR]),')\n/ selectionrate = always\n/ selectionmode = values.att_index;\n</list>')
    }
  }
  write.table(dq0,paste0('att_loc-',as.character(df$iteration),'.txt'),row.names=F,col.names=F,quote=F)
  options("encoding" = "native.enc") # change encoding back to native
  
  
}
