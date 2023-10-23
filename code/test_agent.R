library(tidyverse)
if (sum(stringr::str_detect(Sys.info(), "Alex|alexdombrovski"))>1) {
  setwd("~/code/clock2/")
} else if (sum(stringr::str_detect(Sys.info(), "andypaple"))>1)  {
  setwd('/clock2/')
} else {
  setwd("~/Data_Analysis/clock2")}
source("code/von_mises_basis.R")
source("code/clock2_troll_world.R")
source("code/scepticc.R")
ncenters <- 9 # how many gaussians there are
mean_val <- 10 # mean reward rate
sd_val <- 2 # standard deviation of reward / range of rewards
centers <- sample(seq(0, 360, by = 10), ncenters, replace = FALSE) # line up gaussians here
values <- sample(truncnorm::rtruncnorm(ncenters, a = 0, mean = mean_val, sd = sd_val))
width_sd <- 20 # fixed, how wide are the underlying Gaussians

bump_prominence <- 8 # bump will always be higher, but it will change
bump_value <- mean_val * bump_prominence
bump_center <- sample(seq(0, 360, by = 10), 1, replace = FALSE)

animated_plot <- F

# VM version
# contingency <- vm_circle_contingency(centers = c(centers, bump_center), weights=c(values, bump_value), widths = rep(width_sd, ncenters+1))

# test radians version
centers <- (pi / 180) * centers
width_sd <- (pi / 180) * width_sd
bump_center <- (pi / 180) * bump_center
# 
contingency <- vm_circle_contingency(centers = c(centers, bump_center), weights = c(values, bump_value), widths = rep(width_sd, ncenters + 1), units = "radians")

# a <- rbf$new(value=10, value_sd=1, center=100, width_sd = 20)
# vv <- a$get_tvec()

# contingency <- rbf_set$new(elements = gg)
plot(contingency$get_wfunc())

contingency$get_centers()
contingency$get_weights()



# used in prolific
# tt <- troll_world$new(n_trials=345, values=contingency$get_wfunc(), drift_sd=4) # set up troll world
# tt$apply_flex(high_avg = 1, high_spread = 0, jump_high=FALSE)
# tt$setup_erasure_blocks()
# plot(tt$spread) # prominence of bump vs floor over trials, shows switches, bump drifts in Gaussian random walk

# stable 100-trial contingency for model validation
tt <- troll_world$new(n_trials=300, values=contingency$get_wfunc(), drift_sd=1)
tt$apply_flex(high_avg = 1, high_spread = 0, spread_max = 100, jump_high = FALSE)
tt$setup_erasure_blocks(disappear_clicks = 2, timeout_trials = 1)
sceptic_agent <- scepticc$new(n_basis=12, n_points=360, contingency=tt, beta = 40, epsilon_u = .1)
df <- sceptic_agent$run_contingency()
tt$get_choices()
plot(tt$spread)
df$h <- sceptic_agent$get_entropy_history()
u <- tt$get_choices() %>% mutate(u = trial_type=="erasure" & segment_shown)
df <- df %>% inner_join(u)
sum(df$in_segment)
sum(df$u)

ggplot(df) + geom_line( aes(x=trial, y=h*2)) + geom_point(aes(x = trial, y = choice, alpha = outcome)) + 
  scale_color_viridis_d() + geom_errorbar(aes(x = trial, y = (segment_min + segment_max)/2, ymin = segment_min, ymax = segment_max, color = trial_type)) +
  theme_minimal()
cor.test(as.numeric(df$u), df$h)

ggplot(df) + geom_line( aes(x=trial, y=h*50)) + geom_point(aes(x = trial, y = tt$spread, color = outcome)) + 
  scale_color_viridis_b()

# stable 100-trial contingency for model validation
# tt <- troll_world$new(n_trials=100, values=contingency$get_wfunc(), drift_sd=0)
# tt$apply_flex(high_avg = 1, high_spread = 0, spread_max = 100, jump_high = FALSE)
# plot(tt$spread)
# 
# plot(tt$get_starting_values())
# aa <- tt$get_values_matrix("original", quiet=F) # original values (should be constant over trials)
# aa <- tt$get_values_matrix("drift", quiet=F) # drift alone
# aa <- tt$get_values_matrix("flex", quiet=F) # flex alone
aa <- tt$get_values_matrix("objective", quiet=F) # all manipulations, matrix of expected values, p_reward =0.7 fixed for now

# aa[1:5, 1:10]
# tt$reset_counter()
# tt$get_next_values()[1:10]

loc <- round(tt$get_pvec(), 2)


if (animated_plot) {
  
  for (ii in 1:nrow(aa)) {
    if (tt$erasure_segments$trial_type[ii] == "erasure") {
      ss <- paste(", er:", round(tt$erasure_segments[ii,"segment_min"], 2))
    } else {
      ss <- ""
    }
    plot(aa[ii,], type="l", main=paste("Trial", ii, "epoch", tt$epoch[ii], ss), ylim = range(aa), xaxt='n')
    axis(side = 1, at = seq(1, 360, by=30), labels= loc[seq(1, 360, by=30)])
    Sys.sleep(.2)
    
  }
}
  
#inq_val <- round(t(tt$get_values_matrix()),0)
#write.csv(inq_val,file = '2023-09-12-Design-File-asMatrix.csv')

values <- data.frame(round(t(tt$get_values_matrix())),0)
#values <- values %>% mutate(timepoint = row_number()) %>% rowwise() %>% pivot_longer(cols = starts_with("X"), names_to = "trial")
values <- values %>% mutate(timepoint = row_number()) %>% rowwise() %>% pivot_longer(cols = starts_with("X"), names_to = "trial") %>%
  mutate(trial = extract_numeric(trial)) %>% group_by(trial) %>% summarise(vmax = max(value),
                                                                           vmax_location = timepoint[which.max(value)])
plot(values$vmax_location)
# in progress
#tt$erase_segment(30, trial=10)

# for (i in 1:10) {
#   plot(tt$get_next_values())
#   Sys.sleep(.3)
# }

# placement of erasures
# - to maintain entropy manipulation, preserve top and bottom 5% of the unerased segments
# - erasure indicators fade instantly (one trial) because the values then enter the overall value distribution
# and become a mix of softmax stochastic and uncertainty-driven. Only the first trial is clear on process.
# WIP
#tt$erase_segment()


## completely static contingency
tt <- troll_world$new(n_trials=100, values=contingency$get_wfunc(), drift_sd=0)
#plot(tt$spread)
tt$apply_flex(high_avg = 1, high_spread = 0, jump_high=FALSE)
tt$setup_erasure_blocks(timeout_trials = 4)
plot(tt$get_starting_values())
tt$reset_counter()
# aa <- tt$get_values_matrix("original", quiet=F) # original values (should be constant over trials)
# aa <- tt$get_values_matrix("drift", quiet=F) # drift alone
# aa <- tt$get_values_matrix("flex", quiet=F) # flex alone
aa <- tt$get_values_matrix("objective", quiet=F) # all manipulations

# aa[1:5, 1:10]
# tt$reset_counter()
# tt$get_next_values()[1:10]


## try out sceptic agent
sceptic_agent <- scepticc$new(n_basis=12, n_points=360, contingency=tt)
sceptic_agent$beta <- 5
learning_history <- sceptic_agent$run_contingency(optimize = FALSE)
world <- sceptic_agent$get_contingency()
df <- world$get_choices()
world$erasure_segments


#sceptic_agent$update_weights(tau=pi/2, outcome=100)
# sceptic_agent$emit_choice()
#plot(sceptic_agent$get_choice_probs(), type="l")
#sceptic_agent$get_weights()

# dynamic contingency
ttd <- troll_world$new(n_trials=300, values=contingency$get_wfunc(), drift_sd=5)
ttd$apply_flex(high_avg = 1, high_spread = 0, low_avg = 40, spread_max = 100, jump_high = T)
plot(ttd$get_starting_values())
plot(ttd$spread)
sceptic_agent <- scepticc$new(n_basis=12, n_points=200, contingency=ttd)


sceptic_agent$beta <- 5
# obj_f <- data.frame(choice = seq(0, 2*pi, length.out=360), value = tt$get_starting_values())


learning_history <- sceptic_agent$run_contingency(optimize = FALSE)
world <- sceptic_agent$get_contingency()
df <- world$get_choices()

h <- sceptic_agent$get_entropy_history()
spread <- ttd$spread
v <- ttd$get_values_matrix(type = "objective")
vmax <- apply(v, 1, which.max)
plot(vmax)
df <- cbind(learning_history, h, spread)
# ggplot(obj_f, aes(x=choice, y=value)) + geom_line() +
#   geom_point(data = learning_history, mapping=aes(x=choice, y=outcome, color=trial))
cor.test(df$h, lag(df$spread, 10))
ggplot(df, aes(x=trial, y=outcome, color=choice)) + geom_point() +
  stat_smooth()

ggplot(df) + geom_line( aes(x=trial, y=h*2)) + geom_point(aes(x = trial, y = choice, color = outcome)) + 
  scale_color_viridis_b()

ggplot(df) + geom_line( aes(x=trial, y=h*50)) + geom_point(aes(x = trial, y = ttd$spread, color = outcome)) + 
  scale_color_viridis_b()
plot(lag(df$spread,n = 10), df$h)
cor.test(lag(df$spread,10), df$h)

ggplot(df, aes(x=trial, y=choice)) + geom_point() +
  stat_smooth()


## prototype cubic spline erasure
l <- seq(0, 2*pi, length.out=360)
v <- tt$tvals[1,]


dd <- rep(NA, 10000)
dp <- rep(NA, 10000)
for(ii in 1:10000) {
  vnew <- erase(v, l)
  dd[ii] <- attr(vnew, "delta")
  #dp[ii] <- attr(vnew, "delta_pct")
  dp[ii] <- attr(vnew, "pct_orig")
  #plot(vnew, type="l", main=attr(vnew, "ecent"))
  #Sys.sleep(.5)
}

for(ii in 1:30) {
  vnew <- erase(v, l)
  #dd[ii] <- attr(vnew, "delta")
  #dp[ii] <- attr(vnew, "delta_pct")
  #dp[ii] <- attr(vnew, "pct_orig")
  plot(vnew, type="l", main=attr(vnew, "ecent"))
  Sys.sleep(.5)
}


median(dd)
mean(dd)
sd(dd)
hist(dd)
hist(dp)
summary(dp)
table(dp > 0)

# pi/6 = 30 degrees
# erase <- function(v, loc, width = pi/6, pos = 0.5, pos_sd = 0, v_quantiles = c(.05, .95)) {
#   checkmate::assert_number(pos, lower=0, upper=1)
#   stopifnot(length(v) == length(loc))
#   v_new <- v
#   
#   # left-most point of erased segment
#   low_pos <- sample(seq_along(v), 1) # vector index
#   low_rad <- loc[low_pos] # location of start point in radians based on location
#   high_rad <- (low_rad + width) %% (2*pi) # wrap back onto circle
#   high_pos <- which.min(abs(loc - high_rad))
#   
#   # point within the segment where we park the resampled value
#   pos_q <- truncnorm::rtruncnorm(1, a = 0, b = 1, mean = pos, sd = pos_sd)
#   
#   # quantile of position within erased vector
#   mid_pos <- round(quantile(low_pos:high_pos, pos_q))
#   
#   
#   #mid_rad <- loc[mid_pos]
#   #low_rad <- mid_rad - width/2
#   #high_rad <- mid_rad + width/2
#   #low_pos <- which.min(abs(loc - low_rad))
#   
#   v_new[low_pos:high_pos] <- NA
#   
#   v_low <- quantile(v, v_quantiles[1])
#   v_high <- quantile(v, v_quantiles[2])
#   v_elig <- v[v > v_low & v < v_high]
#   v_point <- sample(v_elig, 1) # sample only from eligible quantiles
#   # place this value sample/point at the mid-point location for interpolation
#   v_new[mid_pos] <- v_point
#   
#   # drop in cubic spline
#   v_new <- zoo::na.spline(v_new)
#   attr(v_new, "ecent") <- mid_pos
#   return(v_new)
# }


erase <- function(v, loc, width = pi/6, v_quantiles = c(.25, .75), r_quantiles = c(.25, .75)) {
  checkmate::assert_number(pos, lower=0, upper=1)
  stopifnot(length(v) == length(loc))
  v_new <- v
  
  # which values are eligible for erasure?
  qs <- quantile(v, r_quantiles)
  pp <- which(v >= qs[1] & v <= qs[2])
  
  # mid-point of erased segment
  mid_pos <- sample(pp, 1) # vector index
  mid_rad <- loc[mid_pos]
  low_rad <- (mid_rad - width/2) %% (2*pi) # location of start point in radians based on location
  high_rad <- (mid_rad + width/2) %% (2*pi) # wrap back onto circle
  low_pos <- which.min(abs(loc - low_rad))
  high_pos <- which.min(abs(loc - high_rad))
  
  # erase segment with NAs
  v_new[low_pos:high_pos] <- NA
  
  # find the value for the midpoint of the erasure based on the eligible quantiles
  v_low <- quantile(v, v_quantiles[1])
  v_high <- quantile(v, v_quantiles[2])
  v_elig <- v[v > v_low & v < v_high]
  v_point <- sample(v_elig, 1) # sample only from eligible quantiles
  # place this value sample/point at the mid-point location for interpolation
  v_new[mid_pos] <- v_point
  
  # use cubic spline interpolation to connect the dots
  v_new <- zoo::na.spline(v_new)
  attr(v_new, "ecent") <- mid_pos
  attr(v_new, "delta") <- v[mid_pos] - v_point
  #attr(v_new, "delta_pct") <- ((v[mid_pos] - v_point)/v[mid_pos]) * 100
  attr(v_new, "pct_orig") <- (v_point/v[mid_pos]) * 100
  return(v_new)
}


library(ggplot2)
ggplot(learning_history, aes(x=trial, y=outcome))

test <- sceptic_agent$optimize_returns()


#system.time(aa <- replicate(100, sceptic_agent$run_contingency()))

hx <- sceptic_agent$get_func_history()
for (ii in 1:nrow(hx)) {
  plot(hx[ii,], type="l", main=paste("Trial", ii, "epoch", tt$epoch[ii]), ylim = range(aa))
  Sys.sleep(.1)
}


# returns history of basis weights
#h <- reshape2::melt(sceptic_agent$get_weight_history(), varnames=c("trial", "basis"))
v <- reshape2::melt(sceptic_agent$get_func_history(), varnames=c("trial", "timestep"))

df <- merge(learning_history, v, by="trial")

# rescale onto same scale
learning_history$outcome_norm <- scales::rescale(learning_history$outcome, to = c(range(v$value)))

library(gganimate)
g <- ggplot(df, aes(x=timestep, y=value)) + geom_line() +
  geom_point(data=learning_history, mapping = aes(x=choice, y=outcome_norm))


g + transition_states(trial) +
  ggtitle("Trial {closest_state}, phase: {tt$epoch[as.integer(closest_state)]}")


for (ii in 1:nrow(h)) {
  #plot(h[ii,], type="l", main=paste("Basis weights, trial", ii, "epoch", tt$epoch[ii]), ylim = range(h))
  # g <- ggplot(df %>% filter(trial==!!ii), aes(x=timestep, y=value)) + geom_line() +
  #    geom_point(data=learning_history[ii,,drop=F], mapping = aes(x=choice, y=outcome_norm))
  # plot(g)
  # Sys.sleep(.1)
}
