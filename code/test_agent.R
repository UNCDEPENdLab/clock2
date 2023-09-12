setwd("~/clock2")
source("code/von_mises_basis.R")
source("code/clock2_troll_world.R")
source("code/scepticc.R")
library(tidyverse)
ncenters <- 9 # how many gaussians there are
mean_val <- 15 # mean reward rate
sd_val <- 2 # standard deviation of reward / range of rewards
centers <- sample(seq(0, 360, by = 10), ncenters, replace = FALSE) # line up gaussians here
values <- sample(truncnorm::rtruncnorm(ncenters, a = 0, mean = mean_val, sd = sd_val))
width_sd <- 5 # fixed, how wide are the underlying Gaussians

bump_prominence <- 6 # bump will always be higher, but it will change
bump_value <- mean_val * bump_prominence
bump_center <- sample(seq(0, 360, by = 40), 1, replace = FALSE)

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


tt <- troll_world$new(n_trials=345, values=contingency$get_wfunc(), drift_sd=5) # set up troll world
tt$apply_flex(high_avg = 1, high_spread = 0)
plot(tt$spread) # prominence of bump vs floor over trials, shows switches, bump drifts in Gaussian random walk
plot(tt$get_starting_values())
# aa <- tt$get_values_matrix("original", quiet=F) # original values (should be constant over trials)
# aa <- tt$get_values_matrix("drift", quiet=F) # drift alone
# aa <- tt$get_values_matrix("flex", quiet=F) # flex alone
aa <- tt$get_values_matrix("objective", quiet=F) # all manipulations, matrix of expected values, p_reward =0.7 fixed for now

aa[1:5, 1:10]
tt$reset_counter()
tt$get_next_values()[1:10]

# for (ii in 1:nrow(aa)) {
#   plot(aa[ii,], type="l", main=paste("Trial", ii, "epoch", tt$epoch[ii]), ylim = range(aa))
#   Sys.sleep(.1)
# }
values <- data.frame(t(tt$get_values_matrix()))

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


## try out sceptic agent
sceptic_agent <- scepticc$new(n_basis=12, n_points=50, contingency=tt)
# sceptic_agent$update_weights(tau=pi/2, outcome=100)
# sceptic_agent$emit_choice()
# sceptic_agent$get_choice_probs()
sceptic_agent$get_weights()
learning_history <- sceptic_agent$run_contingency()


system.time(aa <- replicate(100, sceptic_agent$run_contingency()))

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
