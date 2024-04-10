# 2024-01-03 AndyP
# Testing rng seeds





base_dir <- "~/clock2/code/"
rob_grid <- expand.grid(alpha = c(0.2), gamma = c(0.1),                 # model params
                        beta = c(1), # at very high betas, h and u are decorrelated, no need to test
                        epsilon_u = c(0.9999), # 0.0833 is at chance, low correlation -- not worth testing
                        block_length = c(10), # block length > 15 had higher correlations, not worth testing
                        low_avg = c(10),
                        iteration = c(152),
                        #drift = c(1, 2, 4), bump_prom = c(8, 10, 15),
                        seed = 1)
setwd(base_dir)
source("von_mises_basis.R")
source("clock2_troll_world.R")
source("scepticc.R")
set.seed(rob_grid$iteration)
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
tt$apply_flex(high_avg = 1, high_spread = 0, low_avg = rob_grid$low_avg, spread_max = 100, jump_high = T)
tt$setup_erasure_blocks(disappear_clicks = 2, timeout_trials = 2, block_length = rob_grid$block_length)
sceptic_agent <- scepticc$new(n_basis=12, n_points=200, contingency=tt)
sceptic_agent$alpha <- alpha <- rob_grid$alpha
sceptic_agent$beta <- beta <- rob_grid$beta
sceptic_agent$gamma <- gamma <- rob_grid$gamma
sceptic_agent$epsilon_u <- epsilon_u <- rob_grid$epsilon_u
set.seed(rob_grid$seed)
learning_history <- sceptic_agent$run_contingency(optimize = FALSE)
cc <- round(tt$get_values_matrix(type = 'objective', quiet = F),0)
for (r in 1:nrow(cc)) {
   plot(cc[r,])
   print(r)
   Sys.sleep(1)
}
