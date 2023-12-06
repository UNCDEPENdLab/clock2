# visualize winning contingencies across agent parameters


setwd("~/code/clock2")
source("code/von_mises_basis.R")
source("code/clock2_troll_world.R")
source("code/scepticc.R")
library(tidyverse)
require(utils)

# set up basic contingency
# niterations = 100 # aiming for 100
# prepare input df, starting on a coarse grid
# idf <- expand.grid(alpha = c(.1, .25, .5), gamma = c(.1, .5, .9),                 # model params
#                    beta = c(1, 5), # at very high betas, h and u are decorrelated, no need to test
#                    epsilon_u = c(0.33, 0.99), # 0.0833 is at chance, low correlation -- not worth testing
#                    low_avg = c(10, 20), # higher was generally worse                    # design vars
#                    block_length = c(10), # block length > 15 had higher correlations, not worth testing
#                    #drift = c(1, 2, 4), bump_prom = c(8, 10, 15),
#                    iteration = 10:niterations)


 top100 = c(4931, 5608, 5781, 7480, 5340, 2637, 5871, 5094, 6220, 720, 1269, 5903, 4505, 567, 5173, 1135, 1506, 7323, 7347, 
            283, 1291, 5966, 5655, 1041, 7315, 124, 2534, 4205, 5815, 2619, 430, 7527, 5299, 5438, 1337, 4025, 152, 5213, 6729, 4938, 
            1250, 4563, 2488, 6644, 5288, 1070, 1803, 136, 155, 5646, 727, 6227, 6108, 1639, 5307, 1764, 7933, 6583, 6578, 6896, 1752, 
            6176, 5877, 531, 6252, 6342, 4222, 5969, 6623, 7379, 1246, 7432, 5688, 4175, 5413, 1155, 2825, 7564, 6749, 1464, 5436, 294, 
            4998, 1345, 509, 448, 549, 5643, 6182, 640, 4522, 5293, 7467, 1482, 165, 7590, 2942, 4325, 7751, 163)
iterations <- top100[1:10]
# find more good seeds for contingencies
df <- expand.grid(alpha = c(.1), gamma = c(.3),                 # model params
                  beta = c(10), # at very high betas, h and u are decorrelated, no need to test
                  epsilon_u = c(0.3), # 0.0833 is at chance, low correlation -- not worth testing
                  block_length = c(10), # block length > 15 had higher correlations, not worth testing
                  low_avg = c(10),
                  # iteration = c(1:1000),
                  iteration = iterations,
                  #drift = c(1, 2, 4), bump_prom = c(8, 10, 15),
                  seed = 1)
print(paste0("Evaluating ", nrow(df), " rows"))

for (i in 1:nrow(df)) {
  set.seed(df$iteration[i])
  j = 1
  ncenters <- 9 # how many gaussians there are
  mean_val <- 10 # mean reward rate
  sd_val <- 2 # standard deviation of reward / range of rewards
  centers <- sample(seq(0, 2*pi, by = pi/20), ncenters, replace = FALSE) # line up gaussians here
  values <- sample(truncnorm::rtruncnorm(ncenters, a = 0, mean = mean_val, sd = sd_val))
  width_sd <- 20 # fixed, how wide are the underlying Gaussians
  sanity_checks = F # diagnostic plots inside simulation loop
  ntrials = 300
  cat(sprintf("In loop i: %d, j: %d\n", i, j), file = "run_log.txt", append=T)
  # set up contingency
  bump_prominence <- 10
  bump_value <- mean_val * bump_prominence
  bump_center <- sample(seq(0, 2*pi, by = pi/20), 1, replace = FALSE)
  contingency <- tryCatch(vm_circle_contingency(centers = c(centers, bump_center), weights = c(values, bump_value), widths = rep(width_sd, ncenters + 1), units = "radians"),
                          error = function(e) {
                            print(e)
                            save.image(file=sprintf("contingency_error_state_%d_%d.RData", i, j))
                            return(NULL)
                          })
  tt <- tryCatch(troll_world$new(n_trials=ntrials, values=contingency$get_wfunc(), drift_sd=1),
                 error = function(e) {
                   print(e)
                   save.image(file=sprintf("troll_world_error_state_%d_%d.RData", i, j))
                   return(NULL)
                 })
  tryCatch(tt$apply_flex(high_avg = 1, high_spread = 0, low_avg = df$low_avg[i], spread_max = 100, jump_high = T),
           error = function(e) {
             print(e)
             save.image(file=sprintf("flex_error_state_%d_%d.RData", i, j))
             return(NULL)
           })
  tryCatch(tt$setup_erasure_blocks(disappear_clicks = 2, timeout_trials = 2, block_length = df$block_length[i]),
           error = function(e) {
             print(e)
             save.image(file=sprintf("erasure_error_state_%d_%d.RData", i, j))
             return(NULL)
           })
  if (!is.null(tt)) {message("Generated contingency")} else {message("Trollworld failed")}
  # plot(tt$spread)
  # new seed for the agent from the expanded set
  sceptic_agent <- scepticc$new(n_basis=12, n_points=200, contingency=tt, seed = df$seed[i])
  sceptic_agent$alpha <- alpha <- df$alpha[i]
  sceptic_agent$beta <- beta <- df$beta[i]
  sceptic_agent$gamma <- gamma <- df$gamma[i]
  sceptic_agent$epsilon_u <- epsilon_u <- df$epsilon_u[i]
  learning_history <- tryCatch(sceptic_agent$run_contingency(optimize = FALSE),
                               error = function(e) {
                                 print(e)
                                 save.image(file=sprintf("error_state_%d_%d.RData", i, j))
                                 return(NULL)
                               })
  if (!is.null(learning_history)) {message("Ran contingency")} else {message("Scepticc failed")}
  h <- sceptic_agent$get_entropy_history()
  spread <- tt$spread
  d <- cbind(learning_history, h, spread)
  d <- d %>% inner_join(tt$erasure_segments) %>% mutate(
    u = trial_type == "erasure" & segment_shown,
    a = trial_type == "attention" & segment_shown)
  # browser()
  # if (sanity_checks) {
  #   v <- tt$get_values_matrix(type = "objective")
  #   vmax <- apply(v, 1, which.max)
  #   plot(vmax)
  #   ggplot(d) + geom_line( aes(x=trial, y=h*2)) + geom_point(aes(x = trial, y = choice, color = outcome)) +
  #     scale_color_viridis_b()
  #   ggplot(d) + geom_line( aes(x=trial, y=h*50)) + geom_point(aes(x = trial, y = tt$spread, color = outcome)) +
  #     scale_color_viridis_b()
  # } # for debugging only
  r <- cor(d$h, d$u, use = "complete.obs", method = "spearman")
  # ra <- cor(d$h, d$a, use = "complete.obs", method = "spearman")
  # u_sampled <- sum(d$in_segment)
  # browser()
  
  v <- tt$get_values_matrix(type = "objective")
  vmax <- apply(v, 1, which.max)
  d$vmax <- vmax/180*pi
  # plot(vmax)
  setwd("~/code/clock2/simulations/plots")
  pdf(paste0(df$iteration[i], "_", i, ".pdf"), height = 6, width = 10)
  print(ggplot(d) + geom_line( aes(x=trial, y=h*3)) + geom_point(aes(x = trial, y = choice, color = outcome)) +
          scale_color_viridis_b() + geom_segment(aes(x = trial, xend = trial, y = segment_min, yend = segment_max, alpha = u)) + 
          theme_minimal() + ylab("Location") + geom_text(aes(x = -5, y = 9, label = "Entropy"), angle = 90) +
          geom_text(aes(trial, vmax, label = "$", size = spread), color = "red", alpha = 0.2) +
          geom_label(aes(x = 200, y = 11, label = sprintf(paste0("epsilon_u=", sceptic_agent$epsilon_u, ", gamma=", sceptic_agent$gamma,  ", alpha=", sceptic_agent$alpha, ", beta=", sceptic_agent$beta,  ", r=",round(r,2), ", tot. points=", round(sum(d$outcome),0))))))
  dev.off()
  # Sys.sleep(1)
}
