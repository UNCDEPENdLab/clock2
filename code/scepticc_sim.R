# simulation wrapper for Trollworld and scepticc
# parameters to vary: alpha, beta, gamma, (inside SCEPTIC)
# length of low-entropy epoch, erasure block length (inside Trollworld)
# uncertainty attitude (based on behavioral data, awaiting full SCEPTIC implementation)

# questions: set ntrials to 400 and then cut later if needed?

setwd("~/code/clock2")
source("code/von_mises_basis.R")
source("code/clock2_troll_world.R")
source("code/scepticc.R")
library(tidyverse)
library(parallel)
library(doParallel)
library(doFuture)
plan(multisession)
require(utils)
parallel = F

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

iterations <- c(5888, 8395)
# find more good seeds for contingencies
rob_grid <- expand.grid(alpha = c(.5), gamma = c(.5),                 # model params
                        beta = c(1), # at very high betas, h and u are decorrelated, no need to test
                        epsilon_u = c(0.9), # 0.0833 is at chance, low correlation -- not worth testing
                        block_length = c(10), # block length > 15 had higher correlations, not worth testing
                        low_avg = c(10),
                        # iteration = c(1:1000),
                        iteration = iterations,
                        #drift = c(1, 2, 4), bump_prom = c(8, 10, 15),
                        seed = 1:10)

# since child environments inside workers do not inherit R6 objects, source them once per worker
# thus, need to loop over only a hundred dataframes
# idf_list <- idf %>% group_split(iteration, alpha)
# idf_list <- cross_join(winners, rob_grid) %>% group_split(seed, epsilon_u, iteration)

idf_list <- rob_grid %>% group_split(iteration, epsilon_u, beta, alpha, gamma)


if (parallel) {
# make cluster ----
f <- Sys.getenv('PBS_NODEFILE')
ncores <- detectCores() - 1

nodelist <- if (nzchar(f)) readLines(f) else rep('localhost', ncores)
cat("Node list allocated to this job\n")
print(nodelist)
cl <- makePSOCKcluster(nodelist, outfile='')
print(cl) ##; print(unclass(cl))
registerDoParallel(cl)
parallel::clusterExport(cl, setdiff(ls(), "cl"))

# loop over parameters ----
message("Running simulation")
# registerDoSEQ()
pb <- txtProgressBar(0, max = length(idf_list), style = 3)

}


# results <- foreach(j = seq_len(length(idf_list)), .errorhandling = "remove",
# # results <- foreach(j = seq_len(23), .errorhandling = "remove",
# .packages=c("R6", "tidyverse"), #.export = c("vm_bf", "rbf_set"),
# .combine='rbind') %dopar% {
iterate_sim <- function(df, bump_prominence, ncenters, centers, values, width_sd, i, j) {
  set.seed(df$iteration[i])
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
  results <- as.data.frame(cbind(df[i,], r))
  return(results)
  # df$tt[i] <- tt # for now, don't save the actual contingency
}

# .options.future = list(packages=c("R6", "tidyverse")),
# .combine='rbind') %dofuture% {
# .combine='rbind') %do% {
# suppressMessages(library(pkgs[j], character.only = TRUE))

# completely sequential version
big_d_list <- list()
for (j in seq(length(idf_list))) {
  setwd("~/code/clock2")
  source("code/von_mises_basis.R")
  source("code/clock2_troll_world.R")
  source("code/scepticc.R")
  df <- idf_list[[j]]
  # df$r <- NA
  # df$u_sampled <- NA
  d_list <- list()
  print(paste0(round(j/length(idf_list)*100, 1), ' %'))
  # for (i in seq_len(nrow(df))) {
  for(i in seq(nrow(df))) {
    d <- tryCatch(suppressMessages(iterate_sim(df, bump_prominence, ncenters, centers, values, width_sd, i, j)),
                  error = function(e) {
                    print(e)
                    save.image(file=sprintf("erasure_error_state_%d_%d.RData", i, j))
                    return(NULL)
                  })
    if (is.null(d)) {
      print(paste0("Bad seed, low_avg=", df$low_avg[i], " iteration =", df$iteration[i],
                   " seed=", df$seed[i]))
      break}
    d_list[[i]] <- d
  }
  big_d_list[[j]] <- data.table::rbindlist(d_list)
  setwd("~/code/clock2/simulations")
  data.table::fwrite(data.table::rbindlist(d_list), file = paste0(j, "robust.csv"))
  # {setTxtProgressBar(pb, j)} # update progress bar
  # return(data.table::rbindlist(d_list))
}

big_df <- data.table::rbindlist(big_d_list)
data.table::fwrite(data.table::rbindlist(big_d_list), file = paste0("11_oct_2023_top_20_robust_across_params.csv"))

# stopCluster(cl)

# digest earlier results
# setwd("robust200/")
# files <- list.files(pattern = "robust.csv") 
# df200 <- purrr::map_dfr(files, data.table::fread)
# setwd("../")
# files <- list.files(pattern = "robust.csv")
# df100 <- purrr::map_dfr(files, data.table::fread)
# df_all <- rbind(df100, df200)
# 
# wdf <- df_all %>% summarise(.by = c(low_avg, epsilon_u, iteration, beta, gamma, alpha), mean_r = mean(r)) %>% 
#   filter(beta == 1 & gamma == 0.5 & epsilon_u == .99 & alpha == 0.5) %>%  arrange(mean_r) %>% top_n(-20)
# 
# df_all %>% filter(iteration == 114 & low_avg == 20)
# 
# # examine results
# # df <- big_df
# df100 %>% summarise(.by = c(low_avg, epsilon_u, iteration, beta, gamma), mean_r = mean(r)) %>% 
#   filter(beta == 1 & gamma == 0.5) %>%  arrange(mean_r)
# df200 %>% summarise(.by = c(low_avg, iteration), mean_r = mean(r)) %>% arrange(mean_r) %>% top_n(-20)


df <- big_df
sdf <- df %>% summarise(.by = c(iteration), mean_r = mean(r)) %>% arrange(mean_r) %>% top_n(-20)
write_csv2(sdf, file = "winning_20_beta1_u99_gamma.5_alpha.5_1000_seeds_robust.csv")
# best low_avg, iterations
iterations <- unique(sdf$iteration)

idf <- df %>% summarise(.by = c(low_avg, iteration), mean_r = mean(r)) %>% arrange(mean_r) %>% top_n(-20)

top_50_iterations <- df %>% summarise(.by = c(low_avg, iteration), mean_r = mean(r)) %>% arrange(mean_r) %>% top_n(-50) %>% select(iteration)

# best low_avg, iterations for worst-case scenario (low beta, gamma > .1, epsilon_u = .99)
wdf <- df %>% filter(beta == 1 & gamma > .1 & epsilon_u == .99) %>% 
  summarise(.by = c(low_avg, iteration, gamma), mean_r_worst = mean(r)) %>% arrange(mean_r_worst) %>% top_n(-20) %>% inner_join(idf)
write_csv2(wdf, file = "winning_6_lowavg_seed_robust.csv")



winners <- wdf %>% select(low_avg, iteration)
write_csv2(winners, file = "winning_6_lowavg_seed.csv")


pdf("simulations_plot_worst.pdf", height = 20, width = 20)
ggplot(df, aes(low_avg, r, color = as.factor(iteration), lty = as.factor(block_length))) + geom_line() + facet_grid(alpha + gamma ~ epsilon_u + beta)
dev.off()

pdf("sim_90_iterations.pdf", height = 10, width = 50)
ggplot(sdf, aes(as.factor(iteration), mean_r, color = as.factor(epsilon_u), lty = as.factor(low_avg))) + geom_boxplot() #+ theme(legend.position = "none")
dev.off()
ggplot(df, aes(low_avg, r, color = as.factor(iteration))) + geom_smooth(method = "loess")

pdf("simulations_jitter_worst.pdf", height = 20, width = 20)
ggplot(df, aes(low_avg, r, color = as.factor(iteration), pch = as.factor(block_length))) + geom_jitter() + 
  facet_grid(alpha + gamma ~ epsilon_u + beta) + geom_hline(yintercept = 0, size = 1) + geom_hline(yintercept = 0.25, size = .5) + geom_hline(yintercept = 0.5, size = .25)
dev.off()

ggplot(df, aes(as.factor(low_avg), r, color = as.factor(iteration), lty = as.factor(block_length))) + geom_violin(draw_quantiles = .5)


ggplot(df, aes(as.factor(low_avg), r, color = as.factor(block_length))) + geom_boxplot() #+ facet_wrap(~drift)
ggplot(df, aes(as.factor(block_length), r, color = as.factor(alpha))) + geom_boxplot() + facet_wrap(~as.factor(gamma))
ggplot(df, aes(as.factor(block_length), r, color = as.factor(alpha))) + geom_boxplot() + facet_wrap(~as.factor(beta))
ggplot(df, aes(as.factor(bump_prom), r, color = as.factor(alpha))) + geom_boxplot() #+ facet_wrap(~drift)


ggplot(df, aes(low_avg, r)) + geom_jitter()  

# focus on problematic behaviors, worst-case scenario: block_length of 10-15 and low_avg 10-20 is still best
df <- big_df %>% filter(epsilon_u > 0.8 & beta < 20 & alpha <.2)

m0 <- lme4::lmer(r ~  (1|iteration), df)
summary(m0)


m1 <- lme4::lmer(r ~ (as.factor(low_avg) + epsilon_u + beta + gamma)^2 + (1|iteration), df)
car::Anova(m1, '3')
summary(m1)

lm1 <- lm(r ~ (as.factor(low_avg) + epsilon_u + alpha + beta + gamma + as.factor(iteration))^2, df)
anova(lm1)
summary(lm1)

# look at only the highest epsilon_u
lm2 <- lm(r ~ (as.factor(block_length) + as.factor(low_avg) + beta + gamma + as.factor(iteration))^2, df)
anova(lm2)
summary(lm2)

lm3 <- lm(r ~ (as.factor(block_length) + as.factor(low_avg) + as.factor(iteration)), df)
anova(lm3)
summary(lm3)

