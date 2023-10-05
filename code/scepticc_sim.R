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

# set up basic contingency
ncenters <- 9 # how many gaussians there are
mean_val <- 10 # mean reward rate
sd_val <- 2 # standard deviation of reward / range of rewards
centers <- sample(seq(0, 2*pi, by = pi/20), ncenters, replace = FALSE) # line up gaussians here
values <- sample(truncnorm::rtruncnorm(ncenters, a = 0, mean = mean_val, sd = sd_val))
width_sd <- 20 # fixed, how wide are the underlying Gaussians
sanity_checks = F # diagnostic plots inside simulation loop

ntrials = 300
niterations = 100 # aiming for 100
# prepare input df, starting on a coarse grid
idf <- expand.grid(alpha = c(.05, .1, .25, .5), gamma = c(.1, .25, .5, .9),                 # model params
                   beta = c(1, 5, 20, 40), epsilon_u = c(0.01, 0.0833, 0.33, 0.99), # 0.0833 is at chance
                   low_avg = seq(10, 50, 10),                                               # design vars
                   block_length = c(10, 20, 30), drift = c(2, 4, 8, 16), bump_prom = c(8, 10, 15),
                   iteration = seq(niterations))
# since child environments inside workers do not inherit R6 objects, source them once per worker
# thus, need to loop over only a hundred dataframes
idf_list <- idf %>% group_split(alpha, beta, gamma)

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
pb <- txtProgressBar(0, max = length(idf_list), style = 3)


results <- foreach(j = seq_len(length(idf_list)), 
                   # .packages=c("R6", "tidyverse"), #.export = c("vm_bf", "rbf_set"),
                   # .combine='rbind') %dopar% {
                   # .options.future = list(packages=c("R6", "tidyverse")),
                   # .combine='rbind') %dofuture% {
                     .combine='rbind') %do% {
                     # suppressMessages(library(pkgs[j], character.only = TRUE))
                     # setwd("~/code/clock2")
                     # source("code/von_mises_basis.R")
                     # source("code/clock2_troll_world.R")
                     # source("code/scepticc.R")
                     df <- idf_list[[j]]
                     df$r <- NA
                     # df$u_sampled <- NA
                     for (i in seq_len(nrow(df))) {
                       # set up contingency
                       bump_prominence <- df$bump_prom[i]
                       bump_value <- mean_val * bump_prominence
                       bump_center <- sample(seq(0, 2*pi, by = pi/20), 1, replace = FALSE)
                       contingency <- vm_circle_contingency(centers = c(centers, bump_center), weights = c(values, bump_value), widths = rep(width_sd, ncenters + 1), units = "radians")
                       tt <- troll_world$new(n_trials=ntrials, values=contingency$get_wfunc(), drift_sd=df$drift[i])
                       tt$apply_flex(high_avg = 1, high_spread = 0, low_avg = df$low_avg[i], spread_max = 100, jump_high = T)
                       tt$setup_erasure_blocks(disappear_clicks = 2, timeout_trials = 1)
                       if (!is.null(tt)) {message("Generated contingency")} else {message("Trollworld failed")}
                       # plot(tt$spread)
                       sceptic_agent <- scepticc$new(n_basis=12, n_points=200, contingency=tt)
                       sceptic_agent$alpha <- df$alpha[i]
                       sceptic_agent$beta <- df$beta[i]
                       sceptic_agent$gamma <- df$gamma[i]
                       sceptic_agent$epsilon_u <- df$epsilon_u[i]
                       learning_history <- sceptic_agent$run_contingency(optimize = FALSE)
                       if (!is.null(learning_history)) {message("Ran contingency")} else {message("Scepticc failed")}
                       h <- sceptic_agent$get_entropy_history()
                       spread <- tt$spread
                       d <- cbind(learning_history, h, spread)
                       # if (sanity_checks) {
                       #   v <- tt$get_values_matrix(type = "objective")
                       #   vmax <- apply(v, 1, which.max)
                       #   plot(vmax)
                       #   ggplot(d) + geom_line( aes(x=trial, y=h*2)) + geom_point(aes(x = trial, y = choice, color = outcome)) + 
                       #     scale_color_viridis_b()
                       #   ggplot(d) + geom_line( aes(x=trial, y=h*50)) + geom_point(aes(x = trial, y = tt$spread, color = outcome)) + 
                       #     scale_color_viridis_b()
                       # } # for debugging only
                       df$r[i] <- cor(lag(d$spread,10), d$h, use = "complete.obs")
                       df$u_sampled[i] <- sum(d$in_segment)
                       # df$tt[i] <- tt # for now, don't save the actual contingency
                     }
                     setwd("~/code/clock2/simulations")
                     data.table::fwrite(df, file = paste0(j, ".csv"))
                     {setTxtProgressBar(pb, j)} # update progress bar
                     return(df)
                   }
stopCluster(cl)
# examine results
df <- results
ggplot(df, aes(as.factor(low_avg), r, color = as.factor(alpha))) + geom_boxplot() + facet_wrap(~drift)
ggplot(df, aes(low_avg, r)) + geom_jitter()  
