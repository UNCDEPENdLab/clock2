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
require(utils)

# set up basic contingency
ncenters <- 9 # how many gaussians there are
mean_val <- 10 # mean reward rate
sd_val <- 2 # standard deviation of reward / range of rewards
centers <- sample(seq(0, 360, by = 10), ncenters, replace = FALSE) # line up gaussians here
values <- sample(truncnorm::rtruncnorm(ncenters, a = 0, mean = mean_val, sd = sd_val))
width_sd <- 20 # fixed, how wide are the underlying Gaussians
sanity_checks = F # diagnostic plots inside simulation loop

bump_prominence <- 8 # bump will always be higher, but it will change
bump_value <- mean_val * bump_prominence
# bump_center <- sample(seq(0, 360, by = 10), 1, replace = FALSE)
bump_center <- 130
# VM version
# contingency <- vm_circle_contingency(centers = c(centers, bump_center), weights=c(values, bump_value), widths = rep(width_sd, ncenters+1))

# test radians version
centers <- (pi / 180) * centers
width_sd <- (pi / 180) * width_sd
bump_center <- (pi / 180) * bump_center
# 
contingency <- vm_circle_contingency(centers = c(centers, bump_center), weights = c(values, bump_value), widths = rep(width_sd, ncenters + 1), units = "radians")

ntrials = 400
# prepare input df
idf <- expand.grid(alpha = c(.01, .05, .1, .2, .3, .4, .5), gamma = seq(.1, .9, 10),
                      beta = c(1, 2, 5, 10, 20, 40), u_sample_p = seq(.05, .5, .05), low_avg = seq(10, 50, 10),
                      block_length = c(3, 6, 9, 12, 15, 18, 21), drift = c(2, 4, 8, 16))
# since child environments inside workers do not inherit R6 objects, source them once per worker
# thus, need to loop over only a hundred dataframes
idf_list <- idf %>% group_split(alpha, beta)

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


results <- foreach(j = seq_len(length(idf_list)), .packages=c("R6", "tidyverse"), #.export = c("vm_bf", "rbf_set"),
                 .combine='rbind') %dopar% {
                   {setTxtProgressBar(pb, j)} # update progress bar
                   # suppressMessages(library(pkgs[j], character.only = TRUE))
                   source("code/von_mises_basis.R")
                   source("code/clock2_troll_world.R")
                   source("code/scepticc.R")
                   df <- idf_list[[j]]
                   df$r <- NA
                   df_list <- list()
                   for (i in seq_len(nrow(df)/10)) {
                   # set up contingency
                   tt <- troll_world$new(n_trials=ntrials, values=contingency$get_wfunc(), drift_sd=df$drift[i])
                   tt$apply_flex(high_avg = 1, high_spread = 0, low_avg = df$low_avg[i], spread_max = 100, jump_high = T)
                   # plot(tt$get_starting_values())
                   # plot(tt$spread)
                   sceptic_agent <- scepticc$new(n_basis=12, n_points=200, contingency=tt)
                   sceptic_agent$alpha <- df$alpha[i]
                   sceptic_agent$beta <- df$beta[i]
                   sceptic_agent$gamma <- df$gamma[i]
                   learning_history <- sceptic_agent$run_contingency(optimize = FALSE)
                   h <- sceptic_agent$get_entropy_history()
                   spread <- tt$spread
                   d <- cbind(learning_history, h, spread)
                   if (sanity_checks) {
                   v <- tt$get_values_matrix(type = "objective")
                   vmax <- apply(v, 1, which.max)
                   plot(vmax)
                   ggplot(dd) + geom_line( aes(x=trial, y=h*2)) + geom_point(aes(x = trial, y = choice, color = outcome)) + 
                     scale_color_viridis_b()
                   
                   ggplot(dd) + geom_line( aes(x=trial, y=h*50)) + geom_point(aes(x = trial, y = tt$spread, color = outcome)) + 
                     scale_color_viridis_b()
                   }
                   df$r[i] <- cor(lag(d$spread,10), d$h, use = "complete.obs")
                   df$tt[i] <- tt
                   }
                   return(df)
                 }
  
  # examine results
df <- results
ggplot(df, aes(as.factor(low_avg), r, color = as.factor(alpha))) + geom_boxplot() + facet_wrap(~drift)
ggplot(df, aes(low_avg, r)) + geom_jitter()  
