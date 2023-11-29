library(R6)
library(tidyverse)

<<<<<<< Updated upstream
test_high_uncertainty = TRUE
=======
get_max_erasures = TRUE;
>>>>>>> Stashed changes

if (sum(stringr::str_detect(Sys.info(), "andypapale"))>1)  {
  base_dir <- "~/clock2/code/"
  output_dir <- "~/clock2/simulations"
  
} else {
  base_dir <- "~/code/clock2/code/"
  output_dir <- "~/code/clock2/simulations"
  j <- as.numeric(paste0(Sys.getenv("sourcefilestart")))
  j = 1
  setwd(output_dir)
  files <- list.files(pattern = "grid_")
  df <- data.table::fread(file = files[j])
}

iterate_sim <- function(df, bump_prominence, ncenters, centers, values, width_sd, i, j) {
  if (sum(stringr::str_detect(Sys.info(), "andypapale"))>1 && get_max_erasures==FALSE)  {
    contingency <- vm_circle_contingency(centers = c(centers, bump_center), weights = c(values, bump_value), widths = rep(width_sd, ncenters + 1), units = "radians")
    tt <- troll_world$new(n_trials=ntrials, values=contingency$get_wfunc(), drift_sd=1)
    tt$apply_flex(high_avg = 1, high_spread = 0, low_avg = df$low_avg[i], spread_max = 100, jump_high = T)
    tt$setup_erasure_blocks(disappear_clicks = 2, timeout_trials = 2, block_length = df$block_length[i])
<<<<<<< Updated upstream
    
    if (test_high_uncertainty){
      sceptic_agent <- scepticc$new(n_basis=12, n_points=200, contingency=tt)
      sceptic_agent$alpha <- alpha <- df$alpha[i]
      sceptic_agent$beta <- beta <- df$beta[i]
      sceptic_agent$gamma <- gamma <- df$gamma[i]
      sceptic_agent$epsilon_u <- epsilon_u <- df$epsilon_u[i]
      # new seed for the agent from teh expanded set
      set.seed(df$seed[i])
      sceptic_agent$run_contingency(optimize = FALSE)
    }
    
=======
  } else if (sum(stringr::str_detect(Sys.info(), "andypapale"))>1 && get_max_erasures==TRUE) {
    contingency <- vm_circle_contingency(centers = c(centers, bump_center), weights = c(values, bump_value), widths = rep(width_sd, ncenters + 1), units = "radians")
    tt <- troll_world$new(n_trials=ntrials, values=contingency$get_wfunc(), drift_sd=1)
    tt$apply_flex(high_avg = 1, high_spread = 0, low_avg = df$low_avg[i], spread_max = 100, jump_high = T)
    tt$setup_erasure_blocks(disappear_clicks = 2, timeout_trials = 2, block_length = df$block_length[i])
    sceptic_agent <- scepticc$new(n_basis=12, n_points=200, contingency=tt)
    sceptic_agent$alpha <- alpha <- df$alpha[i]
    sceptic_agent$beta <- beta <- df$beta[i]
    sceptic_agent$gamma <- gamma <- df$gamma[i]
    sceptic_agent$epsilon_u <- epsilon_u <- df$epsilon_u[i]
    set.seed(df$seed[i])
    learning_history <- sceptic_agent$run_contingency(optimize = FALSE)
>>>>>>> Stashed changes
  } else {
    set.seed(df$iteration[i])
    ncenters <- 9 # how many gaussians there are
    mean_val <- 10 # mean reward rate
    sd_val <- 2 # standard deviation of reward / range of rewards
    centers <- sample(seq(0, 2*pi, by = pi/20), ncenters, replace = FALSE) # line up gaussians here
    values <- sample(truncnorm::rtruncnorm(ncenters, a = 0, mean = mean_val, sd = sd_val))
    width_sd <- 20 # fixed, how wide are the underlying Gaussians
    sanity_checks = T # diagnostic plots inside simulation loop
    ntrials = 300
    cat(sprintf("In loop i: %d, j: %d\n", i, j), file = "run_log.txt", append=T)
    # set up contingency
    bump_prominence <- 10
    bump_value <- mean_val * bump_prominence
    bump_center <- sample(seq(0, 2*pi, by = pi/20), 1, replace = FALSE)
    setwd(base_dir)
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
    sceptic_agent <- scepticc$new(n_basis=12, n_points=200, contingency=tt)
    sceptic_agent$alpha <- alpha <- df$alpha[i]
    sceptic_agent$beta <- beta <- df$beta[i]
    sceptic_agent$gamma <- gamma <- df$gamma[i]
    sceptic_agent$epsilon_u <- epsilon_u <- df$epsilon_u[i]
    # new seed for the agent from teh expanded set
    set.seed(df$seed[i])
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
      u = trial_type == "erasure" & segment_shown)
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
    # u_sampled <- sum(d$in_segment)
    # browser()
  }
  if (sum(stringr::str_detect(Sys.info(), "andypapale"))>1)  {
    results <- tt
    return(results)
  } else {
    results <- as.data.frame(cbind(df[i,], r))
    return(results)
    # df$tt[i] <- tt # for now, don't save the actual contingency
  }
}
setwd(base_dir)
source("von_mises_basis.R")
source("clock2_troll_world.R")
source("scepticc.R")


if (sum(stringr::str_detect(Sys.info(), "andypapale"))>1)  {
} else {
  
  d_list <- list()
  for(i in seq(nrow(df))) {
    #for(i in seq(2)) {  # testing 
    d <- tryCatch(suppressMessages(iterate_sim(df, bump_prominence, ncenters, centers, values, width_sd, i, j)),
                  error = function(e) {
                    print(e)
                    save.image(file=sprintf("iterate_sim_error_state_%d_%d.RData", i, j))
                    return(NULL)
                  })
    if (is.null(d)) {
      print(paste0("Bad seed, low_avg=", df$low_avg[i], " iteration =", df$iteration[i],
                   " seed=", df$seed[i]))
      break}
    d_list[[i]] <- d
  }
  out_df <- data.table::rbindlist(d_list)
  setwd(output_dir)
  data.table::fwrite(out_df, file = paste0(j, "_iteration_crc.csv"))
}