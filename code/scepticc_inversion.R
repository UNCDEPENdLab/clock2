# repo_dir <- "~/Data_Analysis/clock2"
repo_dir <- "~/clock2"

setwd(repo_dir)
source(file.path(repo_dir, "code", "von_mises_basis.R"))
source(file.path(repo_dir, "code", "clock2_troll_world.R"))
source(file.path(repo_dir, "code", "scepticc.R"))

library(tidyverse)
library(ggplot2)
library(rstan)

results_dir <- "~/OneDrive/collected_letters/data_projects/clock_2/scepticc"; dir.create(file.path(results_dir), showWarnings = FALSE)

# flags
run_sampling <- F
run_vba <- T


# scepticc inversion in stan
ncenters <- 9 # how many gaussians there are
mean_val <- 10 # mean reward rate
sd_val <- 2 # standard deviation of reward / range of rewards
centers <- sample(seq(0, 360, by = 10), ncenters, replace = FALSE) # line up gaussians here
values <- sample(truncnorm::rtruncnorm(ncenters, a = 0, mean = mean_val, sd = sd_val))
width_sd <- 20 # fixed, how wide are the underlying Gaussians

bump_prominence <- 8 # bump will always be higher, but it will change
bump_value <- mean_val * bump_prominence
bump_center <- sample(seq(0, 360, by = 10), 1, replace = FALSE)

centers <- (pi / 180) * centers
width_sd <- (pi / 180) * width_sd
bump_center <- (pi / 180) * bump_center
# 
contingency <- vm_circle_contingency(centers = c(centers, bump_center), weights = c(values, bump_value), widths = rep(width_sd, ncenters + 1), units = "radians")
ttd <- troll_world$new(n_trials=300, values=contingency$get_wfunc(), drift_sd=5)
ttd$apply_flex(high_avg = 1, high_spread = 0, low_avg = 10, spread_max = 150, jump_high = FALSE)
ttd$setup_erasure_blocks(disappear_clicks = 2, timeout_trials = 1)
# plot(ttd$get_starting_values())
# plot(ttd$spread)

# helper functino to generate multi-subject data
generate_scepticc_data <- function(n=25, contingency = NULL) {
  require(truncnorm)
  if (!inherits(contingency, "troll_world")) stop("contingency must be troll_world object")
  
  seeds <- sample(1:1e4, size=n, replace=FALSE)
  
  alpha_m <- 0.3 # was 0.1
  alpha_sd <- 0.2
  alphas <- rtruncnorm(n, a=0, b=1, mean=alpha_m, sd=alpha_sd)
  
  beta_m <- 1 # was .1, very cool
  beta_sd <- 2
  betas <- rtruncnorm(n, a=0.001, mean=beta_m, sd=beta_sd)
  
  gamma_m <- 0.3
  gamma_sd <- 0.2
  gammas <- rtruncnorm(n, a=0, b=1, mean=gamma_m, sd=gamma_sd)
  
  epsilon_u_m <- 1/12 # 8%
  epsilon_u_sd <- 1/4
  epsilons <- rtruncnorm(n, a=0, b=1, mean=epsilon_u_m, sd=epsilon_u_sd)
  
  epsilon_u_m <- 1/12 # 8%
  epsilon_u_sd <- 1/4
  epsilon_us <- rtruncnorm(n, a=0, b=1, mean=epsilon_u_m, sd=epsilon_u_sd)
  
  epsilon_a_m <- 1/12 # 8%
  epsilon_a_sd <- 1/4
  epsilon_as <- rtruncnorm(n, a=0, b=1, mean=epsilon_a_m, sd=epsilon_a_sd)
  
  dlist <- lapply(seq_len(n), function(i) {
    ag <- scepticc$new(n_basis=12, n_points=360, contingency=contingency,
                       alpha = alphas[i], beta = betas[i], gamma = gammas[i],
                       epsilon_u = epsilon_us[i], epsilon_a = epsilon_as[i], seed=seeds[i])
    
    learning_history <- ag$run_contingency(optimize = FALSE)
    world <- ag$get_contingency()
    df <- world$get_choices()
    return(df)
    
  })
  
  ddf <- tibble(id=1:n, alphas, gammas, betas, epsilon_us, epsilon_as, data=dlist) %>%
    unnest(data)
  
}

nsubj <- 100
#ff_noerasure <- generate_scepticc_data(n = nsubj, contingency = ttd)
ff <- generate_scepticc_data(n = nsubj, contingency = ttd)

pars <- ff %>%
  group_by(id) %>%
  summarise_all(mean)

outcomes <- ff %>%
  select(id, trial, outcome) %>%
  pivot_wider(id_cols = "id", names_from = "trial", values_from = "outcome") %>%
  select(-id) %>%
  as.matrix()

choices <- ff %>%
  select(id, trial, choice_rad) %>%
  pivot_wider(id_cols = "id", names_from = "trial", values_from = "choice_rad") %>%
  select(-id) %>%
  as.matrix()

trial_types <- ff %>%
  select(id, trial, trial_type) %>%
  pivot_wider(id_cols = "id", names_from = "trial", values_from = "trial_type") %>%
  select(-id) %>%
  as.matrix() %>%
  apply(c(1,2), function(x) switch(x, "no erasure"=1, "erasure"=2, "attention"=3))

segment_shown <- ff %>%
  select(id, trial, segment_shown) %>%
  pivot_wider(id_cols = "id", names_from = "trial", values_from = "segment_shown") %>%
  select(-id) %>%
  as.matrix() %>%
  apply(c(1,2), function(x) if (is.na(x)) 0 else x) # replace NA with 0

segment_min <- ff %>%
  select(id, trial, segment_min) %>%
  pivot_wider(id_cols = "id", names_from = "trial", values_from = "segment_min") %>%
  select(-id) %>%
  as.matrix() %>%
  apply(c(1,2), function(x) if (is.na(x)) -99 else x) # replace NA with -99

segment_max <- ff %>%
  select(id, trial, segment_max) %>%
  pivot_wider(id_cols = "id", names_from = "trial", values_from = "segment_max") %>%
  select(-id) %>%
  as.matrix() %>%
  apply(c(1,2), function(x) if (is.na(x)) -99 else x) # replace NA with -99


# round choices in radians onto positions vector/points -- superseded by internal bincode function in stan
choices_to_pbins <- function(choices, bins=120) {
  d_theta <- 2 * pi / bins
  
  # since circle wraps, we want to avoid adding a redundant basis function at 0 vs. at 2*pi
  pvec <- seq(0, 2 * pi - d_theta, length.out = bins)
  choices[] <- Hmisc::cut2(choices, cuts = pvec)
  return(choices)
}

rc <- choices_to_pbins(choices)

# setup inputs to Stan
data_list <- list(
  B = 16, # 12 basis functions
  P = 90, # 90 evaluation points for value
  N = nsubj,
  T = 300,
  Tsubj = rep(300, nsubj),
  sb = 0.3, # SD of basis functions in radians
  sg = 0.3, # SD of generalization function in radians
  outcomes = outcomes,
  choices_rad = choices,
  trial_types = trial_types,
  segment_shown = segment_shown,
  segment_min = segment_min,
  segment_max = segment_max
)
browser()

niter  <- 2000
nwarmup <- 500
nchain <- 4
ncore <- 20
nthin <- 1
adapt_delta <- 0.95
stepsize <- 1
max_treedepth <- 10
options(mc.cores = ncore)

#stanmodel_arg <- file.path(repo_dir, "code", "scepticc_noepsilon.stan")
stanmodel_arg <- file.path(repo_dir, "code", "scepticc.stan")
rm(sm)
sm <- rstan::stan_model(stanmodel_arg) # compile model


setwd(results_dir)

if (run_vba) {
  # VBA: same parameter recovery as MCMC
  ee <- list()
  for (ct in 1) {
    tryCatch(fit_vba <- rstan::vb(object = sm,
                                  data   = data_list,
                                  tol_rel_obj = 0.001) # fullrank diverges 
                                  # init = my_init) # initializing with prior VBA estimates does not speed up estimation
                                  # init = get_posterior_mean(fit_vba))
             )
    tryCatch(ee[[ct]] <- extract(fit_vba))#,
    # algorithm = "meanfield",
    # # eval_elbo = 10,
    # tol_rel_obj = 0.005) # fullrank diverges
    
  }
  saveRDS(ee, file = "~/OneDrive/collected_letters/data_projects/clock_2/vba_output_no_mV_division_by_beta_unfactorized_beta1.RDS")
  
  }


# ee <- readRDS(file = "~/OneDrive/collected_letters/data_projects/clock_2/vba_output_no_mV_division_by_beta_unfactorized.RDS")
ee_vb <- extract(fit_vba)
cor(pars$alphas, colMeans(ee_vb$alpha))
cor(pars$betas, colMeans(ee_vb$beta))
cor(pars$gammas, colMeans(ee_vb$gamma))
cor(pars$epsilon_us, colMeans(ee_vb$epsilon_u))
cor(pars$epsilon_as, colMeans(ee_vb$epsilon_int))

ee1 <- list()
ee1$alpha <- colMeans(rbind(ee[[1]]$alpha, ee[[2]]$alpha, ee[[3]]$alpha, ee[[4]]$alpha))
ee1$beta <- colMeans(rbind(ee[[1]]$beta, ee[[2]]$beta, ee[[3]]$beta, ee[[4]]$beta))
ee1$gamma <- colMeans(rbind(ee[[1]]$gamma, ee[[2]]$gamma, ee[[3]]$gamma, ee[[4]]$gamma))
ee1$epsilon_u <- colMeans(rbind(ee[[1]]$epsilon_u, ee[[2]]$epsilon_u, ee[[3]]$epsilon_u, ee[[4]]$epsilon_u))
ee1$epsilon_int <- colMeans(rbind(ee[[1]]$epsilon_int, ee[[2]]$epsilon_int, ee[[3]]$epsilon_int, ee[[4]]$epsilon_int))

# using previous VBA estimates in case this run crashes
my_init <- function () {list(alpha = ee1$alpha, beta = ee1$beta, gamma = ee1$gamma, epsilon_u = ee1$epsilon_u, epsilon_int = ee1$epsilon_int)}



if(run_sampling) {
  # MCMC
  fit <- rstan::sampling(object  = sm,
                         data    = data_list,
                         #pars    = pars,
                         init    = my_init, # try initialization with vba
                         chains  = nchain,
                         iter    = niter,
                         warmup  = nwarmup,
                         thin    = nthin,
                         control = list(adapt_delta   = adapt_delta,
                                        stepsize      = stepsize,
                                        max_treedepth = max_treedepth))
  ee_s <- extract(fit)
  saveRDS(ee_s, file = "~/OneDrive/collected_letters/data_projects/clock_2/sampling_output_wide_params_no_mV_division_by_beta_unfactorized_vba_init.RDS")
  # ee_s <- readRDS("~/OneDrive/collected_letters/data_projects/clock_2/sampling_output_no_mV_division_by_beta_unfactorized.RDS")
  cor(pars$alphas, colMeans(ee_s$alpha))
  cor(pars$betas, colMeans(ee_s$beta))
  cor(pars$gammas, colMeans(ee_s$gamma))
  cor(pars$epsilon_us, colMeans(ee_s$epsilon_u))
  cor(pars$epsilon_as, colMeans(ee_s$epsilon_int))
}


# cor(pars$epsilon_us - pars$epsilon_as, colMeans(ee[[ct]]$epsilon_int))
# plot(pars$epsilon_us - pars$epsilon_as, colMeans(ee[[ct]]$epsilon_int))

# compile all recovered params for plotting

d <- pars %>% select(alphas, betas, gammas, epsilon_us, epsilon_as) %>% rename(alpha = alphas, beta = betas, gamma = gammas, epsilon_u = epsilon_us,
                                                                               epsilon_att = epsilon_as) %>%
  cbind(as_tibble(cbind(colMeans(ee_s$alpha), colMeans(ee_s$beta), colMeans(ee_s$gamma), colMeans(ee_s$epsilon_u), colMeans(ee_s$epsilon_int)))
        %>% set_names(c("alpha_s", "beta_s", "gamma_s", "epsilon_u_s", "epsilon_int_s"))) %>%
  cbind(as_tibble(cbind((ee1$alpha), (ee1$beta), (ee1$gamma), (ee1$epsilon_u), (ee1$epsilon_int))) 
        %>% set_names(c("alpha_vb", "beta_vb", "gamma_vb", "epsilon_u_vb", "epsilon_int_vb")))
library(ggExtra)
library(cowplot)
library(ggpubr)
p <- ggplot(d, aes(alpha, alpha_s)) + geom_point(alpha = .5) + stat_cor(aes(label = ..r.label..),  r.digits = 2, color = "red") + geom_abline()  + theme_half_open() 
p1 <- ggMarginal(p, type = "density")
p <- ggplot(d, aes(beta, beta_s)) + geom_point(alpha = .5) + stat_cor(aes(label = ..r.label..),  r.digits = 2, color = "red") + geom_abline()  + theme_half_open()
p2 <- ggMarginal(p, type = "density")
p <- ggplot(d, aes(gamma, gamma_s)) + geom_point(alpha = .5) + stat_cor(aes(label = ..r.label..),  r.digits = 2, color = "red") + geom_abline()  + theme_half_open()
p3 <- ggMarginal(p, type = "density")
p <- ggplot(d, aes(epsilon_u - 1/12, epsilon_u_s)) + geom_point(alpha = .5) + stat_cor(aes(label = ..r.label..),  r.digits = 2, color = "red") + geom_abline()  + theme_half_open()
p4 <- ggMarginal(p, type = "density")
p <- ggplot(d, aes(epsilon_att - 1/12, epsilon_int_s)) + geom_point(alpha = .5) + stat_cor(aes(label = ..r.label..),  r.digits = 2, color = "red") + geom_abline()  + theme_half_open()
p4a <- ggMarginal(p, type = "density")



p <- ggplot(d, aes(alpha, alpha_vb)) + geom_point(alpha = .5) + stat_cor(aes(label = ..r.label..),  r.digits = 2, color = "blue") + geom_abline()  + theme_half_open() 
p5 <- ggMarginal(p, type = "density")
p <- ggplot(d, aes(beta, beta_vb)) + geom_point(alpha = .5) + stat_cor(aes(label = ..r.label..),  r.digits = 2, color = "blue") + geom_abline()  + theme_half_open()
p6 <- ggMarginal(p, type = "density")
p <- ggplot(d, aes(gamma, gamma_vb)) + geom_point(alpha = .5) + stat_cor(aes(label = ..r.label..),  r.digits = 2, color = "blue") + geom_abline()  + theme_half_open()
p7 <- ggMarginal(p, type = "density")
p <- ggplot(d, aes(epsilon_u - 1/12, epsilon_u_vb)) + geom_point(alpha = .5) + stat_cor(aes(label = ..r.label..),  r.digits = 2, color = "blue") + geom_abline()  + theme_half_open()
p8 <- ggMarginal(p, type = "density")
p <- ggplot(d, aes(epsilon_att - 1/12, epsilon_int_vb)) + geom_point(alpha = .5) + stat_cor(aes(label = ..r.label..),  r.digits = 2, color = "blue") + geom_abline()  + theme_half_open()
p8a <- ggMarginal(p, type = "density")

p <- ggplot(d, aes(alpha_s, alpha_vb)) + geom_point(alpha = .5) + stat_cor(aes(label = ..r.label..),  r.digits = 2, color = "purple") + geom_abline()  + theme_half_open() 
p9 <- ggMarginal(p, type = "density")
p <- ggplot(d, aes(beta_s, beta_vb)) + geom_point(alpha = .5) + stat_cor(aes(label = ..r.label..),  r.digits = 2, color = "purple") + geom_abline()  + theme_half_open()
p10 <- ggMarginal(p, type = "density")
p <- ggplot(d, aes(gamma_s, gamma_vb)) + geom_point(alpha = .5) + stat_cor(aes(label = ..r.label..),  r.digits = 2, color = "purple") + geom_abline()  + theme_half_open()
p11 <- ggMarginal(p, type = "density")
p <- ggplot(d, aes(epsilon_u_s, epsilon_u_vb)) + geom_point(alpha = .5) + stat_cor(aes(label = ..r.label..),  r.digits = 2, color = "purple") + geom_abline()  + theme_half_open()
p12 <- ggMarginal(p, type = "density")
p <- ggplot(d, aes(epsilon_int_s, epsilon_int_vb)) + geom_point(alpha = .5) + stat_cor(aes(label = ..r.label..),  r.digits = 2, color = "purple") + geom_abline()  + theme_half_open()
p12a <- ggMarginal(p, type = "density")


setwd(results_dir)
pdf("wide_pars_resovery_vb_mcmc4000_noVm_division_by_beta_unfactorized.pdf", 12, 7)
cowplot::plot_grid(as_grob(p1), as_grob(p2), as_grob(p3), as_grob(p4), as_grob(p4a),
                   as_grob(p5), as_grob(p6), as_grob(p7), as_grob(p8), as_grob(p8a),
                   as_grob(p9), as_grob(p10), as_grob(p11), as_grob(p12), as_grob(p12a), nrow = 3) + theme_half_open()
dev.off()

cormat_s = cor(d %>% select(contains("_s")))


cormat_vb = cor(d %>% select(contains("_vb")))

pdf("wide_pars_cormat_vb_noVm_division_by_beta_unfactorized.pdf", 4, 4)
corrplot::corrplot(cormat_vb, method = 'number', type = "lower", diag = F)
dev.off()

pdf("wide_pars_cormat_mcmc4000_noVm_division_by_beta_unfactorized.pdf", 4, 4)
corrplot::corrplot(cormat_s, method = 'number', type = "lower", diag = F)
dev.off()



parVals <- rstan::extract(fit, permuted = TRUE)

# Quick estimation with 1 iteration to get computed quantities for debugging
fit_debug <- rstan::sampling(
  object = sm, data = data_list,
  chains = 1, iter = 2, warmup = 0, thin = nthin,
  control = list(
    adapt_delta = adapt_delta,
    stepsize = stepsize,
    max_treedepth = max_treedepth
  )
)

ee <- extract(fit_debug)

# test position of positions
abs(2 * pi - (ee$p_pvec[1, 60] + ee$p_pvec[1, 1])) < 1e-5
abs(2 * pi - (ee$p_pvec[1, 60] + ee$p_pvec[1, 1])) < 1e-5
all.equal(diff(ee$p_breaks[1, ]))

data_list$segment_min[1, 1]
data_list$segment_max[1, 1]
ee$p_S[1, 1, 1, ]
ee$p_breaks[1, ]

fit_vb <- rstan::vb(
  object = sm,
  data = data_list
)



# examine basis
basis <- extract(fit_debug)$p_Phi[1, , ]
bdf <- reshape2::melt(basis, varnames = c("basis", "point"))
ggplot(bdf, aes(x = point, y = value, color = factor(basis))) + geom_line() # looks good

# compare against generation in R
n_basis <- 12
d_theta <- (2 * pi) / n_basis
basis_sd <- rep(0.3, 12)
weights_0 <- rep(0, 12)

# since circle wraps, we want to avoid adding a redundant basis function at 0 vs. at 2*pi
loc <- seq(0, 2*pi - d_theta, by=d_theta)
vset <- lapply(seq_along(loc), function(ii) {
  vm_bf$new(n_points = 120, center=loc[ii], width_sd=basis_sd[ii], weight=weights_0[ii])
})

bf_set <- rbf_set$new(elements = vset)
basis_r <- t(bf_set$get_basis())

summary(basis - basis_r)
all.equal(basis, basis_r) # TRUE -- YAY!

# examine different approaches to binning data
mv <- 2*pi # maximum value
np <- 120 # number of points
d_theta <- 2*pi/np
pvec <- seq(0, 2 * pi - d_theta, by = d_theta) # positions vector

# 0-359 degrees in radians
## dd <- seq(0, 2 * pi - (2 * pi / 360), length.out = 360)
## dd * 180 / pi

# allow wrap
dd <- seq(0, 2 * pi - 1e-3, length.out = 361)

# look more finely in degrees
dd <- seq(0, 360, by = 0.25)

pvec <- seq(0, 360, length.out=121)
table(.bincode(dd, pvec, include.lowest = FALSE)) # function used by cut() in R

# look at approach of using rounding to bin
mv <- 360
vv <- round((np - 1) * dd / mv) + 1
df <- data.frame(dd, vv)
df %>%
  group_by(vv) %>%
  summarise(mi = min(dd), ma = max(dd), di = max(dd) - min(dd)) %>%
  View()

# c++ friendly -- but not as accurate as .bincode
vv <- round((np - 1) * dd / mv) + 1
# vv2 <- round(dd / (mv - pi / 180) * (np - 1)) + 1

# subtract off overlap
# vv <- round((np - 1) * dd / (mv - (pi / 180))) + 1

rr <- as.integer(Hmisc::cut2(dd, cuts = pvec))
rr2 <- as.integer(cut(dd, 120))
rr3 <- choices_to_pbins(dd, 120)
cbind(vv, rr, rr2, rr3)
cor(vv, rr)

# conclusion: .bincode gives the greatest control by enforcing precise bin boundaries

# Eligibility of each basis function based on location of choice -- looks good!
EB <- extract(fit_debug)$p_EB[1, , ]
edf <- reshape2::melt(EB, varnames=c("point", "basis"))
ggplot(edf, aes(x = point, y = value, color = factor(basis))) + geom_line()

# Resolve complexity of defining evaluation points *centers* in pvec versus bins that span those centers
# Define bins based on evaluation point centers
# Also check .bincode approach against stan-based port/implementation
P <- 120
d_theta = 2*pi/P

pvec <- seq(1/2*d_theta, 2 * pi - 1/2*d_theta, length.out = P)

# we have p+1 cuts, with the last at 2pi
cuts <- c(pvec - 1 / 2 * d_theta, 2 * pi)

dd <- seq(0, 2*pi, length.out=240)
bb <- .bincode(dd, cuts, include.lowest = TRUE)

table(bb)
cbind(dd*180/pi, bb)

choices_binned <- extract(fit_debug)$p_choices[1, , ]

ee <- extract(fit_debug)
ee$p_pvec[1,]
all.equal(ee$p_breaks[1, ], cuts) # TRUE -- GOOD!
which(abs(ee$p_breaks[1, ] - cuts) > .05)

summary(diff(diff(cuts)))
summary(diff(diff(ee$p_breaks[1, ])))

# match
cuts[115:121]
ee$p_breaks[1,115:121]

# compare stan against local cutting
loc <- apply(choices, c(1, 2), .bincode, breaks = cuts, include.lowest = TRUE)
attr(loc, "dimnames") <- NULL
loc[1:10, 1:10]
choices_binned[1:10, 1:10]

# excellent
all.equal(loc, choices_binned) # TRUE

# compare against older cut2 approach above
str(rc)
str(choices_binned)
summary(choices_binned - rc)
cor(as.vector(rc), as.vector(choices_binned))
choices_binned[1:10, 1:10]
rc[1:10, 1:10]
