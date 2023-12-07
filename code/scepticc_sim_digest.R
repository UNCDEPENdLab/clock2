# analyze result of SCEPTICC simulations
# pick the best contingencies

# first: rsync clock2/simulations from CRC
# general goal: pick 100 least-correlated contingencies, stress-test them with e.g. 1000 seeds across parameters, select e.g. top 10 
# in the future, we can recycle the top 100 for different numbers of trials

library(tidyverse)
library(data.table)
repo_dir <- "~/code/clock2"
sim_dir <- "~/code/clock2/simulations/from_crc"
out_dir <- "~/code/clock2/simulations"
setwd(sim_dir)

system(
  paste0("rsync -azP --include '*file*' --exclude '*'  ayd1@htc.crc.pitt.edu:/ihome/adombrovski/ayd1/code/clock2/simulations/"))

files <- list.files(pattern = "*file*") # this is for all contingencies with 100 seeds

library(furrr)
future::plan(multisession)
df <- files %>% future_map_dfr(fread)
dimensions <- names(df %>% select(-r, -seed, -block_length, -low_avg))

# remove NA seeds/r
df <- df %>% filter(!is.na(seed))
# check missingness pattern: no missing cells
# t <- as_tibble(as.data.frame(table(df$iteration, df$seed))) %>% setNames(c("iteration", "seed", "freq"))
# setwd(out_dir)
# pdf("missing_iteration_seeds.pdf", height = 20, width = 30)
# ggplot(t, aes(as.numeric(seed), freq, color = iteration)) + geom_line()
# dev.off()

# low_avg = 10 and timeout = 2 generally work best, stick with these
# sdf <- df %>% summarise(.by = all_of(dimensions), r100 = mean(r)) %>% arrange(alpha, gamma, beta, epsilon_u, block_length, timeout, low_avg, iteration)
sdf <- df %>% summarise(.by = all_of(dimensions), r100 = mean(r), sd100 = sd(r)) %>% arrange(alpha, gamma, beta, epsilon_u, iteration)
sdf <- df %>% summarise(.by = "iteration", r100 = mean(r), sd100 = sd(r), maxr = max(r), minr = min(r), spanr = maxr-minr)
top10 <- sdf %>% filter(abs(r100) < .1 & sd100 < .1 & abs(maxr) < 0.2935 & abs(minr) < 0.28) %>% select(iteration) %>% as.array()


data.table::fwrite(sdf, file = paste0("crc_sim_results_subjects", length(unique(sdf$iteration)), "_iterations.csv"))
# best low_avg, iterations
iterations <- unique(sdf$iteration)

idf <- sdf %>% summarise(.by = c(low_avg, iteration, epsilon_u), r100 = mean(r100), sd100 = mean(sd100)) %>% filter(epsilon_u == .9) %>% arrange(abs(r100), sd100) %>% top_n(-100)
top100 <- unique(idf$iteration)
print(paste(top100, collapse=", "), width = 12) # to paste in plot_winning_contingencies.R
# just to inspect
idf %>% select(iteration, r100, sd100)
# effect of parameters
# all
df %>% summarise(.by = c(gamma, alpha, beta, epsilon_u), r = mean(r)) %>% arrange(r)
# top 100
df %>% filter(iteration %in% top100) %>% summarise(.by = c(gamma, alpha, beta, epsilon_u), r = mean(r)) %>% arrange(r)

##############
# plot top 100
# and top 10
##############
plot_df <- df %>% filter()
ggplot(plot_df, aes(as.factor(iteration), r100, color = gamma, size = beta, alpha = alpha)) + geom_point() + facet_grid(~epsilon_u)
# find seeds with some high correlations
sdf %>% filter(iteration %in% top100[1:10] & r100>0.4) %>% select(iteration) %>% unique()

m1 <- aov(r ~ alpha*gamma*beta*epsilon_u * iteration, df)
library(lsr)
etaSquared(m1)

############
# top 100 stress-tested with 1000 iterations each, 89 completed/11 crashed
# files77 <- list.files(pattern = "*true_top*") # top 100 with 1000 seeds each
# df77 <- files77 %>% map_dfr(fread)
# dimensions <- names(df77 %>% select(-r, -seed))

# low_avg = 10 and timeout = 2 generally work best, stick with these
sdf77 <- df77 %>% summarise(.by = all_of(dimensions), r1000 = mean(r)) %>% arrange(alpha, gamma, beta, epsilon_u, block_length, timeout, low_avg, iteration)
write_csv2(sdf77, file = paste0("crc_sim_results_low_avg10_timeout2_top_", 77, "_iterations.csv"))
# best low_avg, iterations
iterations <- unique(sdf77$iteration)

idf77 <- sdf77 %>% summarise(.by = c(low_avg, iteration, timeout, epsilon_u), r1000 = mean(r1000)) %>% filter(epsilon_u == .9) %>% arrange(r1000)

# compare 100 vs. 1000 seeds
comp_df <- inner_join(idf77, idf)
comp_df %>% select(iteration, r100, r1000) %>% arrange(r1000)
comp_df$diff <- comp_df$r100 - comp_df$r1000
ggplot(comp_df) + geom_point(aes(iteration, diff), color = "red")# + ylim(-0.00001, 0.00001)
hist(comp_df$diff)
# final contingencies (screened through plotting)
final_4 <- idf77 %>% select(iteration, r1000) %>% filter(iteration %in% c(2455, 6820, 8093, 5643)) # 7319 is not in this set
write_csv2(final_4, file = paste0("final_4_seeds_selected.csv"))

# 
# # filtered dataset: 
# comp_sdf <- inner_join(sdf, sdf77) %>% arrange(r1000)
# 
# # bottom line, best parameters:
# # block_length = 10
# # avg_low = 10
# # timeout = 2
# # ntrials = 300
# # best seeds:
# # print(comp_df$iteration[1:10])
# # 8014 5629 5888 6791 8395 6623 4753 2455 8396 6604
# 
# top10_df <- comp_df %>% arrange(r1000) 
# 
# top1 <- sdf %>% filter(low_avg == 15 & iteration == 4935 & timeout == 2)
# ggplot(top1, aes(as.factor(epsilon_u), r, color = gamma, pch = as.factor(alpha), size = as.factor(beta))) + geom_point()
# top1 %>% summarise(.by = epsilon_u, mean_r = mean(r))
# # pretty good!
# 
# top2 <- sdf %>% filter(low_avg == 10 & iteration == 2124 & timeout == 2)
# ggplot(top2, aes(as.factor(epsilon_u), r, color = gamma, pch = as.factor(alpha), size = as.factor(beta))) + geom_point()
# top2 %>% summarise(.by = epsilon_u, mean_r = mean(r))
# # worst 2 rs > .5
# 
# top3 <- sdf %>% filter(low_avg == 10 & iteration == 4998 & timeout == 2)
# ggplot(top3, aes(as.factor(epsilon_u), r, color = gamma, pch = as.factor(alpha), size = as.factor(beta))) + geom_point()
# top3 %>% summarise(.by = epsilon_u, mean_r = mean(r))
# # good!
# 
# top4 <- sdf %>% filter(low_avg == 10 & iteration == 2333 & timeout == 2)
# ggplot(top4, aes(as.factor(epsilon_u), r, color = gamma, pch = as.factor(alpha), size = as.factor(beta))) + geom_point()
# top4 %>% summarise(.by = epsilon_u, mean_r = mean(r))
# # also good, but bigger difference between hi/lo epsilon_u
# 
# 
# 
# top_50_iterations <- df %>% summarise(.by = c(low_avg, iteration), mean_r = mean(r)) %>% arrange(mean_r) %>% top_n(-50) %>% select(iteration)
# 
# # best low_avg, iterations for worst-case scenario (low beta, gamma > .1, epsilon_u = .99)
# wdf <- df %>% filter(beta == 1 & gamma > .1 & epsilon_u == .99) %>% 
#   summarise(.by = c(low_avg, iteration, gamma), mean_r_worst = mean(r)) %>% arrange(mean_r_worst) %>% top_n(-20) %>% inner_join(idf)
# write_csv2(wdf, file = "winning_6_lowavg_seed_robust.csv")
# 
# 
# 
# winners <- wdf %>% select(low_avg, iteration)
# write_csv2(winners, file = "winning_6_lowavg_seed.csv")
# 
# # sdf <- big_df %>% summarise(.by = c(iteration, epsilon_u, alpha, beta, gamma), mean_r = mean(r))
# sdf <- big_df %>% summarise(.by = c(iteration, epsilon_u), mean_r = mean(r))
# 
# pdf("simulations_plot_worst.pdf", height = 20, width = 20)
# ggplot(df, aes(as.factor(epsilon_u), r, color = as.factor(iteration))) + geom_boxplot() #+ facet_grid(alpha + gamma ~ beta)
# dev.off()
# 
# pdf("sim_90_iterations.pdf", height = 10, width = 50)
# ggplot(sdf, aes(as.factor(iteration), mean_r, color = as.factor(epsilon_u), lty = as.factor(low_avg))) + geom_boxplot() #+ theme(legend.position = "none")
# dev.off()
# ggplot(df, aes(low_avg, r, color = as.factor(iteration))) + geom_smooth(method = "loess")
# 
# pdf("simulations_jitter_worst.pdf", height = 20, width = 20)
# ggplot(df, aes(low_avg, r, color = as.factor(iteration), pch = as.factor(block_length))) + geom_jitter() + 
#   facet_grid(alpha + gamma ~ epsilon_u + beta) + geom_hline(yintercept = 0, size = 1) + geom_hline(yintercept = 0.25, size = .5) + geom_hline(yintercept = 0.5, size = .25)
# dev.off()
# 
# ggplot(df, aes(as.factor(low_avg), r, color = as.factor(iteration), lty = as.factor(block_length))) + geom_violin(draw_quantiles = .5)
# 
# 
# ggplot(df, aes(as.factor(low_avg), r, color = as.factor(block_length))) + geom_boxplot() #+ facet_wrap(~drift)
# ggplot(df, aes(as.factor(block_length), r, color = as.factor(alpha))) + geom_boxplot() + facet_wrap(~as.factor(gamma))
# ggplot(df, aes(as.factor(block_length), r, color = as.factor(alpha))) + geom_boxplot() + facet_wrap(~as.factor(beta))
# ggplot(df, aes(as.factor(bump_prom), r, color = as.factor(alpha))) + geom_boxplot() #+ facet_wrap(~drift)
# 
# 
# ggplot(df, aes(low_avg, r)) + geom_jitter()  
# 
# # focus on problematic behaviors, worst-case scenario: block_length of 10-15 and low_avg 10-20 is still best
# df <- big_df %>% filter(epsilon_u > 0.8 & beta < 20 & alpha <.2)
# 
# m0 <- lme4::lmer(r ~  (1|iteration), df)
# summary(m0)
# 
# 
# m1 <- lme4::lmer(r ~ (as.factor(low_avg) + epsilon_u + beta + gamma)^2 + (1|iteration), df)
# car::Anova(m1, '3')
# summary(m1)
# 
# lm1 <- lm(r ~ (epsilon_u + scale(alpha) + scale(beta) + scale(gamma) + as.factor(iteration))^2, df)
# anova(lm1)
# summary(lm1)
# 
# # look at only the highest epsilon_u
# lm2 <- lm(r ~ as.factor(iteration), df)
# anova(lm2)
# summary(lm2)
# 
# lm3 <- lm(r ~ as.factor(iteration), df %>% filter(epsilon_u > .9))
# anova(lm3)
# summary(lm3)
# 
# 
# lm3 <- lm(r ~ (as.factor(block_length) + as.factor(low_avg) + as.factor(iteration)), df)
# anova(lm3)
# summary(lm3)
# 
