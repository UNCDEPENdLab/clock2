# R script for handling cluster queuing for clock2 simulations
library(tidyverse)
library(data.table)
#library(BAMBI)
if (sum(stringr::str_detect(Sys.info(), "Alex|alexdombrovski"))>1) {
  basedir <- "~/code/clock2/code/"
  output_dir <- "~/code/clock2/simulations"
  sbatch_dir <- "~/code/clock2/code/sbatch/"
} else if (sum(stringr::str_detect(Sys.info(), "andypapale"))>1)  {
  basedir <- "~/clock2/code/"
  output_dir <- "~/clock2/simulations"
  sbatch_dir <- "~/clock2/code/sbatch/"
} else {
  basedir <- "~/code/clock2/code/"
  output_dir <- "~/code/clock2/simulations"
  sbatch_dir <- "~/code/clock2/code/sbatch/"
}

test <- F
test_on_mac <- F
animate <- F
subject_params <- T

if(subject_params) {
  top100 = c(4931, 5608, 754, 5781, 7480, 6508, 5340, 2637, 5871, 8056, 5094, 6220, 720, 1269, 5903, 8608, 4505, 567, 5173, 1135, 7174, 1506, 7323, 7347, 283, 1291, 5966, 8592, 5655, 1041, 7315, 124, 2534, 6821, 6658, 4205, 5815, 2619, 430, 7527, 5299, 5438, 1337, 8765, 4025, 152, 6520, 5213, 6729, 4938, 1250, 4563, 2488, 6644, 5288, 1070, 1803, 8706, 136, 7530, 155, 5646, 727, 6227, 6108, 1639, 5307, 1764, 7933, 6583, 6578, 8021, 6896, 1752, 6176, 5877, 531, 7595, 6252, 6342, 4222, 5969, 6623, 7379, 1246, 7432, 6783, 5688, 8448, 4175, 5413, 1155, 2825, 7564, 6749, 1464, 5436, 294, 4998, 868)
iterations <- top100
# get subjects' parameters
setwd(output_dir)
sub_df <- fread("mmclock_fmri_decay_mfx_sceptic_global_statistics.csv") %>%
  select(alpha_transformed, gamma_transformed, beta_transformed) %>% rename(alpha = alpha_transformed, gamma = gamma_transformed, beta = beta_transformed)


grid <- expand.grid(epsilon_u = c(0.3, 0.9), # 0.0833 is at chance, low correlation -- not worth testing
                  block_length = c(10), # block length > 15 had higher correlations, not worth testing
                  low_avg = c(10),
                  # iteration = c(1:1000),
                  iteration = iterations,
                  #drift = c(1, 2, 4), bump_prom = c(8, 10, 15),
                  seed = 1:100) %>% merge(sub_df)
idf_list <- grid %>% group_split(iteration, seed)
} else {
rob_grid <- expand.grid(alpha = c(0.1, 0.2, 0.5), gamma = c(0.1, 0.5, 0.9),                 # model params
                          beta = c(1, 5), # at very high betas, h and u are decorrelated, no need to test
                          epsilon_u = c(0.3, 0.9), # 0.0833 is at chance, low correlation -- not worth testing
                          block_length = c(10), # block length > 15 had higher correlations, not worth testing
                          low_avg = c(10),
                          iteration = 1:10000, # top100,
                        timeout = c(2),
                          #drift = c(1, 2, 4), bump_prom = c(8, 10, 15),
                          seed = 1:100)
idf_list <- rob_grid %>% group_split(iteration)
}

# if (test && test_on_mac) {
#   basedir <- "~/OneDrive - University of Pittsburgh/Momentum_EMA/eeg_data_t_split/"
#   output_dir <- "~/OneDrive - University of Pittsburgh/Momentum_EMA/"
# }
setwd(output_dir)
silent <- F

# set up simulation grid, write files


niterations <- length(idf_list)

#for (f in 1:2) {data.table::fwrite(idf_list[[f]], file = paste0("subjects_grid_", f, ".csv"))}

for (f in 1:niterations) {data.table::fwrite(idf_list[[f]], file = paste0("subjects_grid_", f, ".csv"))}

for (f in 1:niterations) {
# for (f in 1:2) {

  if (!test) {
    system(
      paste0(
        "cd ", sbatch_dir, "; ",
        "sbatch --time=23:59:59 --mem=4g",
        " --export=sourcefilestart=", f,
        " sbatch_clock2_sim.bash"
      )
    )
    #write compute level to temporary file
  }
  if (!silent) {
    cat(
      paste0(
        "cd ", sbatch_dir, "; ",
        "sbatch --time=23:59:59 --mem=4g",
        " --export=sourcefilestart=", f,
        " sbatch_clock2_sim.bash\n"
      )
    )
    }
  }