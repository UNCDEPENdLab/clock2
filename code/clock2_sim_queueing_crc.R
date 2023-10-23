# R script for handling cluster queuing for clock2 simulations
library(tidyverse)
basedir <- "~/code/clock2/code/"
output_dir <- "~/code/clock2/simulations"
sbatch_dir <- "~/code/clock2/code/sbatch/"
test <- F
test_on_mac <- F
# if (test && test_on_mac) {
#   basedir <- "~/OneDrive - University of Pittsburgh/Momentum_EMA/eeg_data_t_split/"
#   output_dir <- "~/OneDrive - University of Pittsburgh/Momentum_EMA/"
# }
setwd(output_dir)
silent <- F

# set up simulation grid, write files
rob_grid <- expand.grid(alpha = c(0.2, 0.5), gamma = c(0.1, 0.5, 0.9),                 # model params
                        beta = c(1, 5), # at very high betas, h and u are decorrelated, no need to test
                        epsilon_u = c(0.3, 0.9), # 0.0833 is at chance, low correlation -- not worth testing
                        block_length = c(10), # block length > 15 had higher correlations, not worth testing
                        low_avg = c(10, 20),
                        iteration = c(223:1222),
                        #drift = c(1, 2, 4), bump_prom = c(8, 10, 15),
                        seed = 1:100)
idf_list <- rob_grid %>% group_split(iteration)
niterations <- length(idf_list)
for (f in 1:niterations) {data.table::fwrite(idf_list[[f]], file = paste0("grid_", f, ".csv"))}


#for (f in 1:niterations) {
  for (f in 1:2) {

  if (!test) {
    system(
      paste0(
        "cd ", sbatch_dir, "; ",
        "sbatch --time=23:00:00 --mem=1g",
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
        "sbatch --time=23:00:00 --mem=1g",
        " --export=sourcefilestart=", f,
        " sbatch_clock2_sim.bash\n"
      )
    )
  }
}
