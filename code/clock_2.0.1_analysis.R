# 2023-05-03 AndyP
# clock 2.0.1 analysis

base_dir <- '/Users/andypapale/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Documents - DNPLskinner' # on skinner

design_file <- '~/clock2/2022-04-25-DesignFile.csv' # in Michael's clock2 repo

df <- read_csv(paste0(base_dir,'/skinner/data/prolific/clock_v2.1_pilot/In_lab_pilot_05-2023/papalea_prosper_eeg_clock_v3_raw_2305011503.csv'))

df <- df %>% filter(trialcode=='dispFeedback_noU' | trialcode=='dispFeedback_U') # the score is calculated here, the variable of interest is rt_shifted

# latency is the latency of the trial (for feedback always 1000, to check latency of choice/RT referenced to starting position filter by trialcode=='experiment_U' | trialcode=='experiment_noU')
# rt_shifted is the position around the clock chosen
# startPos is the starting position range [0 100] for one complete revolution, 12 o'clock is 75=100, 9 o'clock is 50, 6 o'clock is 25, 3 o'clock is 0
# meta_trialCount is trial [1... total number of trials ~ 300]
# trialCount is trial [1 ... trials within block].
# blockCount is block [1-8] blocks are 32-43 trials and contain 1 attentional control block and one uncertainty manipulation block
# ntrials is the number of trials per block (32-43)
# attentionalControl is either "wind" or "cloud" this is the attentional control that subjects learn has no effect
# localUncertainty is either "wind" or "cloud" this is the local uncertainty manipulation that subjects learn changes the number of mushrooms
# uncertaintyBlock is either 0,1,2 within each block, corresponding to the attentionalControl or localUncertainty manipulation
# chooseFog1 is 1 if "cloud" was chosen, something seems to be wrong with this counter 2022-05-03 AndyP, possibly not resetting to 0 after the block is complete?
# chooseWind2 is 1 if "wind" was chosen
# fogPos1 determines the position of the "cloud" stimulus, either [0,1,2,5,10,15,20,25,30,45,50] in 30 degree rotations counter clockwise starting at 12 o'clock for 0
# windPos1 determines the position of the "wind" stimulus, either [0,1,2,5,10,15,20,25,30,45,50] in 30 degree rotations counter clockwise starting at 12 o'clock for 0
# totalEarnings is cumulative earnings in units of 0.5c / mushroom, multiply by 2 to get # mushrooms
# Earnings is earnings on the current trial
# unc_att_start1 is the starting trial of the "cloud" stimulus on each block
# unc_att_start2 is the starting trial fo the "wind" stimulus on each block
# choose_unc_att1 counts the number of times the "cloud" stimulus is chosen (should not exceed 2)
# choose_unc_att2 counts the number of times the "wind" stimulus is chosen (should not exceed 2)
# inc is the number of points on that trial for the given rt_shifted
# rng is the random number generated to determine if the trial is an 'omission' (1 mushroom given).  We should have a fixed probability (freq) of 0.7 for all RT's and trials
# 