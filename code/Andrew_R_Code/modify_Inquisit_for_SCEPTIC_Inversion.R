# 2024-04-29 AndyP
# Modify Inquisit data for SCEPTIC inversion

library(tidyverse)
library(BAMBI)

data_dir <- '/Users/andypapale/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Documents - DNPLskinner/skinner'

df0 <- read_csv(paste0(data_dir,'/data/prolific/clock_v2.1_pilot/prolific_01-09/raw/2024-01-09-Raw-noValueVector.csv'))
df0 <- df0 %>% filter(subject != '6449707327ff66156c264c6f' & subject!='6464b45cf8e8a0f06fe011b9') # remove subjects with corrupted data
df0 <- df0 %>% filter(trialcode == 'feedback')

nsubj <- length(unique(df0$subject))

# get outcomes nsubj x T  (points earned = inc_rg in Inquisit, number of subjects by trial matrix)
outcomes <- df0 %>% 
  select(subject, trial, inc_rg) %>% 
  pivot_wider(id_cols = "subject", values_from = "inc_rg", names_from = "trial") %>% 
  select(!subject) %>%
  data.matrix()

# get choices_rad nsubj x T (choice of subject (rad 0 to 2pi) = pos_shifted in Inquisit (deg), number of subjects by trial matrix)
choices <- df0 %>% 
  mutate(pos_shifted_rad = BAMBI::zero_to_2pi(pos_shifted*pi/180)) %>% 
  select(subject, trial, pos_shifted_rad) %>% 
  pivot_wider(id_cols = "subject", values_from = "pos_shifted_rad", names_from = "trial") %>% 
  select(!subject) %>%
  data.matrix()

# Note Inquisit codes timeout trials as no erasure, and MNH codes timeout trials as erasure, is this a problem?
trial_types <- df0 %>% 
  mutate(trial_types = case_when(trial_type == "no erasure" ~ 1, trial_type == "erasure" ~ 2, trial_type == "attention" ~ 3)) %>%
  select(subject, trial, trial_types) %>% 
  pivot_wider(id_cols = "subject", values_from = "trial_types", names_from = "trial") %>% 
  select(!subject) %>%
  data.matrix()

segment_shown <- df0 %>%
  mutate(segment_shown = case_when(trial_type == "no erasure" ~ 0, trial_type =="erasure" ~ 1, trial_type=="attention" ~ 1)) %>%
  select(subject,trial,segment_shown) %>%
  pivot_wider(id_cols = "subject", values_from = "segment_shown", names_from = "trial") %>% 
  select(!subject) %>%
  data.matrix()


segment_min <- df0 %>% 
  group_by(subject,trial) %>% mutate(segment_min_deg = min(stim_left_deg,stim_right_deg)) %>% ungroup() %>%
  mutate(segment_min = case_when(trial_type != "no erasure" ~ BAMBI::zero_to_2pi(segment_min_deg*pi/180), trial_type == "no erasure" ~ -99)) %>%
  select(subject,trial,segment_min) %>%
  pivot_wider(id_cols = "subject", values_from = "segment_min", names_from = "trial") %>% 
  select(!subject) %>%
  data.matrix()

segment_max <- df0 %>% 
  group_by(subject,trial) %>% mutate(segment_max_deg = max(stim_left_deg,stim_right_deg)) %>% ungroup() %>%
  mutate(segment_max = case_when(trial_type != "no erasure" ~ BAMBI::zero_to_2pi(segment_max_deg*pi/180), trial_type == "no erasure" ~ -99)) %>%
  select(subject,trial,segment_max) %>%
  pivot_wider(id_cols = "subject", values_from = "segment_max", names_from = "trial") %>% 
  select(!subject) %>%
  data.matrix()

# zero_to_2pi wraps 360 to 0, making 360 the 'min' and 30 the 'max', correct this
segment_max[round(segment_max - segment_min,7)==-0.5235988] = 30*pi/180
segment_min[round(segment_max - segment_min,7)==-0.5235988] = 0

data_list <- list(
  B = 16, # 12 basis functions
  P = 90, # 90 evaluation points for value
  N = nsubj,
  T = 240, # fmri
  Tsubj = rep(240, nsubj),
  sb = 0.3, # SD of basis functions in radians
  sg = 0.3, # SD of generalization function in radians
  outcomes = outcomes,
  choices_rad = choices,
  trial_types = trial_types,
  segment_shown = segment_shown,
  segment_min = segment_min,
  segment_max = segment_max
)
