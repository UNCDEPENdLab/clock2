





df1 <- read_csv('/Users/andypapale/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Documents - DNPLskinner/skinner/data/prolific/clock_v2.1_pilot/prolific_01-09/raw/papalea_prosper_clock_2_1_1_seed_152_fmri_raw_2401091855.csv') %>% filter(trialcode != 'getValueVector')
df2 <- read_csv('/Users/andypapale/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Documents - DNPLskinner/skinner/data/prolific/clock_v2.1_pilot/prolific_01-09/raw/papalea_prosper_clock_2_1_1_seed_868_fmri_raw_2401091851.csv') %>% filter(trialcode != 'getValueVector')
df3 <- read_csv('/Users/andypapale/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Documents - DNPLskinner/skinner/data/prolific/clock_v2.1_pilot/prolific_01-09/raw/papalea_prosper_clock_2_1_1_seed_1464_fmri_raw_2401091853.csv') %>% filter(trialcode != 'getValueVector')
df4 <- read_csv('/Users/andypapale/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Documents - DNPLskinner/skinner/data/prolific/clock_v2.1_pilot/prolific_01-09/raw/papalea_prosper_clock_2_1_1_seed_1752_fmri_raw_2401091852.csv') %>% filter(trialcode != 'getValueVector')
df5 <- read_csv('/Users/andypapale/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Documents - DNPLskinner/skinner/data/prolific/clock_v2.1_pilot/prolific_01-09/raw/papalea_prosper_clock_2_1_1_seed_2534_fmri_raw_2401091854.csv') %>% filter(trialcode != 'getValueVector')
df6 <- read_csv('/Users/andypapale/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Documents - DNPLskinner/skinner/data/prolific/clock_v2.1_pilot/prolific_01-09/raw/papalea_prosper_clock_2_1_1_seed_4938_fmri_raw_2401091851.csv') %>% filter(trialcode != 'getValueVector')
df7 <- read_csv('/Users/andypapale/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Documents - DNPLskinner/skinner/data/prolific/clock_v2.1_pilot/prolific_01-09/raw/papalea_prosper_clock_2_1_1_seed_5094_fmri_raw_2401091850.csv') %>% filter(trialcode != 'getValueVector')
df8 <- read_csv('/Users/andypapale/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Documents - DNPLskinner/skinner/data/prolific/clock_v2.1_pilot/prolific_01-09/raw/papalea_prosper_clock_2_1_1_seed_5173_fmri_raw_2401091854.csv') %>% filter(trialcode != 'getValueVector')
df9 <- read_csv('/Users/andypapale/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Documents - DNPLskinner/skinner/data/prolific/clock_v2.1_pilot/prolific_01-09/raw/papalea_prosper_clock_2_1_1_seed_5815_fmri_raw_2401091851.csv') %>% filter(trialcode != 'getValueVector')
df10 <- read_csv('/Users/andypapale/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Documents - DNPLskinner/skinner/data/prolific/clock_v2.1_pilot/prolific_01-09/raw/papalea_prosper_clock_2_1_1_seed_6520_fmri_raw_2401091851.csv') %>% filter(trialcode != 'getValueVector')

df <- rbind(df1,df2,df3,df4,df5,df6,df7,df8,df9,df10)
setwd('/Users/andypapale/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Documents - DNPLskinner/skinner/data/prolific/clock_v2.1_pilot/prolific_01-09/raw/')
write.csv(df,'2023-01-09-Raw-noValueVector.csv')
