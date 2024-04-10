# 2024-01-09 AndyP
# Testing value distributions from Prolific Pilot

#df1 <- read_csv('/Users/andypapale/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Documents - DNPLskinner/skinner/data/prolific/clock_v2.1_pilot/prolific_01-09/raw/papalea_prosper_clock_2_1_1_seed_152_fmri_raw_2401091855.csv')
#df1 <- read_csv('/Users/andypapale/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Documents - DNPLskinner/skinner/data/prolific/clock_v2.1_pilot/prolific_01-09/raw/papalea_prosper_clock_2_1_1_seed_868_fmri_raw_2401091851.csv')
#df1 <- read_csv('/Users/andypapale/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Documents - DNPLskinner/skinner/data/prolific/clock_v2.1_pilot/prolific_01-09/raw/papalea_prosper_clock_2_1_1_seed_1464_fmri_raw_2401091853.csv')
#df1 <- read_csv('/Users/andypapale/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Documents - DNPLskinner/skinner/data/prolific/clock_v2.1_pilot/prolific_01-09/raw/papalea_prosper_clock_2_1_1_seed_1752_fmri_raw_2401091852.csv')
#df1 <- read_csv('/Users/andypapale/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Documents - DNPLskinner/skinner/data/prolific/clock_v2.1_pilot/prolific_01-09/raw/papalea_prosper_clock_2_1_1_seed_2534_fmri_raw_2401091854.csv')
#df1 <- read_csv('/Users/andypapale/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Documents - DNPLskinner/skinner/data/prolific/clock_v2.1_pilot/prolific_01-09/raw/papalea_prosper_clock_2_1_1_seed_4938_fmri_raw_2401091851.csv')
#df1 <- read_csv('/Users/andypapale/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Documents - DNPLskinner/skinner/data/prolific/clock_v2.1_pilot/prolific_01-09/raw/papalea_prosper_clock_2_1_1_seed_5094_fmri_raw_2401091850.csv')
#df1 <- read_csv('/Users/andypapale/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Documents - DNPLskinner/skinner/data/prolific/clock_v2.1_pilot/prolific_01-09/raw/papalea_prosper_clock_2_1_1_seed_5173_fmri_raw_2401091854.csv')
#df1 <- read_csv('/Users/andypapale/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Documents - DNPLskinner/skinner/data/prolific/clock_v2.1_pilot/prolific_01-09/raw/papalea_prosper_clock_2_1_1_seed_5815_fmri_raw_2401091851.csv')
df1 <- read_csv('/Users/andypapale/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Documents - DNPLskinner/skinner/data/prolific/clock_v2.1_pilot/prolific_01-09/raw/papalea_prosper_clock_2_1_1_seed_6520_fmri_raw_2401091851.csv')

#design_file <- '/Users/andypapale/clock2/Inquisit_design_files/Design-Matrix-152.csv'
#design_file <- '/Users/andypapale/clock2/Inquisit_design_files/Design-Matrix-868.csv'
#design_file <- '/Users/andypapale/clock2/Inquisit_design_files/Design-Matrix-1464.csv'
#design_file <- '/Users/andypapale/clock2/Inquisit_design_files/Design-Matrix-1752.csv'
#design_file <- '/Users/andypapale/clock2/Inquisit_design_files/Design-Matrix-2534.csv'
#design_file <- '/Users/andypapale/clock2/Inquisit_design_files/Design-Matrix-4938.csv'
#design_file <- '/Users/andypapale/clock2/Inquisit_design_files/Design-Matrix-5094.csv'
#design_file <- '/Users/andypapale/clock2/Inquisit_design_files/Design-Matrix-5173.csv'
#design_file <- '/Users/andypapale/clock2/Inquisit_design_files/Design-Matrix-5815.csv'
design_file <- '/Users/andypapale/clock2/Inquisit_design_files/Design-Matrix-6520.csv'
d1 <- read_csv(design_file)

df2 <- df1 %>% filter(trialcode == 'getValueVector')
df3 <- df1 %>% filter(trialcode == 'feedback')
rm(df1)
df2$RT <- df2$x_final
df0 <- inner_join(df2, d1,by=c('trial','RT'))

if (!file.exists(paste0('/Users/andypapale/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Documents - DNPLskinner/skinner/data/prolific/clock_v2.1_pilot/prolific_01-09/raw/seed-',df0$seed[1],'-validation'))){
  dir.create(paste0('/Users/andypapale/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Documents - DNPLskinner/skinner/data/prolific/clock_v2.1_pilot/prolific_01-09/raw/seed-',df0$seed[1],'-validation'))
}
setwd(paste0('/Users/andypapale/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Documents - DNPLskinner/skinner/data/prolific/clock_v2.1_pilot/prolific_01-09/raw/seed-',df0$seed[1],'-validation'))
for (iT in 1:240){
  pdf(paste0('seed-',as.character(df0$seed[1]),'-trial-',as.character(iT),'.pdf'),height=12,width=12)
  gg1 <- ggplot(df0 %>% filter(trial==iT), aes(x=RT,y=value)) + geom_point(color='black',size=10) + geom_line(color='black') +
    geom_point(aes(x=RT,y=y_final,color=subject),size=2) + 
    geom_point(data=df3 %>% filter(trial==iT), aes(x=rt_index,y=mag,color=subject),shape=2,size=8) + ylim(c(0,200)) + ggtitle(paste0('seed-',as.character(df0$seed[1]),' trial-',as.character(iT)))
  print(gg1)
  dev.off()
}
