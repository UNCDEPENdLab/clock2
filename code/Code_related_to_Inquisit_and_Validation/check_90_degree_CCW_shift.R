# 2025-02-03 AndyP
# Check 90 degree CCW shift

library(tidyverse)
library(data.table)
library(lme4)
library(BAMBI)
library(pracma)

design_file <- '/Users/andypapale/clock2/Inquisit_design_files/Design-Matrix-152.csv'
d152 <- read_csv(design_file)
df_152 <- df0 %>% filter(seed=='152')
df_152 <- df_152 %>% mutate(resp_theta_deg = case_when(pos_shifted < 90 ~ (pos_shifted + 270),
                                                 pos_shifted >= 90 ~ (pos_shifted - 90)))
df_152_correct <- df_152 %>% mutate(RT = round(resp_theta_deg))
df_152_correct <- inner_join(d152,df_152_correct,by=c('trial','RT'))

hist(df_152_correct$value-df_152_correct$mag,breaks=100)

df_mismatch_152_correct = df_152_correct %>% filter(df_152_correct$value - df_152_correct$mag != 0)

hist(df_mismatch_152_correct$value - df_mismatch_152_correct$mag, breaks=100)

df_152_incorrect <- df_152 %>% mutate(RT = round(pos_shifted))
df_152_incorrect <- inner_join(d152,df_152_incorrect,by=c('trial','RT'))

hist(df_152_incorrect$value-df_152_incorrect$mag,breaks=100)

df_mismatch_152_incorrect <- df_152_incorrect %>% filter(df_152_incorrect$value - df_152_incorrect$mag != 0)

hist(df_mismatch_152_incorrect$value - df_mismatch_152_incorrect$mag, breaks=100)

                  