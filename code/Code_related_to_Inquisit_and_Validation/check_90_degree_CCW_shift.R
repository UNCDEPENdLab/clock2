# 2025-02-03 AndyP
# Check 90 degree CCW shift

library(tidyverse)
library(data.table)
library(lme4)
library(BAMBI)
library(pracma)

design_file <- '/Users/andypapale/clock2/Inquisit_design_files/Design-Matrix-152.csv'
d152 <- read_csv(design_file)
df1 <- df0 %>% filter(seed=='152')
df1 <- df1 %>% mutate(resp_theta_deg = case_when(pos_shifted < 90 ~ (pos_shifted + 270),
                                                 pos_shifted >= 90 ~ (pos_shifted - 90)))
df2 <- df1 %>% mutate(RT = round(resp_theta_deg))
df_correct <- inner_join(d152,df2,by=c('trial','RT'))

hist(df_correct$value-df_correct$mag,breaks=100)

df_mismatch_152_correct = df_correct %>% filter(df_correct$value - df_correct$mag != 0)

df3 <- df1 %>% mutate(RT = round(pos_shifted))
df_incorrect <- inner_join(d152,df3,by=c('trial','RT'))

hist(df_incorrect$value-df_incorrect$mag,breaks=100)

df_mismatch_152_incorrect <- df_incorrect %>% filter(df_incorrect$value - df_incorrect$mag != 0)

                  