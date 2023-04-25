# 2023-04-25 AndyP
# write Inquisit-formatted if-then syntax with data frames to select reward magnitude for clock 3.0
csv_file = '/Users/andypapale/clock2/2022-04-25-DesignFile.csv'
df <- read_csv(csv_file)
df_t <- data.frame(t=c(0,round(seq(from=0,to=5,length.out=360),2),5),ix=seq(from=0,to=361,by=1))
df <- read_csv(csv_file)
nT <- nrow(df_t) # number of timepoints
nV <- ncol(df) # number of values
nR <- nrow(df) # number of trials
# for each pair of values in nT, select an index of nV
options("encoding" = "UTF-8") # encode in UTF-8 as suggested here https://forums.millisecond.com/Topic15777.aspx#15778
iqx_formatted_df <- data.frame(matrix(ncol=1, nrow=nR)) # import csv as data frame

for (iR in 1:nR){
  df0 <- NULL
  for (iT in 1:nT-1){
    if (iT==1){
      df0 <- paste0('if (values.trialCount == ',iR,' ){ \n if (values.rt_shifted > ',df_t$t[iT], ' & values.rt_shifted <= ', df_t$t[iT+1],') {','\n','values.mag = ',df[iR,iT],'; \n')
    } else if (iT>1 & iT<nT-2){
      df0 <- paste0(df0,'else if (values.rt_shifted > ',df_t$t[iT], ' & values.rt_shifted <= ',df_t$t[iT+1],') { \n','values.mag = ',df[iR,iT],'; \n')
    } else if (iT == nT-1){
      df0 <- paste0(df0,'else if (values.rt_shifted < ',df_t$t[iT], ' & values.rt_shifted <= ',df_t$t[iT+1],') { \n','values.mag = ',df[iR,iT],'; \n',' }','\n } ')
    }
  }
  iqx_formatted_df[iR,1] <- df0 
}
write.table(iqx_formatted_df,'contingencies.txt',row.names=F,col.names=F,quote=F)
options("encoding" = "native.enc") # change encoding back to native
