# 2023-04-25 AndyP
# write Inquisit-formatted if-then syntax with data frames to select reward magnitude for clock 3.0
df_t <- data.frame(t=round(seq(from=0,to=5000,length.out=361),0),ix=seq(from=0,to=360,by=1))
nT <- nrow(df_t) # number of timepoints
# for each pair of values in nT, select an index of nV
options("encoding" = "UTF-8") # encode in UTF-8 as suggested here https://forums.millisecond.com/Topic15777.aspx#15778
iqx_formatted_df <- data.frame(matrix(ncol=1, nrow=1)) # import csv as data frame

  df0 <- NULL
  for (iT in 1:nT-1){
    if (iT==1){
      df0 <- paste0('if (values.rt_shifted > ',df_t$t[iT], ' && values.rt_shifted <= ', df_t$t[iT+1],') {','\n','values.mag = ','list.rt_',as.character(df_t$ix[iT]),'.currentvalue','; \n')
    } else if (iT>1 & iT<nT-2){
      df0 <- paste0(df0,'} else if (values.rt_shifted > ',df_t$t[iT], ' && values.rt_shifted <= ',df_t$t[iT+1],') { \n','values.mag = ','list.rt_',as.character(df_t$ix[iT]),'.currentvalue','; \n')
    } else if (iT == nT-1){
      df0 <- paste0(df0,'} else if (values.rt_shifted < ',df_t$t[iT], ' && values.rt_shifted <= ',df_t$t[iT+1],') { \n','values.mag = ','list.rt_',as.character(df_t$ix[iT]),'.currentvalue','; \n',' }')
    }
  }
  iqx_formatted_df[1,1] <- df0 
  write.table(iqx_formatted_df,'contingencies.txt',row.names=F,col.names=F,quote=F)
  options("encoding" = "native.enc") # change encoding back to native
