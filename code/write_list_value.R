# 2023-10-09 AndyP
# write inquisit list.values.appenditem for each rt into value vector

setwd('~/clock2')
# write Inquisit-formatted if-then syntax with data frames to select reward magnitude for clock 3.0
dfr <- data.frame(r=round(seq(from=0,to=359,length.out=360),0))
nT <- nrow(dfr) # number of timepoints
df0 <- NULL
for (iT in 1:nT){
    df0 <- paste0(df0,'list.values.appenditem(list.rt_',as.character(dfr$r[iT]),'.nextvalue);');
}
iqx_formatted_df[1,1] <- df0 
write.table(iqx_formatted_df,'value_list',row.names=F,col.names=F,quote=F)
options("encoding" = "native.enc") # change encoding back to native

dfr <- data.frame(r=round(seq(from=0,to=359,length.out=360),0))
nT <- nrow(dfr) # number of timepoints
df0 <- NULL
for (iT in 1:nT){
  df0 <- paste0(df0,'list.sorted_values.appenditem(list.rt_',as.character(dfr$r[iT]),'.nextvalue);');
}
iqx_formatted_df[1,1] <- df0 
write.table(iqx_formatted_df,'sorted_values_list',row.names=F,col.names=F,quote=F)
options("encoding" = "native.enc") # change encoding back to native