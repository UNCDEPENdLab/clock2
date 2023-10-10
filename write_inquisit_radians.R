# write inquisit radians
# 2023-09-28 AndyP
# write inquisit expressions that encode radians from center of screen
setwd('~/clock2')
# write Inquisit-formatted if-then syntax with data frames to select reward magnitude for clock 3.0
dfr <- data.frame(r=round(seq(from=0,to=360,length.out=361),0))
nT <- nrow(dfr) # number of timepoints
# for each pair of values in nT, select an index of nV
options("encoding" = "UTF-8") # encode in UTF-8 as suggested here https://forums.millisecond.com/Topic15777.aspx#15778
iqx_formatted_df <- data.frame(matrix(ncol=1, nrow=1)) # import csv as data frame

df0 <- NULL
for (iT in 1:nT){
  df0 <- paste0(df0,'/ radians_angle',as.character(dfr$r[iT]),'=rad(',as.character(dfr$r[iT]),')\n','/ heightchange_angle',as.character(dfr$r[iT]),' = sin(expressions.Radians_angle',as.character(dfr$r[iT]),') * expressions.radius_px',
'\n','/ widthchange_angle',as.character(dfr$r[iT]),' = cos(expressions.Radians_angle',as.character(dfr$r[iT]),') * expressions.radius_px','\n')
}
iqx_formatted_df[1,1] <- df0 
write.table(iqx_formatted_df,'radians',row.names=F,col.names=F,quote=F)
options("encoding" = "native.enc") # change encoding back to native

df0 <- NULL
for (iT in 1:nT){
  if (iT==1){
    df0 <- paste0('<list degrees>\n/ items = (',as.character(dfr$r[iT]),', ')
  } else if (iT>1 && iT<nT){
    df0 <- paste0(df0,as.character(dfr$r[iT]),', ')
  } else if (iT==nT){
    df0 <- paste0(df0,as.character(dfr$r[iT]),')\n/ selectionrate = trial\n/ selectionmode = values.deg_index;\n</list>')
  }
}
iqx_formatted_df[1,1] <- df0 
write.table(iqx_formatted_df,'deg_list',row.names=F,col.names=F,quote=F)
options("encoding" = "native.enc") # change encoding back to native
