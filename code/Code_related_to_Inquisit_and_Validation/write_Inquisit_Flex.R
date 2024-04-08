df <- t(read_csv('spread_vector-2455.csv'))
row_prefix = 'spread';
nR <- ncol(df) # these will be the items in each list
options("encoding" = "UTF-8") # encode in UTF-8 as suggested here https://forums.millisecond.com/Topic15777.aspx#15778
iqx_formatted_df <- data.frame(matrix(ncol=1, nrow=nR)) # import csv as data frame
df0 <- NULL
for (iR in 1:nR){ # start at 1st row
    if (iR==2){
      df0 <- paste0('<list ',row_prefix,as.character(iR-1),'>','\n','/ items = ( ',df[iR])
    } else if(iR>2 & iR < nR){
      df0 <- paste0(df0,', ',df[iR])
    }else if (iR==nR){
      df0 <- paste0(df0,',',df[iR],')','\n','/ selectionrate = always \n/ selectionmode=values.trialCount;','\n','</list>')
    }
}

write.table(df0 ,'spread_vector-2455.txt',row.names=F,col.names=F,quote=F)
options("encoding" = "native.enc") # change encoding back to native
