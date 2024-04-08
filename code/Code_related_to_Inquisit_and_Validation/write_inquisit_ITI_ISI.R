# 2023-10-02 AndyP
# Generate ISI/ITI lengths for #clock2

library(tidyverse)
setwd('~/Inquisit Code/EEG_clock/Clock_v2/')

df <- t(read_csv('ISI.csv'))
row_prefix = 'ISI_'
nR <- nrow(df) # these will be separate lists
nC <- ncol(df) # these will be the items in each list
options("encoding" = "UTF-8") # encode in UTF-8 as suggested here https://forums.millisecond.com/Topic15777.aspx#15778
iqx_formatted_df <- data.frame(matrix(ncol=1, nrow=nR)) # import csv as data frame
df0 <- NULL
for (iR in 1:nR){ # start at 1st row
  for (iC in 2:nC){ # start at 2, ignore col names
    if (iC==2){
      df0 <- paste0('<list ',row_prefix,as.character(iR-1),'>','\n','/ items = ( ',df[iR,iC])
    } else if(iC>2 & iC < nC){
      df0 <- paste0(df0,', ',df[iR,iC])
    }else if (iC==nC){
      df0 <- paste0(df0,',',df[iR,iC],')','\n','/ selectionrate = trial \n/ selectionmode=values.isi_index;','\n','</list>')
    }
  }
  iqx_formatted_df[iR,1] <- df0 
}
write.table(iqx_formatted_df,'ISI.txt',row.names=F,col.names=F,quote=F)
options("encoding" = "native.enc") # change encoding back to native

df <- t(read_csv('ITI.csv'))
row_prefix = 'ITI_'
nR <- nrow(df) # these will be separate lists
nC <- ncol(df) # these will be the items in each list
options("encoding" = "UTF-8") # encode in UTF-8 as suggested here https://forums.millisecond.com/Topic15777.aspx#15778
iqx_formatted_df <- data.frame(matrix(ncol=1, nrow=nR)) # import csv as data frame
df0 <- NULL
for (iR in 1:nR){ # start at 1st row
  for (iC in 2:nC){ # start at 2, ignore col names
    if (iC==2){
      df0 <- paste0('<list ',row_prefix,as.character(iR-1),'>','\n','/ items = ( ',df[iR,iC])
    } else if(iC>2 & iC < nC){
      df0 <- paste0(df0,', ',df[iR,iC])
    }else if (iC==nC){
      df0 <- paste0(df0,',',df[iR,iC],')','\n','/ selectionrate = trial \n/ selectionmode=values.iti_index;','\n','</list>')
    }
  }
  iqx_formatted_df[iR,1] <- df0 
}
write.table(iqx_formatted_df,'ITI.txt',row.names=F,col.names=F,quote=F)
options("encoding" = "native.enc") # change encoding back to native