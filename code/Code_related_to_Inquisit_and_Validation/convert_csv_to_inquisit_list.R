# 2023-04-25 Andrew E. Papale
# Email: papalea@pitt.edu
# Convert each row of a csv into an Inquisit List
#' @param csv_file full path to a csv file containing values that you want to convert into Inquisit lists.  Each row will be formatted into a separate list.
#' 
#' @details 
#' The output of this function will be a text file formatted in Inquisit programming language as a static list.  You can either then try to import that file into Inquisit
#' or copy and paste the resulting code into an iqx file.
#' 
#' 
#' @export
#' @author Andrew E. Papale

convert_csv_to_inquisit_list <- function(csv_file = '', row_prefix ='rt_'){

  library(tidyverse)
  df <- read_csv(csv_file)
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
        df0 <- paste0(df0,',',df[iR,iC],')','\n','/ selectionrate = trial \n/ selectionmode=values.index;','\n','</list>')
      }
    }
    iqx_formatted_df[iR,1] <- df0 
  }
  write.table(iqx_formatted_df,'times.txt',row.names=F,col.names=F,quote=F)
  options("encoding" = "native.enc") # change encoding back to native
  
  
  }