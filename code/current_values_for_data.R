



options("encoding" = "UTF-8") # encode in UTF-8 as suggested here https://forums.millisecond.com/Topic15777.aspx#15778
iqx_formatted_df <- data.frame(matrix(ncol=1, nrow=nR)) # import csv as data frame
df0 <- NULL
for (iR in 1:360){ # start at 1st row
      if (iR==1){
        df0 <- paste0('list.curr_values.item(',as.character(iR),'),')
      } else if (iR==nR) {
        df0 <- df0 <- paste0(df0,'list.curr_values.item(',as.character(iR),')')
      } else {
        df0 <- paste0(df0,'list.curr_values.item(',as.character(iR),'),')
      }
}
write.table(df0,'current_values_for_data.txt',row.names=F,col.names=F,quote=F)
options("encoding" = "native.enc") # change encoding back to native