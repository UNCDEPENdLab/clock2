# 2023-12-19 AndyP
# Write Inquisit fmri timings

load('/Users/andypapale/clock2/clock2_fMRIOptITIs_570s_w_trials.rdata')

options("encoding" = "UTF-8")
dq0 <- NULL;
dz0 <- NULL;
dh0 <- NULL;
nR <- 240;

freeze <- c(t(freezemat_squeeze[1,,]))*1000;
ISI <- c(t(ISImat_squeeze[1,,]))*1000;
ITI <- c(t(ITImat_squeeze[1,,]))*1000;

for (iR in 1:nR){
  
  if (iR==1){
    #df0 <- paste0('<list erasure_locations>\n/ items = (',as.character(temp),',');
    dq0 <- paste0('<list ITI>\n/ items = (',as.character(ITI[iR]),',');
    dz0 <- paste0('<list ISI>\n/ items = (',as.character(ISI[iR]),',');
    dh0 <- paste0('<list preClockFreeze>\n/ items = (',as.character(freeze[iR]),',');
  } else if (iR > 1 && iR < nR){
    #df0 <- paste0(df0,as.character(temp),',');
    dq0 <- paste0(dq0,as.character(ITI[iR]),',');
    dz0 <- paste0(dz0,as.character(ISI[iR]),',');
    dh0 <- paste0(dh0,as.character(freeze[iR]),',');
  } else if (iR==nR){
    #df0 <- paste0(df0,as.character(temp),')\n/ selectionrate = always\n/ selectionmode = values.era_loc_index;\n</list>')
    dq0 <- paste0(dq0,as.character(ITI[iR]),')\n/ selectionrate = always\n/ selectionmode = values.trial;\n</list>')
    dz0 <- paste0(dz0,as.character(ISI[iR]),')\n/ selectionrate = always\n/ selectionmode = values.trial; \n</list>')
    dh0 <- paste0(dh0,as.character(freeze[iR]),')\n/ selectionrate = always\n/ selectionmode = values.trial; \n</list>')
  }
}
#write.table(df0,paste0('era_loc-',as.character(rob_grid$iteration),'.txt'),row.names=F,col.names=F,quote=F)
write.table(dq0,paste0('ITI','.txt'),row.names=F,col.names=F,quote=F)
write.table(dz0,paste0('ISI','.txt'),row.names=F,col.names=F,quote=F)
write.table(dh0,paste0('preClockFreeze','.txt'),row.names=F,col.names=F,quote=F)
options("encoding" = "native.enc") # change encoding back to native
