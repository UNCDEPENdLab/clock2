# 2023-12-19 AndyP
# Write Inquisit fmri timings

load('/Users/andypapale/clock2/clock2_fMRIOptITIs_570s_w_trials.rdata')

options("encoding" = "UTF-8")
dq0 <- NULL;
dz0 <- NULL;
dh0 <- NULL;
nR <- 240;

for (iS in 1:100){
  # convert to ms
  freeze <- c(t(freezemat_squeeze[iS,,]))*1000;
  ISI <- c(t(ISImat_squeeze[iS,,]))*1000;
  ITI <- c(t(ITImat_squeeze[iS,,]))*1000;
  
  for (iR in 1:nR){
    
    if (iR==1 && iS==1){
      #df0 <- paste0('<list erasure_locations>\n/ items = (',as.character(temp),',');
      dq0 <- paste0('<list ITI>\n/ items = (',as.character(ITI[iR]),',');
      dz0 <- paste0('<list ISI>\n/ items = (',as.character(ISI[iR]),',');
      dh0 <- paste0('<list preClockFreeze>\n/ items = (',as.character(freeze[iR]),',');
    } else if (iR==nR && iS==100){
      #df0 <- paste0(df0,as.character(temp),')\n/ selectionrate = always\n/ selectionmode = values.era_loc_index;\n</list>')
      dq0 <- paste0(dq0,as.character(ITI[iR]),')\n/ selectionrate = always\n/ selectionmode = values.ITI_index;\n</list>')
      dz0 <- paste0(dz0,as.character(ISI[iR]),')\n/ selectionrate = always\n/ selectionmode = values.ITI_index; \n</list>')
      dh0 <- paste0(dh0,as.character(freeze[iR]),')\n/ selectionrate = always\n/ selectionmode = values.ITI_Index; \n</list>')
    } else{
      #df0 <- paste0(df0,as.character(temp),',');
      dq0 <- paste0(dq0,as.character(ITI[iR]),',');
      dz0 <- paste0(dz0,as.character(ISI[iR]),',');
      dh0 <- paste0(dh0,as.character(freeze[iR]),',');
    }
  }
}

ITI_start_index <- seq(from=1,to=24000,by=240);

df0 <- NULL;
for (iT in 1:100){
  if (iT == 1){
    df0 <- paste0('<list ITI_index_start>\n/ items = (',as.character(ITI_start_index[iT]),',');
  } else if (iT > 1 && iT < 100){
    df0 <- paste0(df0,as.character(ITI_start_index[iT]),',');
  } else if (iT == 100){
    df0 <- paste0(df0,as.character(ITI_start_index[iT]),')\n/ selectionrate = always\n/ poolsize = 100\n/ selectionmode = random;\n</list>')
  }
}

write.table(df0,paste0('ITI_index_start','.txt'),row.names=F,col.names=F,quote=F)
write.table(dq0,paste0('ITI','.txt'),row.names=F,col.names=F,quote=F)
write.table(dz0,paste0('ISI','.txt'),row.names=F,col.names=F,quote=F)
write.table(dh0,paste0('preClockFreeze','.txt'),row.names=F,col.names=F,quote=F)
options("encoding" = "native.enc") # change encoding back to native
