###########
#' estBV
#' Function to estimate breeding values based using pedigree 
#' 
#' Parameters
#' @inputdir is a string of the directory containing the input data files 
#' @stage1file is a string of the file name containing the input data
#' @pedfile is a string of the file name containing the pedigree data
#' @covar is a character for the name of the covariate to be used, if none put NULL
#' @remlF90dir is a character for the name directory containing the REMLF90 programs and parameter files
#' @resultsLabel is a character string that will be included in the results table as metadata
#' 
#' @return a dataframe with the EBVs for each individual
#' #############

estBV<- function(inputdir= "/Volumes/GoogleDrive/My Drive/BV estimation 2019",
                 stage1file= 'Results stage 1 all GRYLD 2019.csv', pedfile= "Pedigree March 21 2019.csv",
                 covar= NULL, remlF90dir= "~/Documents/Execute REML F90", resultsLabel= 'Yield'){
  
  #set directory
  setwd(inputdir) 
  
  #read and format phenotypic data
  dat<- read.csv(stage1file)
  dat<- dat[,c(1:4,match(covar, colnames(dat)))]
  dat[,2]<- as.numeric(dat[,2])
  dat[,1]<- paste("f", dat[,1], sep="")
  dat[,4]<- dat[,4]-mean(dat[,4], na.rm=TRUE) #center Y
  if(!is.null(covar)){
    dat[,5]<- dat[,5]-mean(dat[,5], na.rm=TRUE) #center covariate
  }
  
  #read and format pedigree data
  pd<- read.csv(pedfile)
  pd0<- pd
  pd[,1]<- paste("f", pd[,1], sep="")
  pd[,2]<- paste("f", pd[,2], sep="")
  pd[,3]<- paste("f", pd[,3], sep="")
  pd[which(is.na(pd0[,2])),2]<- 0
  pd[which(is.na(pd0[,3])),3]<- 0
  
  #write files to program directory and and execute program
  setwd(remlF90dir)
  write.table(pd, file='pd.test', col.names=F, sep=" ", row.names=F, quote=F)
  write.table(dat, file='dat.test', col.names=F, sep=" ", row.names=F, quote=F)
  if(!is.null(covar)){
    system2('./bashscript_cov.sh')
  }else{
    system2('./bashscript.sh')
  }
  
  #get read the results 
  sl<- read.table('solutions', header=F, skip=1, fill=TRUE)
  if(!is.null(covar)){
    sl<- sl[which(sl[,2]==3),]
  }else{
    sl<- sl[which(sl[,2]==2),]
  }
  
  #add original ids to results
  ped<- read.table('renadd02.ped', header=F, skip=1, fill=TRUE)
  ped<- ped[match(pd[,1], ped$V10),]
  
  #format the results
  colnames(ped)[1]<- 'id'
  colnames(sl)[3]<- 'id'
  colnames(sl)[4]<- 'blup'
  sl2<- merge(sl, ped, by='id')
  sl3<- sl2[,c('blup', 'V10')]
  colnames(sl3)<- c('EBV', 'mgids')
  sl3$mgids<- gsub("f", "", sl3$mgids)
  rslt<- data.frame(resultsLabel, sl3)
  
  return(rslt)
}


