###########
#' singleTrials
#' Function to perform single trial analysis, including model selection, for a group of studies of various designs
#' models include pedigree relationship and also spatial analysis
#' 
#' Parameters
#' @param dat is a data frame of the studies 
#' @param ped is a data frame of the pedigree file, ordered using pedigremm
#' @param trialvar is the name of the column for the study name
#' @param designvar is the name of the column for the design variable
#' @param trait is the name of the column for the trait
#' @param lwrlim is the minimum value of the trait considered to be a valid value
#' @param uprlim is the maximum value of the trait considered to be a valid value
#' @param idvar is the name of the colum for the id variable
#' @param missingHillsvar is the name of the column for the missing hills variable, if not recorded set to NULL
#' @param rowblk_var is the name of the colum for the row block variable, if not relevant set to NULL
#' @param colblk_var is the name of the column for the column block variable, if not relevant set to NULL
#' @param blk_var is the name of the column for the block variable, if not relevant set to NULL
#' @param colcoord_var is the name of the column for the column coordinates, if not known set to NULL
#' @param rowcoord_var is the name of the column for the row coordinates, if not known set to NULL 
#' @param workspace is the size of the workspace used in asreml for fitting the modesl
#' @param pworkspace is the size of the workspace used in asreml for prediction
#' 
#' @return a list containing 3 object, first the trial data, second a dataframe of resut outputs, third a summary of the models fit
#' #############

singleTrials<- function(dat=dat, ped=ped, trialvar='study', designvar='Design', 
                        trait='GRYLD_MC14', lwrlim=0, uprlim=15,
                        idvar= 'mgid', missingHillsvar= 'missinghills', 
                        rep_var='rep', colblk_var='col', 
                        rowblk_var= 'row', blk_var=NULL, colcoord_var='col', 
                        rowcoord_var='row', workspace='2gb', pworkspace='2gb'){
  require(asreml)
  require(pedigreemm)
  require(nadiv)
  
  #function to modify column name 
  modcol<- function(dat, oldnm, newnm){
    if(is.null(oldnm)){ dat[,newnm]=NA 
    }else{ colnames(dat)[match(oldnm, colnames(dat))]<- newnm}
    return(dat)
  }

  #remove imlausable values 
  if(length(which(dat$TOI<lwrlim))>0){dat[which(dat$TOI>uprlim),'TOI']<- NA}
  if(length(which(dat$TOI<lwrlim))>0){dat[which(dat$TOI<lwrlim),'TOI']<- NA}
  
  #modify column names for rep, design, missing hill, and id variables
  dat<- modcol(dat, missingHillsvar, 'missinghills')
  dat<- modcol(dat, idvar, 'mgid')
  dat<- modcol(dat, rep_var, 'rep')
  dat<- modcol(dat, designvar, 'Design')
  
  #rename column names for blocks
  dat<- modcol(dat, colblk_var, 'colB')
  dat<- modcol(dat, rowblk_var, 'rowB')
  dat<- modcol(dat, blk_var, 'block')
  
  #add rowcoord variables as seperate columns if necessary, or rename them
  if(rowcoord_var==rowblk_var){
    row<- dat[,'rowB']
    dat<- data.frame(dat, row)
  }else{
    dat<- modcol(dat, rowcoord_var, 'row')
  }
  if(colcoord_var==colblk_var){
    col<- dat[,'colB']
    dat<- data.frame(dat, col)
  }else{
    dat<- modcol(dat, colcoord_var, 'col')
  }
  
  #format phenotypic data
  TOI<- dat[,trait]
  dat<- data.frame(dat, TOI)
  dat$mgid<- as.character(dat$mgid)
  
  #set asreml options
  asreml.options(maxit=100,pworkspace=pworkspace, workspace=workspace)
  
  #AIC function
  AIC<- function(m){
    aic<- summary(m)$aic[1]
    return(aic)
  }
  
  #--------------------------------------------start of the trial-by-trial analysis
  ustud<- unique(dat$study)
  for(i in 1:length(ustud)){
    
    #-------subset trial
    ix<- which(dat$study==ustud[i])
    trial<- droplevels.data.frame(dat[ix,])    

    #-------Fix the row-col coordinates if they exist
    hascoord<- TRUE
    if(length(na.omit(trial$row))>0 & length(na.omit(trial$col))>0){
      mnrow<- min(trial$row, na.rm=TRUE)
      mncol<- min(trial$col, na.rm=TRUE)
      mxrow<- max(trial$row, na.rm=TRUE)
      mxcol<- max(trial$col, na.rm=TRUE)
      expecrow<- rep(c(mnrow:mxrow), length(mncol:mxcol))
      expeccol<- sort(rep(c(mncol:mxcol), length(mnrow:mxrow)))
      rc<- data.frame(row=expecrow, col=expeccol)
      trial<- merge(rc, trial, by=c('row', 'col'), all=TRUE, sort=FALSE)
      if(nrow(trial)>nrow(rc) ){
        hascoord=FALSE
      }
    }
    if(hascoord){
      trial$row<- as.factor(trial$row)
      trial$col<- as.factor(trial$col)
    }
    
    #---------get rep info
    if(length(na.omit(trial$row))>0){
      trial$rep<- as.character(trial$rep)
      trial$rep<- as.factor(trial$rep)
      tb<- table(trial$rep)
      repratio<- paste(tb, collapse=":")
      rep_pct_diff<- abs(tb-max(tb))/max(tb)
      if(length(which(rep_pct_diff>0.2))>0){
        rep_problem<- TRUE
      }else{
        rep_problem<- FALSE
      }
    }else{
      rep_problem<- TRUE
    }
    
    #----------check if missing hills should be used as a covariate
    pcthasmiss<- length(which(trial$missinghills>0))/nrow(trial)
    if(pcthasmiss>.2){
      useMHcovar=TRUE #if more than 20% of plots are affected use a covariate
    }else{
      useMHcovar=FALSE
    }
    
    #----------get Ainverse
    ped_sub<- nadiv::prunePed(ped, phenotyped=trial$mgid)
    ainv<- asreml::ainverse(ped_sub)
    assign("ainv", ainv, envir = .GlobalEnv) 

    #set mgid to factor
    trial$mgid<- as.factor(trial$mgid)
    
    #-----assemble the fixed part of the formula
    fxobjs<- list()
    fxobjs[1]<-  'TOI~1'
    if(useMHcovar){
      fxobjs[1]<- paste(fxobjs[1], "missinghills", sep="+") 
    }
    if(!rep_problem){
      fxobjs[2]<- paste(fxobjs[1], "rep", sep="+")
    }
    fxobjs<- lapply(fxobjs, as.formula)
    
    #-----assemble the random part of the formula
    rndobjs<- list()
    rndobjs[1]<-  '~vm(mgid, ainv)'
    
    #if there are row and column blocks
    if(trial[1,designvar]=='P-REP' | trial[1,designvar]=='Row-Column'){
      rndobjs[2]<-  paste(rndobjs[1], 'rowB', sep="+")
      rndobjs[3]<-  paste(rndobjs[1], 'colB', sep="+")
      rndobjs[4]<-  paste(rndobjs[2], 'colB', sep="+")
      
      if(!rep_problem){
        rndobjs[5]<-  paste(rndobjs[1], 'rowB:rep', sep="+")
        rndobjs[6]<-  paste(rndobjs[1], 'colB:rep', sep="+")
        rndobjs[7]<-  paste(rndobjs[5], 'colB:rep', sep="+")
      }
    }
    
    #if the design is alpha-lattice
    if(trial[1,designvar]=='Alpha-Lattice' | trial[1,designvar]=='Augmented Alpha-Lattice'){
      if(!rep_problem){
        rndobjs[3]<-  paste(rndobjs[1], 'block:rep', sep="+")
      }
    }
    rndobjs<- lapply(rndobjs, as.formula)
    
    #-----assemble the error structure part of the formula
    Robjs<- list()
    if(hascoord){
      Robjs[1]<- '~ar1(col):ar1(row)'
      Robjs[2]<- '~idv(col):ar1(row)'
      Robjs[3]<- '~ar1(col):idv(row)'
      Robjs<- lapply(Robjs, as.formula)
    }
    
    #----Fit models and save a list of objects
    modobjs<- list()
    for(j in 1:length(fxobjs)){
      for(k in 1:length(rndobjs)){
        m<- asreml(fixed=fxobjs[[j]], random= rndobjs[[k]], data=trial)
        
        #update until convergence
        counter<- 0
        while(length(which(summary(m)$varcomp[,5]>=1))>0 & counter<11){
          counter<- counter+1
          m<- update(m)
        }
        
        modobjs[[(j*k)]]<- m
      }
    }
    
    #------ fit spatial models if coordinate data are avaliable
    if(hascoord){
      numb<-length(modobjs)
      for(l in 1:length(modobjs)){
        for(n in 1:length(Robjs)){
          m<- modobjs[[l]]
          m<- update(m, residual=Robjs[[n]])
          numb<- numb+1
          
          #update until convergence
          counter<- 0
          while(length(which(summary(m)$varcomp[,5]>=1))>0 & counter<11){
            counter<- counter+1
            m<- update(m)
          }
          modobjs[[numb]]<- m
        }
      }
    }
    
    #select the best model based on AIC
    aicvec<- unlist(lapply(modobjs, AIC))
    model<-modobjs[[which(aicvec==min(aicvec))]]
    
    #get selected model formula info
    resid<- paste(as.character(unlist(model$formulae$residual)), collapse="")
    fixed<- paste(as.character(unlist(model$formulae$fixed))[-c(1:3)], collapse="+")
    if(fixed==""){
      fixed="TOI~1"
    }
    random<- paste(as.character(unlist(model$formulae$random)), collapse="")
    modinfo<- data.frame(study= ustud[i], fixed, random, residual=resid)
    
    #get blups and predicted values
    pv<- predict(model, classify='mgid', levels= list(mgid=trial$mgid))
    predtab<- pv$pvals
    rcf<- model$coef$random
    EBVs<- rcf[grep('mgid', row.names(rcf)),]
    names(EBVs)<- as.character(matrix(unlist(strsplit(names(EBVs), split="_")), nrow=2)[2,])
    
    #get inbreeding coefficients
    p<- pedigreemm::editPed(ped_sub[,2], ped_sub[,3],ped_sub[,1])
    P<- pedigreemm::pedigree(p[,2], p[,3],p[,1])
    I<- pedigreemm::inbreeding(P)
    names(I)<- p[,1]
    
    #calculate reliability
    PEV<- predtab$std.error^2
    vc<- summary(model)$varcomp
    Vg<- vc[match('vm(mgid, ainv)', row.names(vc)),'component']
    Vgs<- Vg*I[predtab$mgid]+1
    Reliability<- 1-PEV/Vgs
    
    #Get table of values
    rslts<- data.frame(study= ustud[i], predtab, EBVs=EBVs[predtab$mgid], Reliability)
    
    #average reliability
    h2<- mean(Reliability)
    
    #add to modinfo
    modinfo<- data.frame(modinfo, h2)
    cat(h2, "\n")
    
    #save trial data as is and results table
    if(i==1){
      trial_all<- trial
      results_all<- rslts
      modinfo_all<- modinfo
    }else{
      trial_all<- rbind(trial_all, trial)
      results_all<- rbind(results_all, rslts)
      modinfo_all<- rbind(modinfo_all, modinfo)
    }
  }#END of trial by trial loop
  return(list(trial_all=trial_all,results_all=results_all,modinfo_all=modinfo_all))
}
