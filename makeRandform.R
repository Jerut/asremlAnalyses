#Function to create response model for the random part
makeRandform<- function(staOutput, baseForm="~study+rep:study+vm(mgid, ainv)"){
  mods<- staOutput[[3]]
  urands<- unique(mods$random)
  urands<- unique(unlist(strsplit(as.character(urands), split="'+'")))
  urands<- urands[-na.omit(c(match('~vm(mgid, ainv)', urands), match('~mgid', urands)))]
  factrs<- gsub(" ", "", strsplit(urands, split="\\+")[[1]])
  factrs<- factrs[-grep('mgid', factrs)]
  randform<- c()
  for(i in 1:length(factrs)){
    ixr<- grep(factrs[i], mods$random)#add conditional rows
    randform<- append(randform, paste(paste('at(trialnumb,', ixr, "):", factrs[i], sep=""), 
                                        collapse=" + "))
    }
  randform<- paste(baseForm, paste(randform, collapse=" + "), sep="+")
  randform<- as.formula(randform)
  return(randform)
}

