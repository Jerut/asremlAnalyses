#Function to create response model for the residual
makeResidform<- function(staOutput){
  mods<- staOutput[[3]]
  uresids<- unique(mods$residual)
  residform<- c()
  for(i in 1:length(uresids)){
    ixR<-which(uresids[i]==mods$residual)
    if(length(ixR>0)){
      if(uresids[i]=='~NULL'){
        residform<- append(residform, paste("dsum(~id(trialnumb):id(Obsno) | trialnumb, levels = c(", 
                                            paste(ixR, collapse=","), "))", sep=""))      
      }else{
        residform<- append(residform, paste("dsum(", uresids[i]," | trialnumb, levels = c(", 
                                            paste(ixR, collapse=","), "))", sep=""))
      }
    }
  }
  residform<-paste(residform, collapse=" + ")
  
  residform<- as.formula(paste("~", residform))
  return(residform)
}
