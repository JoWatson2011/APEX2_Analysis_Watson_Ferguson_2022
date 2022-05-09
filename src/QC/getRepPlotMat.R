getRepPlotMat	<- function(inMat){
  exps			<- sort(unique(inMat[,"Experiment"]))
  logVec		<- vector()
  
  #calculate cumulative probability that there is a site on that peptide
  for(i in 1:nrow(inMat)){
    temp	<- unlist(strsplit(inMat[i,"Phospho (STY) Probabilities"], "\\)"))
    if(length(temp)==0){
      temp <- "A(0"
    }			
    if(temp[1]==""){
      temp <- as.character(c(0,0))
    }
    
    temp	<- unlist(strsplit(temp, "\\("))
    temp	<- temp[seq(2, length(temp),2)]
    
    
    if(sum(as.numeric(temp) > 0.75) ==  inMat[i, "Phospho (STY)"]){
      logVec	<- c(logVec, i)
    }
  }
  
  #subset inMat so it only contains peptides with
  # probability > 0.75 of a mod
  inMat	<- inMat[logVec,]
  
  #modified peptides sequence
  unPeps		<- unique(inMat[,"Modified sequence"])
  repMat		<- matrix(NA, length(unPeps), length(exps)*3)
  rownames(repMat)	<- unPeps
  colnames(repMat)	<- c(paste("ML", exps, sep="_"), paste("HL", exps, sep="_"),  paste("HM", exps, sep="_"))
  
  for(i in 1:nrow(repMat)){
    temp	<- rbind(inMat[inMat[,"Modified sequence"] == rownames(repMat)[i],])
    for(j in 1:length(exps)){
      temp2	<- rbind(temp[temp[,"Experiment"] == exps[j],c("Ratio M/L normalized",
                                                           "Ratio H/L normalized",
                                                           "Ratio H/M normalized")])
      if(nrow(temp2)>0){
        repMat[i,j]				<- median(log2(as.numeric(temp2[,1])), na.rm=T)
        repMat[i,j+length(exps)]	<- median(log2(as.numeric(temp2[,2])), na.rm=T)
        repMat[i,j+2*length(exps)]	<- median(log2(as.numeric(temp2[,3])), na.rm=T)
      }
    }
  }
  
  return(repMat)		
}