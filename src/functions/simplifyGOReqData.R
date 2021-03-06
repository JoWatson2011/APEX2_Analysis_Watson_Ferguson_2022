# Gather data required for simplify() function, to reduce
# semantically similar GO terms

simplifyGOReqData <- function(){
  
  if(length(grep("simplifyGOreqdata.rds", list.files("data/"))) < 1) {

  
  gene2GOi<-AnnotationDbi::select(org.Hs.eg.db,keys(org.Hs.egGO2EG),c("ENTREZID", "SYMBOL"),  "GOALL") %>% filter(ONTOLOGYALL == "BP")
  gene2GO<-unstack(gene2GOi[,c(1,5)])
  
  GOterms <- AnnotationDbi::select(GO.db, keys=gene2GOi$GOALL, columns=c("TERM"), keytype = "GOID")
  
  GO2Gene<-unstack(stack(gene2GO)[2:1])
  
  freq <- sapply(GO2Gene,length)
  freqCutOff <- length(gene2GO)*0.05
  highFreqTerms <- names(freq[freq>freqCutOff])
  
  semData <- godata("org.Hs.eg.db","SYMBOL","BP")
  
  childTerms <- as.list(GOBPCHILDREN)
  
  saveRDS(list(gene2GO = gene2GO,
                        GOterms=GOterms,
                        GO2Gene=GO2Gene,
                        highFreqTerms=highFreqTerms,
                        semData = semData,
                        childTerms=childTerms),
          "data/simplifyGOreqdata.rds")
  
  return(list(gene2GO = gene2GO,
              GOterms=GOterms,
              GO2Gene=GO2Gene,
              highFreqTerms=highFreqTerms,
              semData = semData,
              childTerms=childTerms))
  
  }else{
    return(readRDS("data/simplifyGOreqdata.rds"))
  }
}
