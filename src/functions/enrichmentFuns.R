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


simplifyGO <-  function (GORes, simplifyData = simplifyGOReqData(), simThresh = 0.4, adjpcol) 
{
  if (nrow(GORes) < 1) {
    NA
  }
  else {
    if(!("GOID" %in% colnames(GORes))){
      GORes <- merge(GORes, simplifyData$GOterms, by = "TERM") %>% 
        unique
    }
    sim <- GOSemSim::mgoSim(GORes$GOID, GORes$GOID, semData = simplifyData$semData, 
                            measure = "Rel", combine = NULL)
    sim[is.na(sim)] <- 0
    go1 <- go2 <- similarity <- NULL
    sim.df <- as.data.frame(sim)
    sim.df$go1 <- row.names(sim.df)
    sim.df <- gather(sim.df, go2, similarity, -go1)
    sim.df <- sim.df[!is.na(sim.df$similarity), ]
    sim.df <- sim.df[order(sim.df$similarity, decreasing = T), 
    ]
    sim.df <- sim.df[sim.df$go1 != sim.df$go2, ]
    sim.df <- sim.df[sim.df$similarity > simThresh, ]
    sim.df$remove <- apply(sim.df, 1, function(x) {
      if (x[1] %in% simplifyData$highFreqTerms) {
        return(x[1])
      }
      if (x[2] %in% simplifyData$highFreqTerms) {
        return(x[2])
      }
      else {
        return(NA)
      }
    })
    remove <- na.omit(sim.df$remove)
    sim.df <- sim.df[is.na(sim.df$remove), ]
    if (nrow(sim.df) == 0) {
      GORes.filt <- data.frame()
    }
    else {
      sim.df$go1.pval <- GORes[,adjpcol][match(sim.df$go1, 
                                               GORes$GOID)]
      sim.df$go2.pval <- GORes[,adjpcol][match(sim.df$go2, 
                                               GORes$GOID)]
      for (i in 1:nrow(sim.df)) {
        if (sim.df[i, "go1"] %in% remove) {
          next
        }
        if (sim.df[i, "go2"] %in% remove) {
          next
        }
        go1.pval <- sim.df[i, "go1.pval"]
        go2.pval <- sim.df[i, "go2.pval"]
        if (go1.pval == go2.pval) {
          go1 <- sim.df[i, "go1"]
          go2 <- sim.df[i, "go2"]
          if (go2 %in% simplifyData$childTerms[[go1]]) {
            remove <- c(remove, go2)
            next
          }
          else if (go1 %in% simplifyData$childTerms[[go2]]) 
            remove <- c(remove, go1)
          next
        }
        remove <- c(remove, sim.df[i, which.max(c(go1.pval, 
                                                  go2.pval))])
      }
      GORes.filt <- GORes[!GORes$GOID %in% remove, ]
    }
    return(GORes.filt)
  }
}

calculateEnrichment <- function(genes, enrichr_db, simplify = T, visualise = F, col = "#d2d2d2"){

  df <- enrichr(genes, databases = enrichr_db)
  df_flt <- filter(df[[1]], Adjusted.P.value < 0.05)
  
  if(grepl("GO", enrichr_db)){
    df_flt <-  df_flt %>% 
    separate(Term, c("Term", "GOID"), " \\(GO:") %>%
    mutate(GOID = paste0("GO:", gsub(")", "", GOID)))
    
    if(nrow(df_flt) > 0 & 
       simplify == T &
       grepl("GO_Biological_Process", enrichr_db)){
      colnames(df_flt)[colnames(df_flt) == "Term"] <- "TERM"
      df_flt <- simplifyGO(df_flt, adjpcol = "Adjusted.P.value") 
      colnames(df_flt)[colnames(df_flt) == "TERM"] <- "Term"
    }
  }
  
  if(grepl("Reactome", enrichr_db)){
    df_flt <- df_flt %>% 
    mutate(ID = paste0("R-HSA", gsub(".* R-HSA", "", Term)),
           Term = gsub("Homo sapiens R-HSA.*", "", Term)) 
    if(simplify == T & nrow(df_flt) > 0){
      reactHier_Signl <- readRDS("data/Reactome_Pathways_filtered.csv")
      
      df_flt <- df_flt %>% 
        filter(ID %in% c(reactHier_Signl$from, reactHier_Signl$to))
    }
  }
  
  if(visualise == T){
    gg <- visEnrichment(df_flt, col = col, xlab = enrichr_db) + xlab(enrichr_db)
    
    return(gg)
  }
  
  return(df_flt)
}

visEnrichment <- function(enriched, col = "#d2d2d2", xlab){
  library(ggplot2)
  gg <- enriched %>% 
    arrange(desc(Adjusted.P.value)) %>% 
    mutate(order = row_number()) %>% 
    ggplot(aes(x= reorder(Term, order), y = -log(Adjusted.P.value))) +
    geom_col(fill = col) +
    theme_bw() +
    coord_flip() +
    ylab("-log(FDR)") + 
    theme(plot.title.position = "plot") +
    geom_hline(yintercept = 2.995732, color = "red", linetype = "dotted")
  
  return(gg)
}

saveEnrichmentVis <- function(gg, path){
  library(ggplot2)
  if(length(unique(gg$data$Term)) < 25){
    ggsave(path, gg)
  } else if(length(unique(gg$data$Term)) < 100){
    ggsave(path, gg + theme(axis.text.y = element_text(size = 7)),
           height = 6, width = 8.3, units ="in")
  } else{
    ggsave(path, gg + theme(axis.text.y = element_text(size = 7)),
           height = 11.7, width = 8.3, units ="in")
  }
}