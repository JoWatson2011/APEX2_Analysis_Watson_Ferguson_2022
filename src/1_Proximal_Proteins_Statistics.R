## SETUP
set.seed(1)
library(tidyverse)
library(patchwork)
library(paletteer)
library(enrichR)
library(geneEnrichment)
source("src/functions/createHeatmap.R")
#source("src/functions/enrichmentFuns.R")

proteinGroups <-
  readRDS("data/proteinGroups_Flt.rds") 

experiments <-
  grep("LFQ intensity ", colnames(proteinGroups), value = T) %>%
  gsub("LFQ intensity ", "", .) %>%
  gsub("_40", "", .) %>%
  gsub("_P_", "_", .)

colnames(proteinGroups)[grep("LFQ intensity", colnames(proteinGroups))] <- experiments

# Impute missing values
lfq <- proteinGroups %>%
  dplyr::select(`Gene names`, grep("A_0[123]", colnames(proteinGroups))) %>%
  pivot_longer(cols = -1) %>%
  mutate(grp = gsub("_0[123]", "", name)) %>% 
  group_by(`Gene names`,grp) %>% 
  filter(sum(is.na(value))<2) %>% 
  ungroup() %>% 
  select(-grp) %>% 
  pivot_wider(id_cols = `Gene names`, values_fn=median) %>% 
  as.data.frame()

lfq[,2:ncol(lfq)] <-  
  apply(lfq[,2:ncol(lfq)], 2, function(i) 
    imputeLCMD::impute.QRILC(as.matrix(i), tune.sigma = 1.9)[[1]]
  )

lfq_normGFP <-  lfq %>%
  pivot_longer(-`Gene names`) %>% 
  filter(`Gene names` != "") %>% 
  # group_by(grp, `Gene names`) %>%
  # summarise(value = median(value), .groups = "keep") %>%
  mutate(
    grp = gsub("_0[123]", "", name),
    stim = ifelse(grepl("_C_", grp), "Ctrl", "FGF10")) %>% 
  group_by(stim, `Gene names`) %>% 
  mutate(norm = value/median(value[grepl("_GFPA_", grp)])) %>% 
  filter(!grepl("_GFPA_", grp)) %>% 
  group_by(`Gene names`) %>% 
  group_split() %>% 
  lapply( function(i){
    dat <- i[,c("grp", "norm")]
    av <- aov(norm ~ grp, data = dat)
    p <- summary(av)[[1]][1, 5]
    
    df <- cbind(i, p)
    
    return(df)
  }) %>% 
  do.call(rbind, .) %>% 
  mutate(adjp = p.adjust(p, method = "fdr")) %>% 
  filter(adjp < 0.05) %>% 
  group_by(grp, `Gene names`) %>%
  summarise(value = median(norm), .groups = "keep")


pg_sig <- filter(proteinGroups, `Gene names` %in% lfq_apex_normGFP$`Gene names`)

readr::write_csv(pg_sig, "results/data/APEX_PRO_SIG.csv")
