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

proteinGroups[,experiments] <- 
  apply(proteinGroups[,experiments], 2, function(i) 
    imputeLCMD::impute.QRILC(as.matrix(i), tune.sigma = 1.9)[[1]]
  )

lfq_apex <- proteinGroups %>%
  dplyr::select(`Gene names`, grep("A_0[123]", colnames(proteinGroups))) %>%
  pivot_longer(cols = -`Gene names`) %>% 
  group_by(`Gene names`) %>% 
  filter(sum(is.na(value)) < n())  %>% 
  mutate(`Gene names` = gsub(";.*", "", `Gene names`)) %>% 
  mutate(grp = gsub("_0[123]", "", name))


lfq_apex_normGFP <-  lfq_apex %>% 
  filter(`Gene names` != "") %>% 
  # group_by(grp, `Gene names`) %>%
  # summarise(value = median(value), .groups = "keep") %>%
  mutate(stim = ifelse(grepl("_C_", grp), "Ctrl", "FGF10")) %>% 
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

forhm_normGFP <- lfq_apex_normGFP  %>%  
  group_by(`Gene names`) %>% 
  mutate(value = (value - mean(value)) / sd(value)) %>% 
  pivot_wider(names_from = "grp", values_from = "value") %>% 
  column_to_rownames("Gene names") %>%
  select(`FGFR2-APEX 0` = R2A_R11G_C_A,
         `FGFR2-APEX + FGF10 40` = R2A_R11G_F_A,
         `RAB11-APEX + FGF10 40` = R2M_R11A_F_A
                    )
  #as.matrix()


## TWO WAYS TO DRAW A HEATMAP :)
# Redone for publication.
# Original retained for posperity / reproducibility
set.seed(1)
den_cl_normGFP <- hclust_heatmap(as.matrix(forhm_normGFP), k = 5)
cl <- cbind(
  forhm_normGFP,
  Cluster = den_cl_normGFP
) %>% 
  arrange(Cluster)

# pdf("results/figs/forPaper/FigureS2_ProHeatmap.pdf",
#     width = 40 / 25.4, height = 70 / 25.4, #in. to mm conversion
#     )
pheatmap::pheatmap(as.matrix(cl[,1:3]), #use the matrix without the cluster column
         color=colorRampPalette(c("navy", "white", "red"))(50),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         annotation_row = cl[,"Cluster", drop=F], #data frame of clusters above
         annotation_colors = list(Cluster = RColorBrewer::brewer.pal(max(den_cl_normGFP), "Set2")),
         show_rownames = FALSE,
         fontsize= 6,
         annotation_legend = F
         )
#dev.off() # Will manipulate legend etc. in illustrator


normGFP_cl <- data.frame(gene = names(den_cl_normGFP), cl = den_cl_normGFP) %>% 
  filter(cl %in% c(1,3,4,5)) %>% 
  mutate(cl_name = ifelse(cl == 5, "RE_Profile", 
                          ifelse(cl == 3, "RAB11_unique",
                                 ifelse(cl == 1, "Receptor_profile", 
                                        ifelse(cl == 4, "FGFR2_F_unique", NA)))))


cols <- sapply(1:length(unique(normGFP_cl$cl_name)), function(i){
  paletteer_d("ggsci::category10_d3", length(unique(normGFP_cl$cl_name)))[i]
})
names(cols) <- unique(normGFP_cl$cl_name)

## GOBP Enrichment
normGFP_cl_GOBP <- lapply(unique(normGFP_cl$cl_name), function(i){
  genes <- normGFP_cl[normGFP_cl$cl_name == i, ]$gene
  
  gg <- calculateEnrichment(genes, "GO_Biological_Process_2018", visualise = T, col=cols[i])
  gg <- gg + ggtitle(i)
})

## KEGG Pathway enrichment
normGFP_cl_KEGG <- lapply(unique(normGFP_cl$cl_name), function(i){
  genes <- normGFP_cl[normGFP_cl$cl_name == i, ]$gene
  
  gg <- calculateEnrichment(genes, "KEGG_2019_Human", visualise = T, col=cols[i])
  gg <- gg + ggtitle(i)
})

pg_sig <- filter(proteinGroups, `Gene names` %in% lfq_apex_normGFP$`Gene names`)
pg_sig <- merge(pg_sig, 
                 normGFP_cl[,c("gene", "cl_name")], 
                 by.x = "Gene names", by.y = "gene")
readr::write_csv(pg_sig, "results/data/APEX_PRO_SIG_CLUSTERED.csv")
