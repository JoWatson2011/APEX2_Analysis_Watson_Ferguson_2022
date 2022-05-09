set.seed(1)
library(tidyverse)
library(patchwork)
library(paletteer)
library(enrichR)
library(geneEnrichment)
library(UpSetR)
source("src/functions/createHeatmap.R")

# Import data
sty <-
  readRDS("data/sty_Flt.rds") %>% 
  mutate(`Gene names` = gsub(";.*", "", `Gene names`),
         id_site = paste0(`Gene names`, "_", `Amino acid`, Position))

# Impute missing values in data
sty[,grep("Intensity", colnames(sty))] <-  
  apply(sty[,grep("Intensity", colnames(sty))], 2, function(i) 
    imputeLCMD::impute.QRILC(as.matrix(i), tune.sigma = 1.9)[[1]]
  )

experiments <-
  grep("Intensity ", colnames(sty), value = T) %>%
  gsub("Intensity ", "", .) %>%
  gsub("_40", "", .) %>%
  gsub("_Ph_", "_", .)

colnames(sty)[grep("Intensity ", colnames(sty))] <-
  experiments

# Data in long format
lfq <- sty %>%
  dplyr::select(id_site, grep("A_0[123]", colnames(sty))) %>%
  pivot_longer(cols = -1) %>%
  mutate(grp = gsub("_0[123]", "", name))

### Check overlap of sites identied between experiments
sty_nonimp <-
  readRDS("data/sty_Flt.rds") 
colnames(sty_nonimp)[grep("Intensity ", colnames(sty_nonimp))] <-
  experiments
forUpset <-sty_nonimp %>% 
  select(`Gene names`, grep("A_0[123]", colnames(.))) %>% 
  pivot_longer(cols = -`Gene names`) %>% 
  mutate(grp = gsub("_0[123]", "", name)) %>% 
  group_by(`Gene names`, grp) %>% 
  filter(sum(is.na(value)) < 2) %>% 
  ungroup() %>% 
  select(grp, `Gene names`) %>% 
  unique() %>% 
  mutate(values = 1) %>% 
  pivot_wider(names_from = grp, values_from = values) %>% 
  mutate_if(is.numeric, ~ifelse(is.na(.), 0, 1)) %>% 
  tibble::column_to_rownames("Gene names")

upset(forUpset, 
      nsets = 5,
      sets = colnames(forUpset),
      empty.intersections = "on", 
      order.by = "freq") 


# Differences between GFP control + stimulated

df <- sty[,"id_site", drop = F]
df$GFP_p <- apply(sty[,grep("GFPA_.*_A_0[123]", colnames(sty))],
                  1,
                  function(i){
                    GFPC <- c(i["R2M_GFPA_C_A_01"],
                              i["R2M_GFPA_C_A_02"],
                              i["R2M_GFPA_C_A_03"]
                    )
                    if(class(GFPC)== "list"){
                      GFPC <- do.call(c, GFPC)
                    }
                    GFPF <- c(i["R2M_GFPA_F_A_01"],
                              i["R2M_GFPA_F_A_02"],
                              i["R2M_GFPA_F_A_03"]
                    )
                    if(class(GFPF)== "list"){
                      GFPF <- do.call(c, GFPF)
                    }
                    #
                    t <- try(t.test(GFPC, GFPF), silent = T)
                    #
                    if(class(t) == "try-error"){
                      return(NA)
                    }else{
                      return(t$p.value)
                    }
                  })

nrow(df[df$GFP_p <= 0.05,])

# Normalise to GFP-APEX2 of same time point
lfq_normGFP <-  lfq %>% 
  filter(id_site != "") %>% 
  # group_by(grp, `Gene names`) %>%
  # summarise(value = median(value), .groups = "keep") %>%
  mutate(stim = ifelse(grepl("_C_", grp), "Ctrl", "FGF10")) %>% 
  group_by(stim, id_site) %>% 
  mutate(norm = value/median(value[grepl("_GFPA_", grp)])) %>% 
  filter(!grepl("_GFPA_", grp)) %>% 
  group_by(id_site) %>% 
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
  group_by(grp, id_site) %>%
  summarise(value = median(norm), .groups = "keep")

#Z-score for clustering
forhm_normGFP <- lfq_normGFP  %>%  
  group_by(id_site) %>% 
  mutate(value = (value - mean(value)) / sd(value)) %>% 
  pivot_wider(names_from = "grp", values_from = "value") %>% 
  column_to_rownames("id_site") %>%
  select(`FGFR2-APEX 0` = R2A_R11G_C_A,
         `FGFR2-APEX + FGF10 40` = R2A_R11G_F_A,
         `RAB11-APEX + FGF10 40` = R2M_R11A_F_A
  )
#as.matrix()

# Clustering
set.seed(1)
den_cl_normGFP <- hclust_heatmap(forhm_normGFP, k = 7)

cl <- cbind(
  forhm_normGFP,
  Cluster = den_cl_normGFP
) %>% 
  arrange(Cluster)

# pdf("results/figs/forPaper/STYHeatmap.pdf",
#     width = 40 / 25.4, height = 70 / 25.4, #in. to mm conversion
# )
pheatmap::pheatmap(as.matrix(cl[,1:3]), #use the matrix without the cluster column
                  # color=colorRampPalette(c("navy", "white", "red"))(50),
                   cluster_rows = FALSE,
                   cluster_cols = FALSE,
                   annotation_row = cl[,"Cluster", drop=F], #data frame of clusters above
                   annotation_colors = list(Cluster = RColorBrewer::brewer.pal(max(cl$Cluster), "Set2")),
                   show_rownames = FALSE,
                   fontsize= 6,
                   annotation_legend = F
)
#dev.off() # Will manipulate legend etc. in illustrator

normGFP_cl <- data.frame(id_site = names(den_cl_normGFP), cl = den_cl_normGFP) %>% 
  filter(cl %in% c(4, 5, 6, 7)) %>% 
  mutate(cl_name = ifelse(cl == 4, "RAB11_unique", 
                          ifelse(cl == 5, "RE_profile",
                                 ifelse(cl == 6, "FGFR2_F_unique", 
                                        ifelse(cl == 7, "RE_profile2", NA)))),
         gene = gsub("_.*", "", id_site),
        # R11_interactor = ifelse(gene %in% rab11_interactome$Gene.names, T, F)
        )

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


sty_sig <- filter(sty, id_site %in% lfq_normGFP$id_site) 
sty_sig <- left_join(sty_sig, 
             normGFP_cl[,c("id_site", "cl_name")], 
             by = "id_site") %>% 
  left_join(
    tibble::rownames_to_column(forhm_normGFP, "id_site"),
    by = "id_site"
  )

readr::write_csv(sty_sig, "results/data/APEX_STY_SIG_CLUSTERED.csv")

RE_kegg_forPaper <- calculateEnrichment(
  gsub("_.*", "", 
       normGFP_cl[normGFP_cl$cl_name == "RE_profile", ]$id_site), 
  "KEGG_2019_Human") %>% 
  slice(1:25) %>% 
  arrange(desc(Adjusted.P.value)) %>% 
  mutate(Term = gsub(" Homo sapiens hsa.*", "", Term),
         order = row_number(),
         x = "cluster") %>%
  ggplot(aes(x = x, y = reorder(Term, order))) + 
  #geom_col(fill = "#ABD9E9") + 
  geom_point(aes(color = x, size = Adjusted.P.value)) +
  theme_bw() +
#  xlab("KEGG Pathway") +
 # ylab("-log(FDR)") +
  #coord_flip() +
  guides(size = "none", color = "none") +
  theme(
    legend.key.size = unit(2, "mm"),
    axis.title = element_text(size=6),
    plot.margin = margin(0,0,0,0),
    axis.text.y = element_text(size=5), 		
    axis.text.x = element_text(size = 5, angle = 30, hjust=0, face = "bold"),
    axis.title.y = element_text(size=6, face="bold") ,	
    axis.title.x = element_blank(),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6))
RE_kegg_forPaper
  #       ) + 
  # geom_hline(yintercept = 2.995732, color = "black", 
  #            linetype = "dotted")
# ggsave("results/figs/forPaper/Figure1_STY_RE_pathway.pdf", plot = RE_kegg_forPaper,
#        width = 70, height = 65, unit = "mm", dpi = 300)
ggsave("results/figs/forPaper/STY_ProximalSignalling_KEGG.pdf",
       plot = RE_kegg_forPaper,
       width = 46, height = 75, unit = "mm", dpi = 300)
readr::write_csv(RE_kegg_forPaper$data[,c("Term",
                                          "Adjusted.P.value",
                                          "Genes")],
                 "results/data/SupplementaryTables/4I_APEX-STY-KEGG.csv")


