library(tidyverse)
library(geneEnrichment)
library(VennDiagram)

pg_sig <- readr::read_csv("results/data/APEX_PRO_SIG_CLUSTERED.csv")
sty_sig <- readr::read_csv("results/data/APEX_STY_SIG_CLUSTERED.csv") %>% 
  filter(cl_name != "")

re_cluster <- readr::read_csv("results/data/SupplementaryTables/4I_APEX-STY-KEGG.csv")

phos_prots <- sty_sig[sty_sig$`Gene names` %in%  pg_sig$`Gene names`,]

########
venn <- ggvenn::ggvenn(list("Protein Groups" = na.omit(pg_sig$`Gene names`),
                  "Phosphorylated\nproteins" = na.omit(sty_sig$`Gene names`)),
                  fill_color = c("#D73027", "#4575B4")
) + theme(text = element_text(size = 5))
ggsave("results/figs/ForPaper/STY-PRO_Venn.pdf", venn,
       width = 61, height = 31, unit = "mm", dpi = 300)

pie <- sty_sig %>%
  filter(id_site %in% phos_prots$id_site) %>% 
  mutate(cl_name = ifelse(cl_name != "RE_profile", "Other regulated sites", "RE_profile")) %>% 
  group_by(cl_name) %>% 
  summarise(n = n()) %>%
  select(name = cl_name, n) %>% 
  mutate(
    Fraction = n / sum(n), # Compute percentages
    ymax = cumsum(Fraction), # Compute the cumulative percentages (top of each rectangle)
    ymin = c(0, head(ymax, n=-1)), # Compute the bottom of each rectangle
    labelPosition = (ymax + ymin) / 2, # Compute label position
    Percentage = Fraction * 100, # Compute a good label
    Percentage = format(round(Percentage, 1), nsmall = 1),
    label = paste0(name, "\n (", Percentage, "%)")
  ) %>% 
  ggplot(aes(ymax=ymax, ymin=ymin, xmax=1.2, xmin=0.2, fill=name , show.legend = T)) +
  geom_rect(alpha = 0.6) +
  geom_text( x=3, aes(y=labelPosition, label=label, color=name), size=2, check_overlap = F, nudge_x = -0.2) + # x here controls label position (inner / outer)
  scale_fill_manual(values = c("#ABD9E9", "#4575B4")) +
  scale_color_manual(values = c("#ABD9E9", "#4575B4")) +
  coord_polar(theta="y") +
  xlim(c(-2, 3)) +
  annotate(geom = 'text', x = -1, y = 0, label = paste(nrow(phos_prots), "\nsites")) +
  theme_void() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)#,
        #plot.margin = margin(0.1,0.1,0.1,0.1, "mm")
  )

ggsave("results/figs/forPaper/Figure2_C_PhosProts_clusteredSites.pdf",
       pie,
       width = 50,
       height = 50,
       unit = "mm",
       dpi = 300)

############

KEGG <- calculateEnrichment(unique(phos_prots$`Gene names`), 
                    enrichr_db = "KEGG_2019_Human", visualise = T)
reactome <- calculateEnrichment(unique(phos_prots$`Gene names`), 
                    enrichr_db = "Reactome_2016", visualise = T, simplify = T)
gocc <- calculateEnrichment(unique(phos_prots$`Gene names`), 
                                enrichr_db = "GO_Cellular_Component_2018", visualise = T)

saveEnrichmentVis(KEGG, "results/figs/APEX_PRO_STY/KEGG.tiff")
saveEnrichmentVis(reactome, "results/figs/APEX_PRO_STY/Reactome.tiff")

#### FOR PAPER
RE_kegg_forPaper <- calculateEnrichment(unique(phos_prots$`Gene names`), 
                                         enrichr_db = "KEGG_2019_Human") %>% 
  select(Term,Adjusted.P.value, Genes) %>% 
  mutate(group = "Phosphorylated protein") %>% 
  slice(1:25) %>% 
  rbind(mutate(re_cluster,group = "FGFR2b recycling")) %>% 
  arrange(desc(Adjusted.P.value)) %>% 
  mutate(Term = gsub(" Homo sapiens hsa.*", "", Term),
         order = row_number()) %>%
  ggplot(aes(x = group, y = reorder(Term, order))) + 
  #geom_col(fill = "#ABD9E9") + 
  geom_point(aes(color = group, size = -Adjusted.P.value)) +
  theme_bw() +
  #  xlab("KEGG Pathway") +
  # ylab("-log(FDR)") +
  coord_flip() +
  guides(color = "none") +
  scale_color_manual(values = c(
    "Phosphorylated protein" = "#D73027",
    "FGFR2b recycling" = "#ABD9E9"
  )) +
  scale_size(breaks = seq(-0.05, 0, 0.01),
             limits = c(-0.05, 0)) +
  theme(
    legend.key.size = unit(2, "mm"),
    axis.title = element_text(size=6),
    plot.margin = margin(0,0,0,0),
    axis.text.y = element_text(size=6), 		
    axis.text.x = element_text(size = 6, angle = 30, hjust=1, face = "bold"),
    axis.title.y = element_text(size=6, face="bold") ,	
    axis.title.x = element_blank(),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6))
ggsave("results/figs/forPaper/Figure2_B_PhosProts_pathway.pdf",
       plot = RE_kegg_forPaper, height = 70,
       width = 210, unit = "mm")
  # ggplot(aes(x = reorder(Term, order), y = -log(Adjusted.P.value))) + 
  # geom_col(fill = "#D73027") + 
  # theme_bw() +
  # labs(x = "KEGG Pathway",
  #      y = "-log(FDR)") +
  # coord_flip() +
  # theme(axis.text = element_text(size=5), 		
  #       axis.title = element_text(size=5, face="bold"),
  #       panel.grid = element_blank(),
  #       plot.caption = element_text(size = 6, face = "italic")) + 
  # geom_hline(yintercept = 2.995732, color = "black")

readr::write_csv(RE_kegg_forPaper$data[,c("Term",
                                          "Adjusted.P.value",
                                          "Genes")],
                 "results/data/SupplementaryTables/4L_APEX-PRO-STY-KEGG.csv")

####




mTOR_pw_edges <- readr::read_tsv("results/data/forCytoscape_mTORpw_edges.txt")
mTOR_pw_nodes <- readr::read_tsv("results/data/forCytoscape_mTORpw_nodes.txt") %>% 
  select(name) %>% 
  mutate(type = "protein")


mTOR_pw_prots <- pg_sig[pg_sig$`Gene names` %in% mTOR_pw_nodes$name,c("Gene names", "cl_name")]
mTOR_pw_phos <- sty_sig[sty_sig$`Gene names` %in% mTOR_pw_nodes$name,c("Gene names", "id_site", "cl_name")]

mTOR_pw_edges <- rbind(mTOR_pw_edges,
                       data.frame(from = mTOR_pw_phos$`Gene names`,
                                  to = mTOR_pw_phos$id_site,
                                  action = "local_phosphorylated"),
                       data.frame(from = "RRAGC",
                                  to = "MTOR",
                                  action = "binding/association"),
                       data.frame(from = "TSC1",
                                  to = c("RHEB", "TSC2"),
                                  action = c("inhibition", "binding/association")),
                       data.frame(from = c("GSK3B", "IKBKB", "RPS6KA3", "MAPK1", "MAPK3",
                                           "AKT1", "AKT2", "AKT3"),
                                  to = "TSC1",
                                  action = c("activation phosphorylation",
                                             rep("inhibition phosphorylation", 7))
                                  )
)



mTOR_pw_nodes <- lapply(unique(c(mTOR_pw_edges$from, mTOR_pw_edges$to)),
                        function(i){
                          if(!(i %in% mTOR_pw_nodes$name)){
                            if(grepl("_[STY]", i)){
                              return(c(name = i,
                                       type = "phospho"))
                            }else{
                              return(c(name = i,
                                       type = "protein")
                              )
                            }
                          }
                        }
) %>% do.call(rbind, .) %>% 
  rbind(mTOR_pw_nodes)

mTOR_pw_nodes$cl<- ifelse(mTOR_pw_nodes$name %in% mTOR_pw_prots$`Gene names`, 
                          mTOR_pw_prots$cl_name, 
                          ifelse(mTOR_pw_nodes$name %in% mTOR_pw_phos$id_site,
                                 mTOR_pw_phos$cl_name, 
                                 NA)
)

readr::write_tsv(mTOR_pw_edges, "results/data/forCytoscape_mTORpw_edges.txt")
readr::write_tsv(mTOR_pw_nodes, "results/data/forCytoscape_mTORpw_nodes.txt")



