library(tidyverse)
library(geneEnrichment)
library(VennDiagram)

pg_sig <- readr::read_csv("results/data/APEX_PRO_SIG.csv")
sty_sig <- readr::read_csv("results/data/APEX_STY_SIG_CLUSTERED.csv") %>% 
  mutate(`Gene names` = gsub("_.*","", id_site))
total <- read_csv("results/data/GLOBAL_UPREG_40MIN.csv")
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
  mutate(cl_name = ifelse(is.na(cl_name), "Other regulated sites", "RE_profile")) %>% 
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


KEGG <- calculateEnrichment(unique(phos_prots$`Gene names`), 
                            enrichr_db = "KEGG_2019_Human") %>% 
  slice(1:15) %>% 
  mutate(group = "Phosphorylated protein")
KEGG_STYonly <- calculateEnrichment(
  unique(
    filter(sty_sig, !is.na(cl_name))$`Gene names`
  ), 
  enrichr_db = "KEGG_2019_Human") %>% 
  slice(1:15) %>% 
  mutate(group = "FGFR2b Recycling Proximal Signalling")
KEGG_global <- calculateEnrichment(
  unique(
    total$`Gene names`
  ), 
  enrichr_db = "KEGG_2019_Human") %>% 
  slice(1:15) %>% 
  mutate(group = "Global")
RE_kegg_forPaper <- rbind(KEGG, KEGG_STYonly, KEGG_global) %>% 
  select(Term,Adjusted.P.value, Genes, group) %>%  
  # arrange(desc(Adjusted.P.value)) %>% 
  group_by(Term) %>%
  mutate(n = n()) %>% 
  ungroup() %>% 
  arrange(n) %>%
  mutate(order = 1:n()) %>% 
  mutate(Term = gsub(" Homo sapiens hsa.*", "", Term)) %>% 
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
  ))  +
  scale_size(breaks = seq(-0.05, 0, 0.01),
             limits = c(-0.05, 0),
             range = c(1,3)) +
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
ggsave("KEGG_enrichment.pdf", RE_kegg_forPaper,
       width = 206.453, height = 100, unit = "mm")


