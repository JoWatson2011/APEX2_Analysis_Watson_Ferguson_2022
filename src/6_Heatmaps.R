library(tidyverse)
library(patchwork)

glob <- read_csv("results/data/GLOBAL_UPREG_40MIN.csv")
prox <- read_csv("results/data/SupplementaryTables/APEX_STY_data.csv")
proxpro <- read_csv("results/data/SupplementaryTables/APEX_Proteome_data.csv")

targets <- c(
  "FGFR",
  "PLC",
  "SHC",
  "EEA1",
  "LAMP1",
  "RAB11FIP1",
  "CHMP1",
  "RAB7",
  "RAB25",
  "TFRC",
  "VPS5"
)
sty <- lapply(c("PLCG1",
                "RAB11FIP1",
                "SHC1_Y239",
                "FGFR2_Y657",
                "FGFR2_Y656",
                "RAB25",
                "RAB7",
                "EEA1",
                "CHMP1",
                "TFRC",
                "VPS5"),
              function(gene){
                filter(prox, grepl(gene, id_site))
              }) %>%
  bind_rows() %>% 
  select(id_site, matches("_medianNormalised")) %>% 
  pivot_longer(-id_site) %>% 
  mutate(name = gsub("_medianNormalised", "", name)) %>% 
  filter(!is.na(value))
  # select(id_site, matches("Proximal.*R[123]")) %>% 
  # pivot_longer(-id_site) %>% 
  # mutate(name = gsub("_R[123]", "", name)) %>% 
  # group_by(id_site, name) %>% 
  # summarise(value = median(value, na.rm = T), .groups = "keep") %>% 
  # group_by(id_site) %>% 
  # mutate(value = (value - mean(value, na.rm = T))/sd(value, na.rm=T)) %>% 
  # ungroup() %>% 
  # filter(!is.na(value))
sty_g <- ggplot(sty, aes(x=name, y= id_site, fill=value)) + 
  geom_tile(show.legend = F) +
  scale_x_discrete(limits = c(#"GFP-APEX_UT_Proximal",
                              "FGFR2-APEX 0_Proximal",
                             # "GFP-APEX_FGF10_Proximal",
                              "FGFR2-APEX + FGF10 40_Proximal",
                              "RAB11-APEX + FGF10 40_Proximal"),
                   labels = c(
                    # "GFP-APEX UT",
                     "FGFR2-APEX UT",
                  #   "GFP-APEX + FGF10 40",
                     "FGFR2-APEX + FGF10 40",
                     "RAB11-APEX + FGF10 40"
                   )
  ) +
  theme_minimal() +
  ylab("Phosphorylated site")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        #panel.background = element_rect(fill = "black"),
        text = element_text(size = 8))+
  scale_fill_gradientn(colors = c("darkblue", "white", "darkred"),
                       limits = c(-1.2,1.2))
sty_g  
pro <- lapply(targets,
              function(gene){
                filter(proxpro, grepl(gene, `Gene names`)
                       )
              }) %>%
  bind_rows() %>% 
  select(`Gene names`, ends_with("Normalised")) %>% 
  pivot_longer(-`Gene names`) %>% 
  filter(!grepl("PLC", `Gene names`)) %>% 
  filter(!is.na(value))
pro_g <- ggplot(pro, aes(x=name, y= `Gene names`, fill=value)) + 
  geom_tile() + 
  scale_x_discrete(limits = c("FGFR2-APEX_UT_Proximal_medianNormalised",
                              "FGFR2-APEX_FGF10_Proximal_medianNormalised",
                              "RAB11-APEX_FGF10_Proximal_medianNormalised"),
                   labels = c(
                     "FGFR2-APEX UT",
                     "FGFR2-APEX + FGF10 40",
                     "RAB11-APEX + FGF10 40"
                   )
                   ) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        text = element_text(size = 8)) +
  ylab("Protein") +
  scale_fill_gradientn(colors = c("darkblue", "white", "darkred"),
                       name = "Intensity normalised\nto GFP-APEX2",
                       limits = c(-1.2,1.2))

ggsave("results/figs/endosomalMarkers.pdf", 
       sty_g + pro_g + plot_layout(guides = "collect"))




mtor <- read_csv("results/data/cytoscape/Total_RE_vertices_withSites_mtor.csv")


sty_mtor <- prox %>% 
  filter(id_site %in% mtor$name & 
         (!is.na(Proximal_Regulated_Cluster)|
           !is.na(Global_Regulated))) %>% 
  bind_rows() %>% 
  select(id_site, matches("_medianNormalised")) %>% 
  pivot_longer(-id_site) %>% 
  mutate(name = gsub("_medianNormalised", "", name)) %>% 
  filter(!is.na(value))
# select(id_site, matches("Proximal.*R[123]")) %>% 
# pivot_longer(-id_site) %>% 
# mutate(name = gsub("_R[123]", "", name)) %>% 
# group_by(id_site, name) %>% 
# summarise(value = median(value, na.rm = T), .groups = "keep") %>% 
# group_by(id_site) %>% 
# mutate(value = (value - mean(value, na.rm = T))/sd(value, na.rm=T)) %>% 
# ungroup() %>% 
# filter(!is.na(value))
sty_mtor_g <- ggplot(sty_mtor, aes(x=name, y= id_site, fill=value)) + 
  geom_tile(show.legend = F) +
  scale_x_discrete(limits = c(#"GFP-APEX_UT_Proximal",
    "FGFR2-APEX 0_Proximal",
    # "GFP-APEX_FGF10_Proximal",
    "FGFR2-APEX + FGF10 40_Proximal",
    "RAB11-APEX + FGF10 40_Proximal"),
    labels = c(
      # "GFP-APEX UT",
      "FGFR2-APEX UT",
      #   "GFP-APEX + FGF10 40",
      "FGFR2-APEX + FGF10 40",
      "RAB11-APEX + FGF10 40"
    )
  ) +
  theme_minimal() +
  ylab("Phosphorylated site")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        #panel.background = element_rect(fill = "black"),
        text = element_text(size = 8))+
  scale_fill_gradientn(colors = c("darkblue", "white", "darkred"),
                       limits = c(-1.2,1.2))
pro_mtor <- proxpro %>% 
  filter(`Gene names` %in% mtor$name & 
           !is.na(Regulated)) %>% 
  bind_rows() %>% 
  select(`Gene names`, ends_with("Normalised")) %>% 
  pivot_longer(-`Gene names`) %>% 
  mutate(name = gsub("_medianNormalised", "", name)) %>% 
  filter(!is.na(value))
# select(id_site, matches("Proximal.*R[123]")) %>% 
# pivot_longer(-id_site) %>% 
# mutate(name = gsub("_R[123]", "", name)) %>% 
# group_by(id_site, name) %>% 
# summarise(value = median(value, na.rm = T), .groups = "keep") %>% 
# group_by(id_site) %>% 
# mutate(value = (value - mean(value, na.rm = T))/sd(value, na.rm=T)) %>% 
# ungroup() %>% 
# filter(!is.na(value))
pro_mtor_g <- ggplot(pro_mtor, aes(x=name, y= `Gene names`, fill=value)) + 
  geom_tile() +
  scale_x_discrete(limits = c(#"GFP-APEX_UT_Proximal",
    "FGFR2-APEX_UT_Proximal",
    # "GFP-APEX_FGF10_Proximal",
    "FGFR2-APEX_FGF10_Proximal",
    "RAB11-APEX_FGF10_Proximal"),
    labels = c(
      # "GFP-APEX UT",
      "FGFR2-APEX UT",
      #   "GFP-APEX + FGF10 40",
      "FGFR2-APEX + FGF10 40",
      "RAB11-APEX + FGF10 40"
    )
  ) +
  theme_minimal() +
    ylab("Protein")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        #panel.background = element_rect(fill = "black"),
        text = element_text(size = 8))+
  scale_fill_gradientn(colors = c("darkblue", "white", "darkred"),
                       name = "Intensity normalised\nto GFP-APEX2",
                       limits = c(-1.2,1.2))


ggsave("results/figs/revisions_mTORquant.tiff",
       sty_mtor_g + pro_mtor_g + plot_layout(guides = "collect") )


prox %>% 
  filter(id_site %in% mtor$name) %>% 
  select(id_site, matches("^[FR].*_R[123]")) %>% 
  pivot_longer(-id_site) %>% 
  mutate(name = gsub("_R[123]", "", name)) %>% 
  group_by(id_site, name) %>% 
  summarise(value = median(value, na.rm = T), .groups = "keep") %>% 
  group_by(id_site) %>% 
  mutate(value = scale(value)) %>% 
  mutate(value = (value - mean(value, na.rm = T))/sd(value, na.rm=T)) %>% 
  ungroup() %>% 
  ggplot(aes(x= name, y = id_site, fill=value)) +
  geom_tile() +
  scale_fill_gradient2() +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
