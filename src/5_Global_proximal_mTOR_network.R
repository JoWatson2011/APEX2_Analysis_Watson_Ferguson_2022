library(tidyverse)
library(ggraph)
library(igraph)
library(enrichR)
library(EnrichmentBrowser)
library(org.Hs.eg.db)

sty_sig <- readr::read_csv("results/data/APEX_STY_SIG_CLUSTERED.csv")
total <- read_csv("results/data/GLOBAL_UPREG_40MIN.csv")
string <- readRDS("data/STRINGexpmtgene_lowconf.rds") %>% 
  dplyr::select(-ID) %>% 
  na.omit()

hsa <- getGenesets(org = "hsa", db = "kegg", cache = TRUE, return.type="list")
mtor_pathways <- c(
  grep("mTOR", names(hsa), value = T),
  grep("utophagy", names(hsa), value = T),
  grep("AMPK", names(hsa), value = T)
)

mtor_pathway_genes <- sapply(mtor_pathways, function(i){
  data.frame(
    pathway = i,
    symbol = unique(
      AnnotationDbi::select(org.Hs.eg.db,
                            hsa[[i]],
                            "SYMBOL",
                            "ENTREZID")$SYMBOL
    )
  )
}, simplify = F, USE.NAMES = T ) %>% 
  bind_rows() %>% 
  mutate(pathway = gsub("hsa[0-9]*_", 
                        "",
                        pathway)
  )

total_v <- total %>% 
  dplyr::select(v= `Gene names`) %>% 
  unique() %>% 
  group_by(v) %>% 
  unique() %>% 
  na.omit() %>% 
  mutate(set = "Global")

RE_v <- sty_sig %>% 
  filter(cl_name == "RE_profile") %>% 
  dplyr::select(v =`Gene names`) %>% 
  na.omit() %>% 
  unique() %>% 
  mutate(set = "Local"
  )

comb_v <- rbind(RE_v, total_v) %>% 
  group_by(v) %>% 
  summarise(set = paste(set, collapse = " & "), .groups = "keep") %>% 
  ungroup() %>% 
  mutate(process = ifelse(v %in% mtor_pathway_genes$symbol,
                          "mTOR",""),
         nodetype = "prot"
  )

el <- string %>% 
  filter(protein1 %in% comb_v$v &
           protein2 %in% comb_v$v) %>% 
  dplyr::select(-experimental) %>% 
  mutate(edgetype = "PPI") %>% 
  rbind(
    data.frame(
      protein1 = sty_sig[sty_sig$`Gene names` %in% comb_v$v &
                           sty_sig$cl_name == "RE_profile",]$id_site,
      protein2 = sty_sig[sty_sig$`Gene names` %in% comb_v$v &
                           sty_sig$cl_name == "RE_profile",]$`Gene names`,
      edgetype = "site-prot"),
    data.frame(
      protein1 = total[total$`Gene names` %in% comb_v$v,]$id,
      protein2 = total[total$`Gene names` %in% comb_v$v,]$`Gene names`,
      edgetype = "site-prot") 
  ) %>% 
  filter(!is.na(protein1) & !is.na(protein2))
v_new <-  data.frame(
  v = unique(el$protein1[!el$protein1 %in% comb_v$v])
) %>% 
  mutate(#name = gsub("^.*_", "", name),
    set = ifelse(v %in% 
                   sty_sig$id_site[sty_sig$cl_name == "RE_profile"] &
                   v %in% total$id,
                 "Local & Global",
                 ifelse(v %in% 
                          sty_sig$id_site[sty_sig$cl_name == "RE_profile"],
                        "Local",
                        ifelse(v %in% total$id,
                               "Global",
                               "")
                 )
    ),
    process = "",
    nodetype = "site"
  ) %>% 
  rbind(comb_v) %>% 
  na.omit()

nw <- BioNet::largestComp(
  graph_from_data_frame(el, 
                        directed = F,
                        vertices = v_new)
)

nw_prots <- graph_from_data_frame(el,
                                  directed = F,
                                  vertices = v_new)
nw_prots <- induced.subgraph(nw_prots, V(nw_prots)[V(nw_prots)$nodetype == "prot"])

ggraph(induced.subgraph(nw, V(nw)[V(nw)$nodetype == "prot"])) + 
  geom_edge_link() + geom_node_point(aes(colour = process)) 

flt_nw <- BioNet::largestComp(
  induced.subgraph(nw, V(nw)[V(nw)$process == "" | 
                               V(nw)$nodetype == "site"])
)

ggraph(flt_nw,
       layout = "nicely") + 
  geom_edge_link() +
  geom_node_label(aes(label = name, fill = set)) 

mtor_nw <- BioNet::largestComp(
  induced.subgraph(nw, 
                   V(nw)[(V(nw)$process == "mTOR") |
                           V(nw)$nodetype == "site"])
)

#Export to cytoscape for visualisation
write_csv(as_data_frame(mtor_nw), 
          "results/data/forCytoscape_mTOR_network_edgelist.csv")
write_csv(as_data_frame(mtor_nw, "vertices"), 
          "results/data/forCytoscape_mTOR_network_vertices.csv")
