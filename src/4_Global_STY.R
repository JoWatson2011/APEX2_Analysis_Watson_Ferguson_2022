set.seed(1)
library(tidyverse)
library(patchwork)
library(paletteer)
library(enrichR)
library(geneEnrichment)
source("src/functions/createHeatmap.R")

######
# IMPORT &
# EXTRACT COLNAMES
######
sty <-
  readRDS("data/sty_Flt.rds") %>% 
  mutate(`Gene names` = gsub(";.*", "", `Gene names`),
         id_site = paste0(`Gene names`, "_", `Amino acid`, Position))
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

lfq <- sty %>%
  dplyr::select(id_site, grep("T_0[123]", colnames(sty))) %>%
  pivot_longer(cols = -1) %>%
  mutate(grp = gsub("_0[123]", "", name))

# Do the bait proteins have an effect on quantification of global phosphorylation?
#   NULL HYPOTHESIS: Quantification of sites in the 0 min and 40 min conditions should not change depending on bait protein.
# We want to ACCEPT the null hypothesis (p > 0.05)

total_sty <- sty %>%
  mutate(id = paste0(`Gene names`, "_", `Amino acid`, Position)) %>% 
  select("id", grep("T_0", colnames(.))) 


total_sty$C.p <- apply(total_sty[,2:16], 1, function(i){
  GFPC <- c(i["R2M_GFPA_C_T_01"],
            i["R2M_GFPA_C_T_02"],
            i["R2M_GFPA_C_T_03"]
  )
  if(class(GFPC)== "list"){
    GFPC <- do.call(c, GFPC)
  }
  R2AC <- c(i["R2A_R11G_C_T_01"],
            i["R2A_R11G_C_T_02"],
            i["R2A_R11G_C_T_03"]
  )
  if(class(R2AC)== "list"){
    R2AC <- do.call(c, R2AC)
  }
  
  t <- try(t.test(GFPC, R2AC), silent = T)
  
  if(class(t) == "try-error"){
    return(NA)
  }else{
    return(t$p.value)
  }
})

total_sty$C.p <- p.adjust(total_sty$C.p)
boxplot(total_sty$C.p)

total_sty$F.p <- apply(total_sty[,2:16], 1, function(i){
  dat <- data.frame(cond = c(rep("GFPF", 3),
                             rep("R2AF", 3),
                             rep("R11F", 3)
  ),
  intensity = i[c(
    grep("GFPA_F", names(i)),
    grep("R2A_R11G_F", names(i)),
    grep("R11A_F", names(i))
  )]
  )
  
  oneway <- aov(intensity ~ cond, data = dat)
  p <- summary(oneway)[[1]][1, 5]
  
  
  return(p)
  
})

total_sty$F.p <- p.adjust(total_sty$F.p)

n <- pivot_longer(total_sty[,c("id", "C.p", "F.p")],cols = -id) %>% 
  group_by(name) %>%  
  filter(value < 0.05) %>% 
  summarise(n = n())

plot_p <- pivot_longer(total_sty[,c("id", "C.p", "F.p")],
                       cols = -id) %>% 
  ggplot(aes(x = name, y = -log(value))) +
  geom_boxplot(outlier.size = 0.2) +
  #geom_jitter(size = 0.2, alpha = 0.7) +
  #geom_violin(aes(fill = name), show.legend = F)+
  geom_hline(yintercept = -log(0.05)) +
  ggrepel::geom_text_repel(data = n, 
                           aes(x = name, label = n), 
                           y = -log(0.02), size= 2,
                           direction = "x", color = "red") +
  scale_x_discrete(labels = c("GFP-APEX UT\n~ FGRR2-APEX UT`,\nt-test", 
                              "GFP-APEX_40`\n~ FGFR2-APEX_40`\n~ RAB11-APEX_40`,\nOne way Anova")) +
  ylab("-log(FDR)") +
  #ggtitle("Do Global samples from UT and FGF10 40` differ due to APEX tags?") +
  theme_bw() +
  theme(text = element_text(size =6),
        axis.title.x = element_blank())

which(total_sty$C.p <= 0.05)
sum(total_sty$F.p <= 0.05)

ggsave("results/figs/forPaper/global_APEX2_effect.pdf", plot_p,
       width = 50, height = 50,
       units = "mm", dpi = "print")

plot_p$data %>% filter(value > 0.05) 

## Differences between time points
total_sty <- sty %>% 
  mutate(id = paste0(`Gene names`, "_", `Amino acid`, Position)) %>% 
  select(id, `Gene names`, grep("_T_0", experiments, value = T)) %>% 
  pivot_longer(cols = -c(id, `Gene names`),
               names_to = "experiment",
               values_to = "intensity") %>% 
  mutate(experiment = ifelse(grepl("_C_", experiment), "Control", "FGF10")) %>% 
  group_by(id, experiment) %>% 
  mutate(experiment = if(any(experiment == "Control")){
    paste0(experiment, "_", 1:6)
  }else{
    paste0(experiment, "_", 1:9)
  }) %>% 
  pivot_wider(names_from = experiment, values_from = intensity, values_fn = mean) 

write_csv(total_sty, "results/data/STY_Total_sum.csv")

# Upload to Perseus

total_sty_volcano <- readr::read_tsv("results/data/perseus_Control-FGF10_totalsty.txt")
total_sty_sig <- filter(total_sty_volcano, Significant == "+")

g <- total_sty_volcano %>% 
  mutate(`Gene names`= gsub(";.*", "", `Gene names`)) %>% 
  mutate(label = ifelse(`Gene names` == "FGFR2",
                        "FGFR2",
                        NA)) %>% 
  ggplot(aes(y = `-Log(P-value)`, x = Difference)) +
  geom_point(size = 0.1) +
  gghighlight::gghighlight(Significant == "+",
                           label_params = list(
                             size = 5 / (14/5),
                             label.padding = 0.1,
                             box.padding = 0.1),
                           use_direct_label = F
  ) +
  geom_point(aes(color = label, size =label)) +
  scale_color_manual(values = c(FGFR2 = "#ABD9E9"), na.value = "black", guide = "none") +
  scale_size_manual(values = c(FGFR2 = 0.5), na.value = 0.1, guide = "none") +
  xlab("Difference (Control - FGF10)") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(), 
        axis.line = element_line("black"),
        # Text sizes may need modifying based on fig. sizes
        axis.text = element_text(size=5), 		
        axis.title = element_text(size=5, face="bold"),
        plot.margin = margin(0,0,0,0,"cm")
  )
ggsave("results/figs/forPaper/Figure2_D_TotalVolcano_FGFRlabelled.pdf", 
       g, width = 4.5, height = 5, unit = "cm", dpi = 300)


export <- total_sty_sig %>% 
  filter(Difference < 0 ) %>% 
  arrange(desc(`-Log(P-value)`)) 

library(geneEnrichment)

kegg <- calculateEnrichment(unique(export$`Gene names`),"KEGG_2019_Human", visualise = T)
react <- calculateEnrichment(unique(export$`Gene names`),"Reactome_2016", visualise = T)

kegg_forPaper <- calculateEnrichment(unique(export$`Gene names`),"KEGG_2019_Human") %>% 
  arrange(desc(Adjusted.P.value)) %>% 
  slice(1:25) %>% 
  mutate(Term = gsub(" Homo sapiens hsa.*", "", Term),
         order = row_number()) %>%
  ggplot(aes(x = reorder(Term, order), y = -log(Adjusted.P.value))) + 
  geom_col(fill = "#ABD9E9") + 
  theme_bw() +
  xlab("KEGG Pathway") +
  ylab("-log(FDR)") +
  coord_flip() +
  theme(axis.text = element_text(size=5), 		
        axis.title = element_text(size=5, face="bold")) + 
  geom_hline(yintercept = 2.995732, color = "black")
ggsave("results/figs/forPaper/Total-Up_pathway.pdf", plot = kegg_forPaper,
       width = 80, height = 80, unit = "mm", dpi = 300)

readr::write_csv(export, "results/data/GLOBAL_UPREG_40MIN.csv")

apex_sty_sig <- readr::read_csv("results/data/APEX_STY_SIG_CLUSTERED.csv")

venn <- ggvenn::ggvenn(list(Total = total_sty_sig$id, 
                            `APEX:\nFGFR2\nrecycling profile` = 
                              apex_sty_sig[apex_sty_sig$cl_name == "RE_profile",]$id_site),
                       fill_color = c("#ABD9E9", "#D73027"),
                       text_size = 2,
                       set_name_size = 1,
                       show_percentage = F)
ggsave("results/figs/forPaper/Global_Proximal-overlap.pdf",
       width = 50, height = 50, units = "mm", dpi = 300)

total_apex_RE <- apex_sty_sig[apex_sty_sig$cl_name == "RE_profile",] %>% 
  filter(id_site %in% total_sty_sig$id)

justTotal <- unique(export$`Gene names`)
total_apex_RE <- unique((apex_sty_sig[apex_sty_sig$cl_name == "RE_profile",] %>% 
                           filter(id_site %in% total_sty_sig$id))$`Gene names`)
justApex <- unique(filter(apex_sty_sig, !(id_site %in% total_sty_sig$id) &
                            cl_name == "RE_profile")$`Gene names`)
genes <- list(Global = justTotal,
              `Local` = justApex)
tmp <- lapply(names(genes), function(i){
  calculateEnrichment(genes[[i]],
                      "KEGG_2019_Human",
                      visualise = F) %>% 
    arrange(desc(Adjusted.P.value)) %>% 
    slice(1:25) %>% 
    mutate(group = i)
}) %>% 
  bind_rows() %>% 
  ggplot(aes(x= group, y= Term)) +
  geom_point(aes(color = group, size = -Adjusted.P.value)) +
  scale_color_manual("", values = c("#4575B4", "#D73027")) +
  #  geom_tile(aes(fill = Adjusted.P.value)) +
  scale_x_discrete(position = "top") +
  scale_size("FDR", range =c(1,3)) +
  ylab("KEGG Pathway") +
  guides(color = "none") +
  theme_bw() +
  theme(
    # plot.background = element_blank(),
    #  panel.background = element_rect(fill = "black"),
    #panel.grid = element_line("grey", linetype="dashed"), 
    # panel.grid = element_blank(),
    # legend.key = element_blank(),
    #   axis.line = element_line("black"),
    # Text sizes may need modifying based on fig. sizes
    plot.margin = margin(0,0,0,0,"mm"),
    axis.text.y = element_text(size=6), 		
    axis.text.x = element_text(size = 6, angle = 30, hjust=0, face = "bold"),
    axis.title.y = element_text(size=6, face="bold") ,	
    axis.title.x = element_blank(),
    legend.key.size = unit(0.25, "mm"),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6)
  )
ggsave("results/figs/forPaper/Global_proximal_pathwayHeatmap.pdf",
       tmp, width = 90, height = 106, unit = "mm", dpi = 300)
readr::write_csv(mutate(tmp$data,
                        group = ifelse(group == "Local", "Proximal", "Global")),
                 "results/data/SupplementaryTables/Global-STY-KEGG.csv")
