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

experiments <-
  grep("Intensity ", colnames(sty), value = T) %>%
  gsub("Intensity ", "", .) %>%
  gsub("_40", "", .) %>%
  gsub("_Ph_", "_", .)

colnames(sty)[grep("Intensity ", colnames(sty))] <-
  experiments

# Impute data
lfq <- sty %>% 
  dplyr::select(id_site, `Gene names`, grep("T_0[123]", colnames(sty))) %>%
  pivot_longer(cols = -c(id_site, `Gene names`)) %>%
  mutate(grp = gsub("_0[123]", "", name)) %>% 
  group_by(id_site,grp) %>% 
  filter(sum(is.na(value))<2) %>% 
  ungroup() %>% 
  select(-grp) %>% 
  pivot_wider(values_fn=median)
lfq[,3:ncol(lfq)] <-  
  apply(lfq[,3:ncol(lfq)], 2, function(i) 
    imputeLCMD::impute.QRILC(as.matrix(i), tune.sigma = 1.9)[[1]]
  )
total_sty <- as.data.frame(lfq)


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


total_sty$C.p <- apply(total_sty[,3:17], 1, function(i){
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

total_sty$F.p <- apply(total_sty[,3:17], 1, function(i){
  dat <- data.frame(
    cond = factor(
      c(rep("GFPF", 3),
        rep("R2AF", 3),
        rep("R11F", 3))
    ),
    intensity = as.numeric(i[c(
      grep("GFPA_F", names(i)),
      grep("R2A_R11G_F", names(i)),
      grep("R11A_F", names(i))
    )])
  )
  
  oneway <- aov(intensity ~ cond, data = dat)
  p <- summary(oneway)[[1]][1, 5]
  
  
  return(p)
  
})

total_sty$F.p <- p.adjust(total_sty$F.p)

n <- pivot_longer(total_sty[,c("id_site", "C.p", "F.p")],cols = -id_site) %>% 
  group_by(name) %>%  
  filter(value < 0.05) %>% 
  summarise(n = n())

plot_p <- pivot_longer(total_sty[,c("id_site", "C.p", "F.p")],
                       cols = -id_site) %>% 
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


## Differences between time points
total_sty <- total_sty %>% 
  pivot_longer(cols = -c(id_site, `Gene names`),
               names_to = "experiment",
               values_to = "intensity") %>% 
  mutate(experiment = ifelse(grepl("_C_", experiment), "Control", "FGF10")) %>% 
  group_by(id_site, experiment) %>% 
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


readr::write_csv(export, "results/data/GLOBAL_UPREG_40MIN.csv")

apex_sty_sig <- readr::read_csv("results/data/APEX_STY_SIG_CLUSTERED.csv")

venn <- ggvenn::ggvenn(list(Total = total_sty_sig$id_site, 
                            `APEX:\nFGFR2\nrecycling profile` = 
                              apex_sty_sig[!is.na(apex_sty_sig$cl_name),]$id_site),
                       fill_color = c("#ABD9E9", "#D73027"),
                       text_size = 2,
                       set_name_size = 1,
                       show_percentage = F)
ggsave("results/figs/forPaper/Global_Proximal-overlap.pdf",
       width = 50, height = 50, units = "mm", dpi = 300)
