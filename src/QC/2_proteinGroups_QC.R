library(dplyr)
library(tidyr)
library(gplots)
library(data.table)
library(patchwork)
# library(readr)
library(RColorBrewer)
library(ggplot2)


# Read proteinGroups
proteinGroups <- readRDS("data/proteinGroups.RDS")
experiment_cols <- grep("LFQ intensity .*_P_",
                        colnames(proteinGroups),
                        value = T)
# Convert Intensity columns to 'numeric'
proteinGroups[,experiment_cols] <- apply(proteinGroups[,experiment_cols],
                                         2,
                                         as.numeric)

# Summary of variables proteinGroups will be filtered on.
# After running the script this information will be stored in the
# variable named report.
report <- vector()
report["Total Proteins Identified"] <- nrow(proteinGroups)
report["Potential Contaminants"] <- sum(proteinGroups$`Potential contaminant` == "+")
report["Reverse"] <- sum(proteinGroups$Reverse == "+")
report["Only Identified by Site"] <- sum(proteinGroups$`Only identified by site` == "+")
report["No Quantitative Data (incomplete cases)"] <- proteinGroups %>%
  select(grep("LFQ intensity" , experiment_cols, value = T)) %>%
  filter(rowSums(. == 0) == ncol(.)) %>% nrow()


# Filter proteinGroups.txt to remove contaminants, reverse.
proteinGroups_flt <- proteinGroups %>%
  filter(`Potential contaminant` != "+",
         `Reverse` != "+",
         `Only identified by site` != "+",
         `Razor + unique peptides` > 1,
         `Unique + razor sequence coverage [%]` >= 5 
         )

# Add info to report
report["Proteins remaining following filtering"] <- nrow(proteinGroups_flt)
report["Incomplete cases in filtered dataset"] <- proteinGroups_flt %>%
  select(all_of(experiment_cols)) %>%
  filter(rowSums(. == 0) == ncol(.)) %>%
  nrow()

# Filter rows with no quantitative information
proteinGroups_flt_cc <- proteinGroups_flt[rowSums(proteinGroups_flt[experiment_cols] > 0) >= 1,]

#Collapse gene/protein names
proteinGroups_flt_cc <- proteinGroups_flt_cc %>% 
  mutate(`Gene names` = gsub(";.*", "", `Gene names`),
         `Majority protein IDs` = gsub(";.*", "", `Majority protein IDs`),
         `Protein names` = gsub(";.*", "", `Protein names`))

# Transform
proteinGroups_log <- apply(proteinGroups_flt_cc[,experiment_cols], 2, function(i){
  tmp <- ifelse(i == 0, NA, i)
  tmp2 <- log10(tmp)
  # tmp3 <- limma::normalizeBetweenArrays(tmp2, "cyclicloess")
  return(tmp2)
})

# Normalise
proteinGroups_norm <- cbind(
  proteinGroups_flt_cc[,!(names(proteinGroups_flt_cc) %in% c(experiment_cols))],
  limma::normalizeBetweenArrays(proteinGroups_log, "quantile")
)

# Visualise normalilsation
proIntUnnorm <- proteinGroups_flt_cc %>%
  select(`id`, all_of(experiment_cols)) %>% 
  pivot_longer(cols = -c(`id`),
               names_to = "experiment",
               values_to = "intensity") %>% 
  mutate(
    intensity = ifelse(intensity == 0, NA, intensity), 
    intensity = log10(intensity),
    rep = substr(experiment,
                 nchar(experiment)-1,
                 nchar(experiment))) %>% 
  ggplot(aes(x = intensity, color = rep, group = experiment)) +
  geom_density() +
  theme(legend.position = 'none')

proIntNorm <- proteinGroups_norm %>% 
  select(`id`, all_of(experiment_cols)) %>% 
  pivot_longer(cols = -c(`id`),
               names_to = "experiment",
               values_to = "intensity") %>% 
  mutate(rep = substr(experiment,
                      nchar(experiment)-1,
                      nchar(experiment))) %>% 
  filter(intensity < 10) %>% 
  ggplot(aes(x = intensity, color = rep, group = experiment)) +
  geom_density() +
  theme(legend.position = 'none') +
  ggtitle("", subtitle = "Normalised Intensities")

proIntUnnorm / proIntNorm

ggsave("results/figs/QC/proNormCurve.tiff",
      proIntUnnorm / proIntNorm,
      width = 210, height = 297, units = "mm")

saveRDS(proteinGroups_norm, "data/proteinGroups_Flt.rds")
write.csv(proteinGroups_norm, "data/proteinGroups_Flt.csv")

#### FIGURES
# Pie chart to visualise propoprtion of missing values across all experiments.
pie(table(proteinGroups_flt_cc[,experiment_cols] == 0 |
            is.na(proteinGroups_flt_cc[,experiment_cols])),
    labels = c("Not NA", "NA"),
    main = "Missing Ratios in Filtered Data")

#Heatmap correlation
cors <- proteinGroups %>%
  select(c(grep("_A_", colnames(.)),
           grep("_T_", colnames(.)))
         ) %>% 
  cor(use = "pairwise.complete.obs")
rownames(cors) <- colnames(cors) <- gsub("LFQ intensity ",
                                         "",
                                         colnames(cors))
cor_hm <- cors %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("Exp1") %>% 
  pivot_longer(-Exp1, names_to = "Exp2") %>% 
  # mutate(Exp1 = factor(Exp1, levels = rownames(cors)),
  #        Exp2, factor(Exp2, levels = colnames(cors))
  #        ) %>% 
  ggplot(aes(x = Exp1, y = Exp2, fill = value)) +
  geom_tile() +
  geom_text(aes(label =round(value,2)),size = 2) +
  theme(axis.text.x = element_text(size = 5, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 5, angle = 45),
        axis.title = element_blank(),
        legend.key.size = unit(.25, "cm"),
        legend.margin = margin(0,0,0,0, "cm"),
        plot.margin = margin(0,0,0,0, "cm"),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5)) +
  scale_x_discrete(limits=rownames(cors)) + 
  scale_y_discrete(limits=colnames(cors)) +
  scale_fill_gradientn(colours = c("#2C7BB6",
                                          "#ABD9E9",
                                          "#FFFFBF",
                                          "#FDAE61",
                                          "#D7191C"),
                       limits = c(0,1),
                       breaks = seq(0,1, 0.2)
                       
  )


ggsave("results/figs/forPaper/gg_PRO_CorHM.pdf", cor_hm)
### PCA
pca <- prcomp(t(na.omit(proteinGroups_norm[ , experiment_cols])
                )
              )
pca_df <- as.data.frame(pca$x)
pca_df$run <- gsub("LFQ intensity ", "",
                   gsub("_P", "", rownames(pca_df)
                        )
                   ) %>%
   substr(0, nchar(.) - 3)
# pca_df$rep <- gsub("LFQ intensity ", "", rownames(pca_df)) %>%
#   substr(nchar(.) - 1, nchar(.))
pca_df$bait <- ifelse(grepl("R2A", pca_df$run), "FGFR2", 
                      ifelse(grepl("R11A", pca_df$run), "RAB11", "GFP")
                      )
eigs <- pca$sdev^2

pc1<-signif(100*(eigs[1] / sum(eigs)), 4)
pc2<-signif(100*(eigs[2] / sum(eigs)), 4)


pca_g <- ggplot(pca_df, aes(PC1, PC2, color= run, shape = bait)) +
  #ggforce::geom_mark_ellipse(aes(group = run, fill = bait), linetype = 0, alpha = 0.2) +
  geom_point(size = 3, alpha = 0.7)+
  scale_x_continuous(paste("PC1 (", pc1, "%)", sep=""))+
  scale_y_continuous(paste("PC2 (", pc2, "%)", sep=""))+
  scale_color_brewer("Run", palette = "Paired") +
  guides(color = "none", shape = "none") +
  #scale_fill_discrete(type = c("#D73027", "#4575B4", "#838B8B")) +
  theme(
    plot.title = element_text(size=20, face="bold",hjust = 0.5),
    panel.background = element_blank(),
    panel.grid = element_line("grey", linetype="dashed"), 
    legend.key = element_blank(),
    axis.line = element_line("black"),
    axis.text = element_text(size=6), 
    axis.title = element_text(size=6, face="bold")
  ) 
ggsave("results/figs/forPaper/proPCA.pdf", pca_g, 
       width = 60, height = 60, units = "mm", dpi = 300)

###
# Number identified in
# each experiment
###

proNo <- proteinGroups_flt_cc %>% 
  select(`id`, all_of(experiment_cols)) %>% 
  pivot_longer(cols = -`id`,
               names_to = "experiment",
               values_to = "intensity") %>% 
  filter(intensity != 0) %>% 
  mutate(experiment = gsub("LFQ intensity ", "", experiment)) %>% 
  mutate(rep = substr(experiment, nchar(experiment) - 2, nchar(experiment)),
         experiment = substr(experiment, 0, nchar(experiment) - 3)
  ) %>%
  group_by(experiment,rep) %>% 
  summarise(n = n(), .groups = "keep") %>%
  group_by(experiment) %>% 
  summarise(sd = sd(n, na.rm = T),
            n = mean(n), .groups = "keep") %>% 
  ungroup() %>% 
  unique() %>% 
  mutate(apex = ifelse(grepl("_A_", experiment),
                       "APEX", "Total"),
         experiment = gsub("_[A|T]", "", experiment)) %>% 
  ggplot(aes(x = experiment, y = n, fill = experiment,
             ymax = n+sd, ymin = n-sd
  )) +
  geom_col() +
  facet_wrap(~ apex, ncol = 1) +
  scale_fill_discrete(type = c("#838B8B",
                               "#D73027",
                               "#838B8B",
                               "#ABD9E9",
                               "#4575B4")
                      ) +
  geom_errorbar(width = 0.5) +
    theme(
      panel.background = element_blank(),
      panel.grid = element_line("grey", linetype="dashed"), 
      legend.key = element_blank(),
      axis.line = element_line("black"),
      # Text sizes may need modifying based on fig. sizes
      axis.text = element_text(size=12), 		
      axis.title = element_text(size=14, face="bold"),
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    ggtitle("Quantified Protein Groups")

ggsave("results/figs/forPaper/proNumQuant.pdf", proNo,
       width = 7, height = 7)



