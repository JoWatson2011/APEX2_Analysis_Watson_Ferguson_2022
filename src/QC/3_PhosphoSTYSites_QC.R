library(dplyr)
library(data.table)
library(readr)
library(ggplot2)
library(tidyr)
library(patchwork)

###
#Import Data
###
sty <- readRDS("data/phosphoSTY.RDS")
experiment_cols <- grep("Intensity",
                        colnames(sty),
                        value = T)

###
# Filter
###
sty[,experiment_cols] <- apply(sty[,experiment_cols],
                                         2,
                                         as.numeric)

# Summary of variables Phospho (STY)Sites.txt will be filtered on.
# After running the script this information will be stored in the
# variable named report.
report <- vector()
report["Total Sites Identified"] <- nrow(sty)
report["Potential Contaminants"] <- sum(sty$`Potential contaminant` == "+")
report["Reverse"] <- sum(sty$Reverse == "+")
report["Localisation probability < 0.75"] <- sum(sty$`Localization prob` < 0.75)
report["No detected phosphorylation (incomplete cases)"] <-
  sty %>%
  select(all_of(experiment_cols)) %>%
  filter(rowSums(. == 0 ) == ncol(.)) %>% nrow()

# Filter Phospho (STY)Sites.txt to remove contaminants, reverse, and low localisation
# probability sites. Localisation probability threshold can be changed here.
sty_flt <- sty %>%
  filter(`Potential contaminant` != "+",
         `Reverse` != "+", 
         `Localization prob` > 0.75)

# Add info to report
report["Sites remaining following filtering"] <- nrow(sty_flt)
report["Incomplete cases in filtered dataset"] <- sty_flt %>% 
  select(all_of(experiment_cols)) %>%
  filter(rowSums(. == 0 ) == ncol(.)) %>%
  nrow()

# Filter rows with no quantitative info
sty_flt_cc <- sty_flt[rowSums(sty_flt[experiment_cols] > 0) >= 1,] 

###
# Normalise
#
# Proximal and global are normalised separately.
# See plot below
###
proximal_cols <- grep("_A_", experiment_cols, value = T)
global_cols <- grep("_T_", experiment_cols, value = T)

sty_log <- apply(sty_flt_cc[,experiment_cols], 2, function(i){
  tmp <- ifelse(i == 0, NA, i)
  tmp2 <- log10(tmp)
  return(tmp2)
})

sty_norm <- cbind(
  sty_flt_cc[,!(names(sty_flt_cc) %in% c(proximal_cols, global_cols))],
  limma::normalizeQuantiles(sty_log[,proximal_cols]),
  limma::normalizeQuantiles(sty_log[,global_cols])
) 

styIntUnnorm <- sty_flt_cc %>%
  select(id, all_of(experiment_cols)) %>% 
  mutate_all(~ ifelse(. == 0, NA, .)) %>% 
  mutate_at(vars(experiment_cols), log2) %>% 
  pivot_longer(cols = -c(id),
               names_to = "experiment",
               values_to = "intensity") %>% 
  mutate(rep = substr(experiment,
                      nchar(experiment)-1,
                      nchar(experiment))) %>% 
  ggplot(aes(x = intensity, color = rep, group = experiment)) +
  geom_density() +
  theme(legend.position = 'none')

styIntNorm <- sty_norm %>% 
  select(id, all_of(experiment_cols)) %>% 
  pivot_longer(cols = -c(id),
               names_to = "experiment",
               values_to = "intensity") %>% 
  mutate(rep = substr(experiment,
                      nchar(experiment)-1,
                      nchar(experiment))) %>% 
  ggplot(aes(x = intensity, color = rep, group = experiment)) +
  geom_density() +
  theme(legend.position = 'none') +
  ggtitle("", subtitle = "Normalised Intensities")

styInt_Unnorm_f <- sty_flt_cc %>% 
  select(experiment_cols) %>% 
  mutate_all(~ ifelse(. == 0, NA, .)) %>% 
  mutate_all(~ log2(.)) %>% 
  pivot_longer(cols = everything(),
               names_to = "experiment",
               values_to = "intensity") %>% 
  mutate(rep = substr(experiment, nchar(experiment)-1, nchar(experiment)),
         apex = ifelse(grepl("_T_", experiment), "Total", "Apex")) %>% 
  ggplot(aes(x = intensity, color = rep, group = experiment)) +
  geom_density() +
  theme(legend.position = 'none') +
  ggtitle("Phospho (STY) Intensities",
          subtitle = "Unnormalised Intensities") +
  facet_grid(~ apex)


ggsave("results/figs/QC/styNormCurve.tiff",
       styInt_Unnorm_f / styIntUnnorm / styIntNorm, 
       width = 210, height = 297, units = "mm")

ggsave("results/figs/forPaper/Suppl_PhosphoIntensities.pdf",
       styInt_Unnorm_f /  styIntNorm  &
         theme(text = element_text(size = 5), 
               plot.title = element_blank(),
               panel.background = element_blank()),
       width = 46,
       height = 70,
       unit = "mm",
       dpi = 300)


saveRDS(sty_norm, "data/sty_Flt.rds")
readr::write_csv(sty_norm, "data/sty_Flt.csv")

# Pie chart to visualise categorical characteristics.
ph_loc_pie <- sty %>% 
  select(id, `Localization prob`, all_of(experiment_cols)) %>% 
  # mutate(across(2:length(phospho), as.numeric)) %>% 
  pivot_longer(-c(id,`Localization prob`)) %>%
  filter(value >0) %>% 
  mutate(APEX = ifelse(grepl("_A_", name), "Proximal", "Global")) %>% 
  mutate(
    name = 
      ifelse(
        `Localization prob` > 0.75, "Class I", "Class II"
      )
  ) %>% 
  select(-value) %>% 
  unique() %>% 
  group_by(name,APEX) %>% 
  summarise(n=n(), .groups = "keep") %>% 
  ungroup()

ph_loc_pie_prox <- ph_loc_pie[ph_loc_pie$APEX=="Proximal",] %>% 
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
  scale_fill_manual(values = c("#D73027", "#4575B4", "#ABD9E9")) +
  scale_color_manual(values = c("#D73027", "#4575B4", "#ABD9E9")) +
  coord_polar(theta="y") +
  xlim(c(-2, 3)) +
  theme_void() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)#,
        #plot.margin = margin(0.1,0.1,0.1,0.1, "mm")
  )
ph_loc_pie_global <- ph_loc_pie[ph_loc_pie$APEX=="Global",] %>% 
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
  scale_fill_manual(values = c("#D73027", "#4575B4", "#ABD9E9")) +
  scale_color_manual(values = c("#D73027", "#4575B4", "#ABD9E9")) +
  coord_polar(theta="y") +
  xlim(c(-2, 3)) +
  theme_void() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)#,
        #plot.margin = margin(0.1,0.1,0.1,0.1, "mm")
  )


ph_aa_pie <- sty %>% 
  select(id, `Amino acid`, all_of(experiment_cols)) %>% 
  # mutate(across(2:length(phospho), as.numeric)) %>% 
  pivot_longer(-c(id,`Amino acid`)) %>%
  filter(value >0) %>% 
  mutate(APEX = ifelse(grepl("_A_", name), "Proximal", "Global")) %>% 
  select(-value) %>% 
  unique() %>% 
  group_by(`Amino acid`,APEX) %>% 
  summarise(n=n(), .groups = "keep") %>% 
  ungroup() %>% 
  rename(name = `Amino acid`)

ph_aa_pie_prox <- ph_aa_pie[ph_aa_pie$APEX=="Proximal",] %>% 
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
  scale_fill_manual(values = c("#D73027", "#4575B4", "#ABD9E9")) +
  scale_color_manual(values = c("#D73027", "#4575B4", "#ABD9E9")) +
  coord_polar(theta="y") +
  xlim(c(-2, 3)) +
  theme_void() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)#,
        #plot.margin = margin(0.1,0.1,0.1,0.1, "mm")
  )
ph_aa_pie_global <- ph_aa_pie[ph_aa_pie$APEX=="Global",] %>% 
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
  scale_fill_manual(values = c("#D73027", "#4575B4", "#ABD9E9")) +
  scale_color_manual(values = c("#D73027", "#4575B4", "#ABD9E9")) +
  coord_polar(theta="y") +
  xlim(c(-2, 3)) +
  theme_void() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)#,
        #plot.margin = margin(0.1,0.1,0.1,0.1, "mm")
  )

ggsave("results/figs/forPaper/Suppl_phosphoPies.pdf",
       ph_loc_pie_global  + ph_loc_pie_prox +
         ph_aa_pie_global  + ph_aa_pie_prox,
  width = 90,
  height = 56,
  unit = "mm",
  dpi = 300)

###
# Correlation between conditions
###
cors <- sty %>% 
  select(c(grep("_A_", colnames(.)),
           grep("_T_", colnames(.)))
  ) %>% 
  cor(use = "pairwise.complete.obs")
rownames(cors) <- colnames(cors) <- gsub("Intensity ", "", colnames(cors))

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

ggsave("results/figs/forPaper/gg_STY_CorHM.pdf", cor_hm)

###
# Calculate PCA
###
pca <- prcomp(t(na.omit(sty_norm[ , experiment_cols])
                )
              )
pca_df <- as.data.frame(pca$x)
pca_df$run <- gsub("Intensity ", "",
                   gsub("_Ph", "", rownames(pca_df)
)) %>%
  substr(0, nchar(.) - 3)
# pca_df$rep <- gsub("Intensity ", "", rownames(pca_df)) %>%
#   substr(nchar(.) - 1, nchar(.))
pca_df$bait <- ifelse(grepl("R2A", pca_df$run), "FGFR2", 
                      ifelse(grepl("R11A", pca_df$run), "RAB11", "GFP")
)
pca_df$biotin <- ifelse(grepl("_A_", pca_df$run), 
                        "APEX",
                        "Total")
eigs <- pca$sdev^2

pc1<-signif(100*(eigs[1] / sum(eigs)), 4)
pc2<-signif(100*(eigs[2] / sum(eigs)), 4)


pca_g <- ggplot(pca_df, aes(PC1, PC2, color= run, shape = bait)) +
  #ggforce::geom_mark_ellipse(aes(group = run, fill = bait), linetype = 0, alpha = 0.2) +
  geom_point(size = 3, alpha = 0.7)+
  scale_x_continuous(paste("PC1 (", pc1, "%)", sep=""))+
  scale_y_continuous(paste("PC2 (", pc2, "%)", sep=""))+
  scale_color_brewer("Run", palette = "Paired") +
 # guides(color = "none") +
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
ggsave("results/figs/forPaper/styPCA.pdf", pca_g, 
       width = 60, height = 60, units = "mm", dpi = 300)

