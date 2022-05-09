library(data.table)
library(tidyverse)
library(patchwork)



# Plot precursor error, score, intensity
evidence <- fread(input = "data/evidence.txt",
                  select = c("Mass error [ppm]", "Intensity", "Score",
                             "Sequence", "Gene names",
                             "Experiment",
                             "Phospho (STY)",
                             "Phospho (STY) site IDs"))


x <- evidence %>%
  filter(!is.na(Score)) %>% 
  mutate(ome = ifelse(grepl("_Ph_", Experiment), "Phosphopeptides", "Peptides"),
         apex = ifelse(grepl("_T_", Experiment), "Total", "APEX"))
med <- x %>% 
  group_by(ome, apex) %>% 
  summarise(med = median(Score), .groups = "keep")
andromScore <- ggplot(x, aes(Score)) + 
  geom_histogram(binwidth=5) +
  theme_classic() +
  scale_y_continuous(expand=expansion(mult = c(0, .2)), name = "Counts") +
  scale_x_continuous(name = "Andromeda Score",
                     limits= c(0, 350),
                     expand= expansion(mult = c(0, .2))) +
  facet_wrap(apex~ome) +
  geom_text(data = med, x = 300, y = 50000, size = 2, aes(label = paste("Median = ", med)))
ggsave("results/figs/forPaper/suppl_andromScore.pdf", 
       andromScore + theme(axis.text = element_text(size=6),
                           axis.title = element_text(size=6), 
                           strip.text = element_text(size=6)), 
       width = 85, height = 54, unit = "mm", dpi = 300)

intensities_ph <- evidence %>% 
  select(Experiment, Intensity) %>% 
  filter(grepl("_Ph_", Experiment)) %>% 
  mutate(apex = ifelse(grepl("_A_", Experiment), "APEX", "Total"),
         expmt = gsub("_[AT]_", "_", Experiment),
         order = ifelse(apex == "APEX", 2, 1))

intensities_ph_count <- ggplot(intensities_ph, aes(x = log10(Intensity), fill = apex)) +
  geom_histogram(data = filter(intensities_ph, apex == "Total"), 
                 aes(x = log10(Intensity)),
                 binwidth=0.2,
                 alpha = 0.9) +
  geom_histogram(data = filter(intensities_ph, apex == "APEX"),
                 aes(x = log10(Intensity)),
                 binwidth=0.2,
                 alpha = 0.9) +
  theme_classic() +
  theme(legend.title = element_blank()) +
  scale_y_continuous(expand=expansion(mult = c(0, .02)), name = "Counts") +
  scale_x_continuous(name = "Intensities (Log10)",
                     expand=expansion(mult = c(0, .02))) +
  ggtitle("Phospho Intensities")

intensities_ph_fct <- ggplot(intensities_ph, aes(fill = apex)) +
  #geom_histogram(binwidth=0.2, position = 'identity', alpha = 0.8) +
  geom_histogram(data = filter(intensities_ph, apex == "Total"), 
                 aes(x = log10(Intensity)),
                 binwidth=0.2,
                 alpha = 0.9) +
  geom_histogram(data = filter(intensities_ph, apex == "APEX"),
                 aes(x = log10(Intensity)),
                 binwidth=0.2,
                 alpha = 0.9) +
  theme_classic() +
  theme(strip.text.x = element_text(size = 6.5),
        legend.title = element_blank()) +
  scale_y_continuous(expand=expansion(mult = c(0, .02)), name = "Counts") +
  scale_x_continuous(name = "Intensities (Log10)",
                     expand=expansion(mult = c(0, .02))) +
  facet_wrap(~ expmt, ncol = 3)

ggsave("results/figs/QC/styIntensitiesHist.tiff",
       intensities_ph_count/ intensities_ph_fct +
         plot_layout(guides = 'collect'), 
       width = 210, height = 297, units = "mm")




intensities_prot <- evidence %>% 
  select(Experiment, Intensity) %>% 
  filter(grepl("_P_", Experiment)) %>% 
  mutate(apex = ifelse(grepl("_A_", Experiment), "APEX", "Total"),
         expmt = gsub("_[AT]_", "_", Experiment),
         order = ifelse(apex == "APEX", 2, 1))

intensities_prot_count <- ggplot(intensities_prot, aes(x = log10(Intensity), fill = apex)) +
  geom_histogram(data = filter(intensities_prot, apex == "Total"), 
                 aes(x = log10(Intensity)),
                 binwidth=0.2,
                 alpha = 0.9) +
  geom_histogram(data = filter(intensities_prot, apex == "APEX"),
                 aes(x = log10(Intensity)),
                 binwidth=0.2,
                 alpha = 0.9) +
  theme_classic() +
  theme(legend.title = element_blank()) +
  scale_y_continuous(expand=expansion(mult = c(0, .02)), name = "Counts") +
  scale_x_continuous(name = "Intensities (Log10)",
                     expand=expansion(mult = c(0, .02))) +
  ggtitle("Proteome Intensities")

intensities_prot_fct <- ggplot(intensities_prot, aes(fill = apex)) +
  #geom_histogram(binwidth=0.2, position = 'identity', alpha = 0.8) +
  geom_histogram(data = filter(intensities_prot, apex == "Total"), 
                 aes(x = log10(Intensity)),
                 binwidth=0.2,
                 alpha = 0.9) +
  geom_histogram(data = filter(intensities_prot, apex == "APEX"),
                 aes(x = log10(Intensity)),
                 binwidth=0.2,
                 alpha = 0.9) +
  theme_classic() +
  theme(strip.text.x = element_text(size = 6.5),
        legend.title = element_blank()) +
  scale_y_continuous(expand=expansion(mult = c(0, .02)), name = "Counts") +
  scale_x_continuous(name = "Intensities (Log10)",
                     expand=expansion(mult = c(0, .02))) +
  facet_wrap(~ expmt, ncol = 3)

ggsave("results/figs/QC/proIntensitiesHist.tiff",
  intensities_prot_count/ intensities_prot_fct +
    plot_layout(guides = 'collect'), 
       width = 210, height = 297, units = "mm")

library(patchwork)
ggsave("results/figs/forPaper/Suppl_histIntensities.pdf",
  (intensities_prot_count / intensities_ph_count) +plot_layout(guides ="collect") &
    scale_fill_manual(values = c("Total" = "#4575B4", "APEX" = "#D73027")) &
    theme(text = element_text(size = 5),
        panel.background = element_blank(),
        plot.margin = margin(0, 0.2, 0, 0.2, "mm")),
  width = 45,
  height = 71,
  unit = "mm",
  dpi = 300
)
