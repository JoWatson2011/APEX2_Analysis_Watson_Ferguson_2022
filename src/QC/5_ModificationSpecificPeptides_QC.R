library(data.table)
library(tidyverse)

modsp <- readRDS("data/modificationSpecificPeptides.rds")
sty <- readRDS("data/phosphoSTY.rds")

# Filter out potential contaminants and reverse
modsp_flt <- modsp %>% filter(`Potential contaminant` != "+", 
                              `Reverse` != "+")

# Subset 'experiments' for phospho
phospho <- colnames(modsp_flt)[grepl("_Ph_", colnames(modsp_flt)) &
                                  grepl("Intensity ", colnames(modsp_flt))]

noSTY <- modsp_flt %>% 
  select(id, `Phospho (STY)`, all_of(phospho)) %>% 
  mutate(across(2:length(phospho), as.numeric)) %>% 
  pivot_longer(-c(id,`Phospho (STY)`)) %>%
  filter(value >0) %>% 
  mutate(APEX = ifelse(grepl("_A_", name), "Proximal", "Global")) %>% 
  group_by(`Phospho (STY)`,APEX) %>% 
  summarise(n=n(), .groups = "keep") %>% 
  filter(`Phospho (STY)` != 0) %>% 
  mutate(
    name = 
      ifelse(
        `Phospho (STY)` == 1, "Single",
        ifelse(
          `Phospho (STY)` == 2,"Double",
          ifelse(
            `Phospho (STY)` > 2,">Double", NA
          )
        )
      )
  ) %>% 
  group_by(name, APEX) %>% 
  summarise(n=sum(n), .groups = "keep") %>% 
  ungroup()


gg_global <- noSTY[noSTY$APEX=="Global",] %>% 
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

gg_proximal <- noSTY[noSTY$APEX=="Proximal",] %>% 
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

library(patchwork)
ggsave("results/figs/forPaper/Suppl_phosphoPie_noMods.pdf",
       gg_global+gg_proximal,
       width = 56,
       height = 56,
       unit = "mm",
       dpi = 300)

