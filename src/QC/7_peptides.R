library(data.table)
library(tidyverse)
peptides <- readRDS("data/peptides.rds")

peptides_flt <- peptides %>%
  filter(`Potential contaminant` != "+", `Reverse` != "+") %>% 
  select(id, grep("Experiment", colnames(peptides)))

tmp <- peptides_flt %>% 
  pivot_longer(cols = -id) %>% 
  mutate(rep = substr(name, 31, 32),
         grp = substr(name, 12, 30),
         name = substr(name, 12, 32),
         ome = ifelse(grepl("_Ph_", name), "Phospho", "Prot")) %>% 
  group_by(id, name) %>%
  filter(any(is.na(value))) %>%
  group_by(grp, name, rep, ome) %>% 
  summarise(n = n()) %>% 
  ungroup()


ggplot(tmp, aes(x = name, y = n, fill = grp)) + 
  geom_bar(stat = "identity") +
  scale_y_continuous(name = "Peptides identified") +
  scale_x_discrete(name = NULL) +
  scale_fill_discrete(guide = F) +
  labs(title = "Identified peptides: peptides.txt") +
  facet_wrap(~ome, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.title = element_blank())
  
id_count <- data.frame(experiment = colnames(peptides_flt[,-1]),
                       id = apply(peptides_flt[,-1], 2, function(x) sum(!is.na(x))),
                       colour = substr(colnames(peptides_flt[,-1]), 12,19))

g <- ggplot(id_count, aes(x=experiment, y=id, fill = colour)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(name = "Peptides identified") +
  scale_x_discrete(name = NULL) +
  labs(title = "Identified peptides: peptides.txt") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.title = element_blank())

ggsave(plot = g, filename = "results/figs/identified_peps.tif", device = "tiff")
