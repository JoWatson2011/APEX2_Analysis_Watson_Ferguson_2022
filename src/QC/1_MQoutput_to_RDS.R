# This script takes the MQ outputs and compresses
# the data so it is quickly read into R.
# These are the data available in this repo.

library(dplyr)
library(data.table)

# Read experiment names
experiments <- (fread(input = "data/summary.txt") %>% 
                #, select = "Experiment") %>%
                  unique %>%
                  #filter(grepl("STY", Experiment)) %>% 
                  mutate(Ligand = substr(Experiment, 5, 9),
                         Timepoint = substr(Experiment, 16, 17),
                         rep = substr(Experiment, 19, 20)) %>% 
                  group_by(Ligand) %>% 
                  arrange(Timepoint, .by_group=T))$Experiment
experiments_cols <- unique(experiments[experiments != ""])

# Use experiment names to get column names for phospho (STY)Sites.txt
experiments_sty <- sapply(experiments_cols[grepl("_Ph_", experiments_cols)], function(i){
  paste0("Intensity ", i)
}, USE.NAMES = F)

sty <- fread(input = "data/Phospho (STY)Sites.txt",
             select = c("id", "Protein", "Protein names",
                        "Gene names",	"Amino acid", "Position",
                        "Sequence window", "Reverse", "Potential contaminant",
                        "Localization prob", "Mod. peptide IDs", "Fasta headers",
                        "Number of Phospho (STY)",
                        experiments_sty
             ), integer64 = "numeric"
             ) %>% as.data.frame()

saveRDS(sty, "data/PhosphoSTY.RDS")


# Use experiment names to get column names for proteinGroups.txt
experiments_pro <- sapply(experiments_cols[grepl("_P_", experiments_cols)], function(i){
  paste0("LFQ intensity ", i)
}, USE.NAMES = F)

pg <- fread(input = "data/proteinGroups.txt",
            select = c("id", "Majority protein IDs", "Sequence length",
                       "Gene names", "Protein names", "Reverse",
                       "Potential contaminant", "Fasta headers",
                       "Razor + unique peptides",
                       "Unique + razor sequence coverage [%]",
                       "Only identified by site",
                       experiments_pro), integer64 = "numeric"
) %>% as.data.frame()
saveRDS(pg, "data/proteinGroups.RDS")

# Same for modificationSpecificPeptides
modsp <- fread(input = "data/modificationSpecificPeptides.txt",
               select = c("Sequence", "Mass", "Mass Fractional Part",	
                          "Protein Groups", "Proteins",	"Gene Names",
                          "Protein Names", "Unique (Groups)",
                          "Unique (Proteins)", "Phospho (STY)",
                          paste("Experiment", experiments_cols),
                          paste("Intensity", experiments_cols),
                          "Reverse", "Potential contaminant",
                          "id", "Phospho (STY) site IDs"),
               integer64 = "numeric") %>% 
  as.data.frame()

saveRDS(modsp, "data/modificationSpecificPeptides.rds")

#peptides.txt

peptides <- fread(input = "data/peptides.txt",
                  select = c("id", "Sequence", "Length", "Mass", "Proteins", 
                             "Potential contaminant", "Reverse",
                             paste("Experiment", experiments_cols)),
                  integer64 = "numeric"
) %>% as.data.frame()


saveRDS(peptides, "data/peptides.rds")
