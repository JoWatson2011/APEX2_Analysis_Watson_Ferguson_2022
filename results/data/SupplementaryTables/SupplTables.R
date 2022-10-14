library(dplyr)
library(tidyr)

## STY
sty <-
  readRDS("data/sty_Flt.rds") %>% 
  select(-matches("med")) %>% 
  mutate(`Gene names` = gsub(";.*", "", `Gene names`),
         id_site = paste0(`Gene names`, "_", `Amino acid`, Position))

sty_sig <- readr::read_csv("results/data/APEX_STY_SIG_CLUSTERED.csv") %>% 
  mutate(Proximal_Regulated_Cluster = ifelse(
    !is.na(cl_name), "FGFR2_recycling","+"
    #ifelse(cl_name == "RE_profile2", "RE_profile", cl_name)
  )
  ) %>% 
  select(id_site, Proximal_Regulated_Cluster)
tot_sig <- readr::read_csv("results/data/GLOBAL_UPREG_40MIN.csv") %>% 
  select(Global_regulated = Significant, id_site = id_site)


merge_sty <- left_join(sty, sty_sig , by = "id_site") %>% 
  left_join(tot_sig, by = "id_site") 

# merge_sty %>% 
#   select(id_site, matches("Global"), matches("Proximal")) %>% 
#   pivot_longer(-id_site) %>% 
#   mutate(
#     group = ifelse(grepl("Global", name), "Global", "Proximal")
#   ) %>% 
#   filter(!is.na(value)) %>% 
#   group_by(id_site, group, value) %>% 
#   filter(n()>1)


sty_intcol <- grep("_[AT]", colnames(merge_sty), value = T)
# sty_normcol <- c("FGFR2-APEX 0",
#                  "FGFR2-APEX + FGF10 40",
#                  "RAB11-APEX + FGF10 40")
colnames(merge_sty)[colnames(merge_sty) %in% sty_intcol] <- paste0(
  ifelse(grepl("GFPA", sty_intcol), "GFP-APEX",
         ifelse(grepl("R11A", sty_intcol), "RAB11-APEX",
                ifelse(grepl("R2A", sty_intcol), "FGFR2-APEX", "")
         )
  ),
  "_",
  ifelse(grepl("_C_", sty_intcol), "UT", "FGF10"),
  "_",
  ifelse(grepl("_T_", sty_intcol), "Global", "Proximal"),
  "_R",
  substr(sty_intcol, nchar(sty_intcol), nchar(sty_intcol))
)
# colnames(sty_exp)[colnames(sty_exp) %in% sty_normcol] <- paste0(
#   sty_normcol,
#   "_Proximal_medianNormalised"
# )

#readr::write_csv(sty_exp, "results/data/SupplementaryTables/APEX_STY_data.csv")
readr::write_csv(merge_sty, "results/data/SupplementaryTables/SuppTable_6_SEP22.csv")


## Protein groups
pg <- readRDS("data/proteinGroups_Flt.rds") %>% 
  select(-contains("_med"))

pg_sig <- readr::read_csv("results/data/APEX_PRO_SIG.csv") %>% 
  mutate(Regulated = "+") %>% 
  select(id, Regulated)


merge_pg <- left_join(pg, pg_sig, by = "id")
pro_intcol <- grep("_[AT]_", colnames(merge_pg), value = T)
# pro_normcol <- grep("_A$", colnames(pro_exp), value = T)

colnames(merge_pg)[colnames(merge_pg) %in% pro_intcol] <- paste0(
  ifelse(grepl("GFPA", pro_intcol), "GFP-APEX",
         ifelse(grepl("R11A", pro_intcol), "RAB11-APEX",
                ifelse(grepl("R2A", pro_intcol), "FGFR2-APEX", "")
         )
  ),
  "_",
  ifelse(grepl("_C_", pro_intcol), "UT", "FGF10"),
  "_",
  ifelse(grepl("_T_", pro_intcol), "Global", "Proximal"),
  "_R",
  substr(pro_intcol, nchar(pro_intcol), nchar(pro_intcol))
)
# colnames(pro_exp)[colnames(pro_exp) %in% pro_normcol] <- paste0(
#   ifelse(grepl("R11A", pro_normcol), "RAB11-APEX",
#          ifelse(grepl("R2A", pro_normcol), "FGFR2-APEX", "")
#   ),
#   "_",
#   ifelse(grepl("_C_", pro_normcol), "UT", "FGF10"),
#   "_Proximal_medianNormalised"
# )

# readr::write_csv(pro_exp, "results/data/SupplementaryTables/APEX_Proteome_data.csv")
readr::write_csv(merge_pg, "results/data/SupplementaryTables/SuppTable_5_SEP22.csv")
