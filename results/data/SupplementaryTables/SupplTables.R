library(dplyr)
library(tidyr)
sty <-
  readRDS("data/sty_FltMed.rds") %>% 
  select(-matches("med")) %>% 
  mutate(`Gene names` = gsub(";.*", "", `Gene names`),
         id_site = paste0(`Gene names`, "_", `Amino acid`, Position))

pg <- readRDS("data/proteinGroups_FltMed.rds") %>% 
  select(-contains("_med"))


sty_sig <- readr::read_csv("results/data/APEX_STY_SIG_CLUSTERED.csv") %>% 
  mutate(cl_name = ifelse(
    cl_name == "RE_profile", "FGFR2_Recycling","+"
    #ifelse(cl_name == "RE_profile2", "RE_profile", cl_name)
  )
   )# %>% 
  # select(id, id_site, 2:8, 11, 13:14, Cluster = cl_name, 
  #        grep("_A_", colnames(.)))
#readr::write_csv(sty_sig, "results/data/SupplementaryTables/APEX_PhosphoSTY.txt")

pg_sig <- readr::read_csv("results/data/APEX_PRO_SIG_CLUSTERED.csv") %>% 
  mutate(cl_name = ifelse(
    cl_name == "RE_profile", "FGFR2_Recycling", "+"
  )) #%>% 
  # select(id, `Gene names`, 3:5,9:10,
  #        grep("_A", colnames(.)),
  #        Cluster = cl_name)


tot_sig <- readr::read_csv("results/data/GLOBAL_UPREG_40MIN.csv") 

# tot <- readr::read_csv("results/data/STY_Total_sum.csv") %>% 
#   left_join(sty[,c(54,2:3,5:7,10:13)],
#             by = c("id" = "id_site")) %>% 
#   select(1:2,18:26,grep("Control", colnames(.)), grep("FGF", colnames(.))) %>% 
#   mutate(Regulated = ifelse(id %in% tot_sig$id, "+", ""))
#readr::write_csv(tot, "results/data/SupplementaryTables/APEX_Global.txt")


sty_exp <- sty %>% 
  left_join(sty_sig[,c(44:48)], by = "id_site") %>% 
  mutate(Proximal_Regulated_Cluster = ifelse(is.na(cl_name), "", cl_name)) %>% 
  left_join(tot_sig[,c("id", "Significant")], by = c("id_site" = "id")) %>% 
  mutate(Global_Regulated = ifelse(is.na(Significant), "", Significant)) %>%
  select(-c("Reverse", "Potential contaminant","Localization prob",
            "Mod. peptide IDs","Fasta headers",
            "Number of Phospho (STY)", "Significant",
            "cl_name")) %>% 
  rename(Uniprot = Protein)

sty_intcol <- grep("_[AT]", colnames(sty_exp), value = T)
sty_normcol <- c("FGFR2-APEX 0",
                 "FGFR2-APEX + FGF10 40",
                 "RAB11-APEX + FGF10 40")
colnames(sty_exp)[colnames(sty_exp) %in% sty_intcol] <- paste0(
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
colnames(sty_exp)[colnames(sty_exp) %in% sty_normcol] <- paste0(
  sty_normcol,
  "_Proximal_medianNormalised"
)

pro_exp <- pg_sig %>% 
  mutate(Regulated = ifelse(is.na(cl_name), "", cl_name)) %>% 
  select(-c("Sequence length", "Razor + unique peptides",
            "Unique + razor sequence coverage [%]",
            "cl_name")) 
pro_intcol <- grep("_[AT]_", colnames(pro_exp), value = T)
pro_normcol <- grep("_A$", colnames(pro_exp), value = T)

colnames(pro_exp)[colnames(pro_exp) %in% pro_intcol] <- paste0(
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
colnames(pro_exp)[colnames(pro_exp) %in% pro_normcol] <- paste0(
  ifelse(grepl("R11A", pro_normcol), "RAB11-APEX",
         ifelse(grepl("R2A", pro_normcol), "FGFR2-APEX", "")
  ),
  "_",
  ifelse(grepl("_C_", pro_normcol), "UT", "FGF10"),
  "_Proximal_medianNormalised"
)

readr::write_csv(pro_exp, "results/data/SupplementaryTables/APEX_Proteome_data.csv")
readr::write_csv(sty_exp, "results/data/SupplementaryTables/APEX_STY_data.csv")
