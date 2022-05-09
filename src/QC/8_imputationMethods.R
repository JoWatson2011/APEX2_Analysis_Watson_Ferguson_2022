library(tidyverse)
library(patchwork)

sty <- readRDS("data/sty_FltMed.rds")
proteinGroups <- readRDS("data/proteinGroups_FltMed.rds")

experiment_sty <- gsub("_0[123]", "_med", readRDS("data/experiments_sty.RDS"))
experiment_pro <- gsub("_0[123]", "_med", readRDS("data/experiments_pro.RDS"))

imp <- lapply(seq(1.2, 2.5, 0.1), function(i){
  sty %>%
    select(all_of(experiment_sty)) %>%
    pivot_longer(cols = everything(),
                 names_to = "experiments",
                 values_to = "intensities") %>%  
    mutate(imput =  imputeLCMD::impute.QRILC(as.matrix(intensities), tune.sigma = i)[[1]],
           Imputed = is.na(intensities),
           tune = i)
}) %>% do.call(rbind, .)


imp_g <- ggplot(imp, aes(x = imput, fill = Imputed)) +
  geom_histogram(data = filter(imp, Imputed == T),
                 aes(x = imput),
                 #binwidth=0.2,
                 alpha = 0.9) +
  geom_histogram(data = filter(imp, Imputed == F),
                 aes(x = imput),
                 #binwidth=0.2,
                 alpha = 0.9) +
  theme(legend.position = "none") +
  ggtitle("Imputation from distribution of all median intensities",
          "Modified ") +
  facet_wrap(~ tune) +
  xlab("Normalised Intensities")
ggsave("results/figs/imputation/sty_impMedfromTotal.tiff", imp_g,
       width = 7, height = 7, units = "in")

imp_med_fct <- sty %>%
  pivot_longer(cols = all_of(experiment_sty),
               names_to = "experiments",
               values_to = "intensities") %>%
  mutate(Imputed = is.na(intensities)) %>%
  group_by(experiments) %>%
  mutate(imput =  imputeLCMD::impute.QRILC(as.matrix(intensities), tune.sigma = 2)[[1]])


imp_med_fct_g <- ggplot(imp_med_fct, aes(x = imput, fill = Imputed)) +
  geom_histogram(data = filter(imp_med_fct, Imputed == T),
                 aes(x = imput),
                 #binwidth=0.2,
                 alpha = 0.9) +
  geom_histogram(data = filter(imp_med_fct, Imputed == F),
                 aes(x = imput),
                 #binwidth=0.2,
                 alpha = 0.9) +
  facet_wrap(~ experiments) +
  theme(legend.position = "none") +
  ggtitle("Imputation from distribution of each condition median intensities") +
  xlab("Normalised Intensities")
ggsave("results/figs/imputation/sty_impMedfromCond.tiff", imp_med_fct_g,
       width = 7, height = 7, units = "in")



######
# proteinGroups
######
imp <- proteinGroups %>%
  select(all_of(experiment_pro)) %>%
  pivot_longer(cols = everything(),
               names_to = "experiments",
               values_to = "intensities") %>%  
  mutate(imput =  imputeLCMD::impute.QRILC(as.matrix(intensities), tune.sigma = 1.5)[[1]],
         Imputed = is.na(intensities))
imp_g <- ggplot(imp, aes(x = imput, fill = Imputed)) +
  geom_histogram(data = filter(imp, Imputed == T),
                 aes(x = imput),
                 #binwidth=0.2,
                 alpha = 0.9) +
  geom_histogram(data = filter(imp, Imputed == F),
                 aes(x = imput),
                 #binwidth=0.2,
                 alpha = 0.9) +
  theme(legend.position = "none") +
  ggtitle("Imputation from distribution of all median intensities") +
  xlab("Normalised Intensities")
ggsave("results/figs/imputation/proteinGroups_impMedfromTotal.tiff", imp_g,
       width = 7, height = 7, units = "in")

imp_med_fct <- proteinGroups %>%
  pivot_longer(cols = experiment_pro,
               names_to = "experiments",
               values_to = "intensities") %>%
  mutate(Imputed = is.na(intensities)) %>%
  group_by(experiments) %>%
  mutate(imput =  imputeLCMD::impute.QRILC(as.matrix(intensities), tune.sigma = 1.5)[[1]])


imp_med_fct_g <- ggplot(imp_med_fct, aes(x = imput, fill = Imputed)) +
  geom_histogram(data = filter(imp_med_fct, Imputed == T),
                 aes(x = imput),
                 #binwidth=0.2,
                 alpha = 0.9) +
  geom_histogram(data = filter(imp_med_fct, Imputed == F),
                 aes(x = imput),
                 #binwidth=0.2,
                 alpha = 0.9) +
  facet_wrap(~ experiments) +
  theme(legend.position = "none") +
  ggtitle("Imputation from distribution of each condition median intensities") +
  xlab("Normalised Intensities")
ggsave("results/figs/imputation/proteinGroups_impMedfromCond.tiff", imp_med_fct_g,
       width = 7, height = 7, units = "in")


