

# libraries ---------------------------------------------------------------

library(tidyverse)
library(tidymodels)
library(ggrepel)

# load data ---------------------------------------------------------------

base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/objects/"
load(paste0(base_dir, "20230228_185samples_19kgenes_vst_qSVA123567_MP_RNAbatch_Race_resids.RDS"))

# get drugs only, convert to df -------------------------------------------

n_samples <- df_vsd_regress %>% pull(sample) %>% unique %>% length

psych_drugs <- c("anti_anxiety", "anti_depressant", "antipsychotics", "mood_stabilizers", 
                 "anti_epileptics", "anti_histamines", "anticholinergics", "other_psychotropic_drug")
rec_drugs <- c("smoker", "nicotine", "cocaine", "alcohol", "opioids", "cannabis", 
               "major_stimulants", "non_psychiatric")

df_drugs <- df_vsd_regress %>% 
  ungroup %>% 
  select(sample, all_of(c(psych_drugs, rec_drugs))) %>% 
  distinct() %>% 
  
  # combine drug classes that are essentially the same
  mutate(anti_histamines = ifelse(anti_histamines == 1 | anticholinergics == 1, 1, 0),
         nicotine = ifelse(nicotine == 1 | smoker == 1, 1, 0)
  ) %>% 
  dplyr::select(-anticholinergics, -smoker)

# count NAs for each drug
prop_na <- (df_drugs %>% 
  dplyr::select(-sample) %>% 
  is.na %>% 
  colSums) / n_samples

# Plot NAs by sample
df_drugs %>% 
  pivot_longer(2:ncol(.), names_to = "drug", values_to = "value") %>% 
  mutate(dx = str_remove(sample, ".*_") %>% factor(levels = c("control", "bipolar", "mdd", "schizo")),
         value = factor(value)) %>% 
  arrange(dx, sample) %>% 
  mutate(sample = factor(sample, levels = unique(.$sample))) %>% 
  
  ggplot(aes(x = drug, y = sample, fill = value)) +
  geom_tile() +
  geom_tile(aes(x = 0, y = sample, fill = dx)) +
  scale_fill_manual(values = c("black", "white", "red", "yellow", "green", "blue")) +
  coord_flip() +
  theme(axis.text.x = element_blank())

# remake df_drugs for analysis
df_drugs2 <- df_drugs %>% 
  
  # remove drugs with > 5% missing
  dplyr::select(sample, all_of(names(prop_na[prop_na < 0.1]))) %>% 
  
  # for other NAs, turn to 0 (assume patient did not use)
  mutate_if(is.numeric, ~ifelse(is.na(.x), 0, .x))


# PCA ---------------------------------------------------------------------

# RUN PCA
drugs_pca <- df_drugs2 %>% 

  # run PCA
  as.data.frame %>% 
  column_to_rownames("sample") %>% 
  prcomp(center = TRUE, scale = TRUE)

# GET PC SCORES FOR EACH DRUG
df_drugs_pca <- drugs_pca$rotation %>% 
  as.data.frame() %>% 
  rownames_to_column("drug") %>% 
  as_tibble()

# PLOT
df_drugs_pca %>% 
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point() +
  geom_hline(aes(yintercept = 0), lty = 2) +
  geom_vline(aes(xintercept = 0), lty = 2) +
  geom_text_repel(aes(label = drug))


# factor analysis ---------------------------------------------------------

drug_data <- df_drugs2 %>% 
  as.data.frame %>% 
  column_to_rownames("sample")

library(psych)

fa <- fa(r = drug_data,
         nfactors = 4,
         rotate = "varimax")
summary(fa)
fa$loadings

# KMO test
# * 0.00 to 0.49 unacceptable
# * 0.50 to 0.59 miserable
# * 0.60 to 0.69 mediocre
# * 0.70 to 0.79 middling
# * 0.80 to 0.89 meritorious
# * 0.90 to 1.00 marvelous

KMO(drug_data) # remove cocaine (only one patient used), alcohol as its own

KMO(drug_data %>% dplyr::select(-cocaine)) # nicotine, cannabis, alcohol not great

KMO(drug_data %>% dplyr::select(-c(cocaine, nicotine, cannabis, alcohol)))
drug_data_filt <- drug_data %>% dplyr::select(-c(cocaine, nicotine, cannabis, alcohol))

# GET EIGENVALUES
ev <- eigen(cor(drug_data_filt)) # get eigenvalues
ev$values
scree(drug_data_filt, pc=FALSE)

factor_analysis <- factanal(drug_data,
                            factors = 2,
                            rotation = "varimax")

fa.diagram(factor_analysis$loadings)

df_loadings <- factor_analysis$loadings[,1:4] %>% 
  as.data.frame() %>% 
  rownames_to_column("drug") %>% 
  as_tibble()




