

# libraries ---------------------------------------------------------------



library(tidyverse)
library(patchwork)
library(janitor)
library(RColorBrewer)
library(RCCA)
library(biomaRt)
library(tidymodels)  
library(ggrepel)
library(raveio)
library(ggpubr)
library(ggchicklet)



# set theme for plots -----------------------------------------------------


theme_set(theme_bw() +
            theme(plot.title = element_text(size = 18),
                  axis.title = element_text(size = 15),
                  axis.text = element_text(size = 12),
                  strip.text = element_text(size = 15),
                  legend.title = element_text(size = 15),
                  legend.text = element_text(size = 12)))




# data --------------------------------------------------------------------


  soft_power <- 3
  base_dir <- "/Users/work/Library/CloudStorage/OneDrive-NationalInstitutesofHealth/projects/WGCNA/acsg_wgcna/"
  prefix <- "20230228_185samples_19kgenes_vst_qSVA123567_MP_RNAbatch_Race_resids"

# LOAD OBJECTS  
  load(paste0(base_dir, "objects/", prefix, ".RDS")) # df_vsd_regress
  load(paste0(base_dir, "objects/", prefix, "_SIGNED_SFT", soft_power, "_MODS.RDS")) # df_modules

# IDENTIFY MODULE SET OF INTEREST  
  df_modules_filt <- df_modules %>% 
    filter(min_size == 40 & cut_height == 0.97) %>% 
    unite("mod_set", c(sft, min_size, cut_height), sep = "_") %>% 
    arrange(mod_set, module)

# COVARIATES
  load(paste0(base_dir, "data/covariates/185_all_covariates_clean.Rdata"))

# ADD DRUG CLASSES TO COVARIATES
  drugs_of_abuse <- c("smoker", "nicotine", "alcohol", "cannabis", "cocaine", "opioids", "major_stimulants")
  mood_meds <- c("anti_anxiety", "mood_stabilizers", "anti_epileptics")
  other_psych_meds <- c("anticholinergics", "anti_depressant", "antipsychotics" )
  other_drugs <- c("anti_histamines", "non_psychiatric", "other_psychotropic_drug")
  
  df_covariates_all_clean_numeric <- df_covariates_all_clean_numeric %>% 
    dplyr::select(-bmi, -manner) %>% 
    
    mutate(
      
      drugs_of_abuse = ifelse(df_covariates_all_clean_numeric %>% 
                                dplyr::select(all_of(drugs_of_abuse)) %>% 
                                rowSums(na.rm = TRUE) != 0, 1, 0),
      mood_meds = ifelse(df_covariates_all_clean_numeric %>% 
                           dplyr::select(all_of(mood_meds)) %>% 
                           rowSums(na.rm = TRUE) != 0, 1, 0),
      other_psych_meds = ifelse(df_covariates_all_clean_numeric %>% 
                                  dplyr::select(all_of(other_psych_meds)) %>% 
                                  rowSums(na.rm = TRUE) != 0, 1, 0),
      other_drugs = ifelse(df_covariates_all_clean_numeric %>% 
                             dplyr::select(all_of(other_drugs)) %>% 
                             rowSums(na.rm = TRUE) != 0, 1, 0)
      
    ) %>% dplyr::select(-all_of(c(drugs_of_abuse, mood_meds, other_psych_meds, other_drugs)))
  
# GENE SYMBOL TO ENSEMBL GENE ID MAPPING
  load(paste0(base_dir, "objects/ensembl_id_to_gene_symbol.Rdata")) # df_ensembl_to_symbol
  load(paste0(base_dir, "objects/df_hsapiens_genome.RDS")) # df_hsapiens_genome


# define matrices ---------------------------------------------------------------


# Select covariates based on significant covariates from 4.1: 
# anti-anxiety, opioids, brain_weight, age_death

# DEFINE MATRICES  
  df_grcca_x <- df_vsd_regress %>% 
    dplyr::select(ensembl_gene_id, sample, resids) %>% 
    left_join(df_modules_filt) %>% 
    arrange(module) 
  
  X.mat <- df_grcca_x %>% 
    pivot_wider(id_cols = sample, names_from = ensembl_gene_id, values_from = resids) %>% 
    as.data.frame %>% 
    column_to_rownames("sample")
  
  x_group <- df_grcca_x %>% 
    dplyr::select(ensembl_gene_id, module) %>% 
    distinct() %>% 
    pull(module) %>% 
    as.numeric

# ADD COVARIATES AND CONFOUNDERS
  
  technical_covariates <- c("mapped_percent", "five_prime_three_prime_bias", "rin_acsg", 
                            "rn_aextraction_batch", "library_batch", "pmi", "g_cpercent", 
                            "ph")
  
  confounds <- c("age_death", "brain_weight", "height", "weight", "race")
  
  # deal with NAs? remove covariates at a certain threshold of NAs?
  df_covariates_all_clean_numeric %>% 
    dplyr::select(-all_of(technical_covariates)) %>% 
    dplyr::select(-all_of(confounds)) %>% 
    mutate(case = ifelse(multi_nomial_dx == 1, 0, 1)) %>% 
    dplyr::select(sample, case, everything(), -multi_nomial_dx) %>% 
    is.na %>% colSums()
  
  # remove BMI because we already have height & weight
  # marital status unknown = 0
  # manner of death or suicide unknown = 0
  # if drugs are unknown, set to 0 (only using information if they definitively had drug in system)
  
  df_covariates_all_clean_numeric_noNA <- df_covariates_all_clean_numeric %>% 
    dplyr::select(-all_of(technical_covariates)) %>% 
    mutate_if(is.numeric, ~ ifelse(is.na(.x), 0, .x))
  
  # from univariate:
  # [1] "anti_epileptics"         "nicotine"                "alcohol"                
  # [4] "non_psychiatric"         "race"                    "opioids"                
  # [7] "anticholinergics"        "manner"                  "cocaine"                
  # [10] "other_psychotropic_drug" "antipsychotics"  
  
  covariates <- c("gender", "source", "marital_status", "suicide",
                  "drugs_of_abuse", "mood_meds", "other_psych_meds", "other_drugs")
  
  YC.mat <- df_covariates_all_clean_numeric_noNA %>% 
    mutate(case = ifelse(multi_nomial_dx == 1, 0, 1)) %>% 
    dplyr::select(-multi_nomial_dx) %>% 
    dplyr::select(sample, case, all_of(covariates), all_of(confounds)) %>% 
    column_to_rownames("sample")
  
  Y.mat <- YC.mat %>% dplyr::select(-all_of(confounds))
  colnames(Y.mat)
  C.mat <- YC.mat %>% dplyr::select(all_of(confounds))
  colnames(C.mat)
  
# EXPORT FOR AGOSTON'S TOOLBOX
  
  cca_dir <- paste0(base_dir, "RCCA_toolkit/GENES_caseControl_4drugClass_5confound/")
  dir.create(cca_dir)
  cca_data_dir <- paste0(cca_dir, "data/")
  dir.create(cca_data_dir)
  
  write.table(X.mat, file = paste0(cca_data_dir, "X.txt"))
  write.table(Y.mat, file = paste0(cca_data_dir, "Y.txt"))
  write.table(C.mat, file = paste0(cca_data_dir, "C.txt"))
  write.table(x_group, file = paste0(cca_data_dir, "XGroup.txt"))

# GRCCA LABELS FILES
  
  df_labels_x <- tibble(ensembl_gene_id = colnames(X.mat)) %>% 
    left_join(df_modules_filt) %>% 
    left_join(df_hsapiens_genome %>% 
                dplyr::select(gene_id, gene_name) %>% 
                distinct %>% 
                dplyr::rename("ensembl_gene_id" = "gene_id")) %>% 
    mutate(Label = row_number(), .before = 1) %>% 
    dplyr::select(Label, module, ensembl_gene_id, gene_name) %>% 
    dplyr::rename_all(~c("Label", "Category", "Ensembl_id", "HGNC_id")) %>% 
    mutate(HGNC_id = ifelse(is.na(HGNC_id), Ensembl_id, HGNC_id))
  
  df_labels_y <- tibble(Label = colnames(Y.mat))
  
  write.csv(df_labels_x, paste0(cca_dir, "data/LabelsX.csv"), row.names = FALSE)
  write.csv(df_labels_y, paste0(cca_dir, "data/LabelsY.csv"), row.names = FALSE)
  
  
  dim(X.mat)
  dim(Y.mat)
  dim(C.mat)
  length(x_group)  
  
  
  


