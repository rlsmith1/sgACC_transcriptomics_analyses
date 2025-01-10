
################################################################################

# Load libraries, data, and functions for gene-level analyses

################################################################################

base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
prefix <- "08Mar2024_GENES_qSVAgeSexRaceGC"

### LIBRARIES ###

library(tidyverse)
library(tidymodels)
library(tidytext)

library(janitor)

library(patchwork)
library(ggrepel)
library(ggh4x)
library(ggpubr)
library(cowplot)
library(RColorBrewer)
library(pals)

library(biomaRt)
library(clusterProfiler)
library(msigdbr)
library(org.Hs.eg.db)
library(EnsDb.Hsapiens.v79) # hg38
library(ensembldb)
library(DESeq2)

library(rhdf5)

library(readxl)
library(writexl)

### DATA ###

load(paste0(base_dir, "objects/", prefix, ".RDS")) # df_vsd_regress
load(paste0(base_dir, "objects/14Nov2023_covariates.Rdata")) # df_covariates, df_covariates_numeric (generated in clean_covariates.R)
load(paste0(base_dir, "objects/akula_et_al_DE_results.RDS")) # df_akula_res
load(paste0(base_dir, "objects/drug_MCA_results.Rdata")) # df_ind_loadings, df_var_loadings

### Ensembl ID to gene symbol mappings ###
ensembl_ids <- df_vsd_regress %>% dplyr::select(-sample) %>% colnames
df_ensembl_to_symbol <- mapIds(
    EnsDb.Hsapiens.v79, 
    keys = ensembl_ids,
    column = "SYMBOL", 
    keytype = "GENEID", 
    multiVals = "first"
) %>% 
    enframe(name = "ensembl_gene_id", value = "gene_symbol") %>% 
    mutate(gene_symbol = ifelse(is.na(gene_symbol), ensembl_gene_id, gene_symbol))

### SOURCE FUNCTIONS ###

## for plots
source(paste0(base_dir, "scripts/full_analysis_scripts/functions/plot_functions.R"))

## hypergeometric
source(paste0(base_dir, "scripts/full_analysis_scripts/functions/hypergeometric_function.R"))

### GROUP COVARIATES BASED ON TYPE ###

technical_covariates <- c("source", "pmi", "ph", "pmi_confidence", "max_rin", "max_rine", "mapped_percent", "five_prime_three_prime_bias",
                          "rna_extraction_batch", "library_batch", "gc_percent")
biological_covariates <- c("age_death", "sex_at_birth", "race", "marital_status", "education", "manner_death", "suicide", "brain_weight",
                           "height", "weight", "bmi")
drug_covariates <- c("smoker", "nicotine_cotinine", "alcohol", "sedative_hypnotic_anxiolitics", "opioids", "cannabinoids",
                     "major_stimulants_cocaine_included", "minor_stimulants", "anticholinergics", "antidepressants",
                     "anti_epileptics", "anti_histamines", "antipsychotics", "mood_stabilizers", "non_psychiatric",
                     "other_psychotropic_drug", "benzos")
