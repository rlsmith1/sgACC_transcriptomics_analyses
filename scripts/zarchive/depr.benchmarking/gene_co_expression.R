
########################################################################################

# Benchmark WGCNA and WTCNA results against the literature

########################################################################################

# setup -------------------------------------------------------------------

## LIBRARIES
library(tidyverse)
library(patchwork)
library(janitor)
library(RColorBrewer)
library(tidymodels)
library(tidytext)
library(stm)
library(ggwordcloud)
library(ggrepel)
library(readxl)
library(writexl)
library(clusterProfiler)

## SET DIRECTORIES
base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
prefix <- "08Mar2024_GENES_qSVAgeSexRaceGC"

## IDENTIFY MODULE SET OF INTEREST
soft_power <- 3
minimum_size <- 40
tree_cut_height <- 0.98

## LOAD DATA OBJECTS 
load(paste0(base_dir, "objects/hsapiens_genome_v110.RDS")) # df_hsapiens_genome
load(paste0(base_dir, "objects/14Nov2023_covariates.Rdata")) # df_covariates, df_covariates_numeric (generated in clean_covariates.R)
load(paste0(base_dir, "objects/", prefix, ".RDS")) # df_vsd_regress
load(paste0(base_dir, "objects/", prefix, "_SIGNED_SFT", soft_power, "_MODS.RDS")) # df_modules

## FILTER FOR MODULE SET OF INTEREST  
df_modules_filt <- df_modules %>% 
    filter(min_size == minimum_size & cut_height == tree_cut_height) %>% 
    unite("mod_set", c(sft, min_size, cut_height), sep = "_") %>% 
    arrange(mod_set, module) %>% 
    
    # rename modules to specify gene-level
    mutate(module = paste0("geneM", module) %>% factor(levels = paste0("geneM", levels(module)))
    )

## ENSEMBL ID TO GENE SYMBOL MAPPINGS
df_ensembl_to_symbol <- df_hsapiens_genome %>% 
    dplyr::select(gene_id, gene_name) %>% 
    distinct() %>% 
    dplyr::rename("gene_symbol" = "gene_name", "ensembl_gene_id" = "gene_id")

## PLOT THEME
theme_set(theme_bw() +
              theme(plot.title = element_text(size = 12),
                    axis.title = element_text(size = 12),
                    axis.text = element_text(size = 10),
                    strip.text = element_text(size = 10),
                    legend.title = element_text(size = 10),
                    legend.text = element_text(size = 10)
              )
)
module_colors <- df_modules_filt %>% dplyr::select(module, color) %>% distinct %>% deframe

# Link dx number to group abbreviation
group <- c(0, 1, 2, 3)
names(group) <- c("BD", "Control", "MDD", "SCZ")

# DRUG MCA RESULTS
load(paste0(base_dir, "objects/drug_MCA_results.Rdata")) # df_ind_loadings, df_var_loadings

# DEFINE COVARIATE TYPES
technical_covariates <- c("source", "pmi", "ph", "pmi_confidence", "max_rin", "max_rine", "mapped_percent", "five_prime_three_prime_bias",
                          "rna_extraction_batch", "library_batch", "gc_percent")
biological_covariates <- c("age_death", "sex_at_birth", "race", "marital_status", "education", "manner_death", "suicide", "brain_weight",
                           "height", "weight", "bmi")
drug_covariates <- c("smoker", "nicotine_cotinine", "alcohol", "sedative_hypnotic_anxiolitics", "opioids", "cannabinoids",
                     "major_stimulants_cocaine_included", "minor_stimulants", "anticholinergics", "antidepressants",
                     "anti_epileptics", "anti_histamines", "antipsychotics", "mood_stabilizers", "non_psychiatric",
                     "other_psychotropic_drug", "benzos")

# ADD DRUG MCA DIMENSIONS TO COVARIATE DF
df_covariates_mca <- df_covariates_numeric %>% 
    left_join(
        df_ind_loadings %>% 
            dplyr::rename_at(vars(contains("dim")), ~ str_replace(.x, "dim", "MC")),
        by = join_by(sample)
    ) %>% 
    
    # remove drug covariates (since we have the MCs)
    dplyr::select(-any_of(drug_covariates))


# run PCA on each module --------------------------------------------------


# RUN PCA FOR EACH MODULE
df_mods_pca <- df_modules_filt %>% 
    
    # combin with sample expression information
    left_join(df_vsd_regress %>% 
                  pivot_longer(2:ncol(.), names_to = "ensembl_gene_id", values_to = "resids"),
              by = join_by(ensembl_gene_id)
    ) %>% 
    
    # create transcript x sample expression matrix for each module
    pivot_wider(id_cols = c(ensembl_gene_id, module), 
                names_from = sample, 
                values_from = resids) %>% 
    group_by(module) %>% 
    nest() %>% 
    
    # run PCA on transcript x sample expression matrix
    mutate(pca = map(.x = data,
                     .f = ~ .x %>% 
                         as.data.frame %>%
                         column_to_rownames("ensembl_gene_id") %>% 
                         prcomp(center = TRUE, scale = TRUE)
    )) %>% 
    
    # combine PCA results
    mutate(sample_pcs = map(.x = pca,
                            .f = ~.x$rotation %>% 
                                as.data.frame %>% 
                                rownames_to_column("sample") %>% 
                                as_tibble() %>% 
                                clean_names() %>% 
                                pivot_longer(matches("*[0-9]"), 
                                             names_to = "pc", 
                                             values_to = "pc_score") %>% 
                                
                                left_join(
                                    
                                    summary(.x) %>% 
                                        .$importance %>% 
                                        as.data.frame %>% 
                                        rownames_to_column("metric") %>% 
                                        as_tibble() %>% 
                                        pivot_longer(2:ncol(.), 
                                                     values_to = "value", 
                                                     names_to = "pc") %>% 
                                        mutate(pc = tolower(pc)) %>% 
                                        pivot_wider(id_cols = pc, 
                                                    names_from = metric, 
                                                    values_from = value) %>% 
                                        clean_names()
                                    
                                ) %>% 
                                filter(proportion_of_variance > 0.01)
    )
    ) %>% 
    arrange(module) %>% 
    unnest(cols = c(sample_pcs)) %>% 
    dplyr::select(module, sample, pc, proportion_of_variance, pc_score) %>% 
    
    # combin sample PC scores with sample covariate information
    left_join(df_covariates_mca %>% 
                  pivot_longer(2:ncol(.), 
                               names_to = "covariate", 
                               values_to = "covariate_val"),
              by = "sample")

## EXTRACT PC1 AS EIGENGENE
df_kme <- df_mods_pca %>% 
    
    # pull eigengene value for each sample
    filter(pc == "pc1") %>% 
    pivot_wider(id_cols = c(module, sample, pc_score), 
                names_from = covariate, 
                values_from = covariate_val) %>% 
    dplyr::rename("kme" = "pc_score") %>% 
    mutate(dx = case_when(
        
        dx == 0 ~ "BD",
        dx == 1 ~ "Control",
        dx == 2 ~ "MDD",
        dx == 3 ~ "SCZ"
        
    ) %>% 
        factor(levels = c("Control", "BD", "MDD", "SCZ")),
    .before = 4)


# LINEAR REGRESSION OF MES AND DIAGNOSIS
df_lm_dx <- df_kme %>% 
    group_by(module) %>% 
    nest() %>% 
    mutate(
        data = map(
            .x = data,
            .f = ~ .x %>% 
                mutate(kme_resids = lm(kme ~ brain_weight) %>% residuals + median(kme))
        ),
        lm = map(
            .x = data,
            .f = ~ lm(kme_resids ~ dx + MC1 + MC2 + MC3 + MC4 + MC5 + MC6 + MC7 + MC8, data = .x)
        ),
        lm_sum = map(
            .x = lm,
            .f = ~ summary(.x)
        ),
        lm_res = map(
            .x = lm,
            .f = ~ tidy(.x)
        )
    ) %>% 
    unnest(cols = c(lm_res)) %>% 
    clean_names %>% 
    filter(term != "(Intercept)") %>% 
    mutate(p_adj = p.adjust(p_value, method = "BH"),
           term = str_remove(term, "dx"))



# pull genes that are part of that module & plot results ------------------


df_lm_dx %>% 
    filter(p_value < 0.05)

df_modules_filt %>% 
    filter(module == "geneM9")



