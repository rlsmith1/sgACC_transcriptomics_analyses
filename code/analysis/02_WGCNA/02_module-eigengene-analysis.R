
#==============================================================================#
# Identify WGCNA modules that are associated with psychiatric diagnosis
# according to canonical univariate module eigengene analysis
#==============================================================================#

## This script runs PCA on each module expression matrix
## And correlates all module PCs with all covariates (including MCA dimensions)
## Extract PC1 results for module eigengene results


## Setup
base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
source(paste0(base_dir, "scripts/load_WGCNA_res.R"))


# Data format ---------------------------------------------------------------

## Link dx to group number (in the numeric covariates tibble)
dx_group <- c(
    "BD" = 0,
    "Control" = 1,
    "MDD" = 2,
    "SCZ" = 3
)

## Combine known covariates with toxicology MCA dimensions
df_covariates_mca <- df_covariates_numeric %>% 
    left_join(
        df_ind_loadings %>% 
            dplyr::rename_at(vars(contains("dim")), ~ str_replace(.x, "dim", "MC")),
        by = join_by(sample)
    ) %>% 
    
    # remove drug covariates (since we have the MCs)
    dplyr::select(-any_of(drug_covariates)) %>% 
    
    # bring sample to column 1
    dplyr::select(sample, everything())



# Run PCA on each module -------------------------------------------------


## Run PCA on samples x gene expression matrix for each module
df_mods_pca <- df_modules_filt %>% 
    dplyr::select(module, ensembl_gene_id) %>% 
    
    # combine with sample expression information
    left_join(df_vsd_regress %>% 
                  pivot_longer(2:ncol(.), names_to = "ensembl_gene_id", values_to = "resids"),
              by = join_by(ensembl_gene_id)
    ) %>% 
    
    # create sample x gene expression matrix for each module
    group_by(module) %>% 
    nest() %>% 
    mutate(
        data = map(
            .x = data,
            .f = ~ .x %>% 
                pivot_wider(id_cols = sample, 
                            names_from = ensembl_gene_id, 
                            values_from = resids)
        )
    ) %>% 
    
    # run PCA on gene x sample expression matrix
    mutate(pca = map(
        .x = data,
        .f = ~ .x %>% 
            as.data.frame %>%
            column_to_rownames("sample") %>% 
            prcomp(center = TRUE, scale = TRUE)
    )
    )


## Extract sample scores in PC space
df_mods_pca_covariates <- df_mods_pca %>% 
    
    # extract PC scores
    mutate(sample_pcs = map(
        .x = pca,
        .f = ~.x$x %>% 
            as.data.frame %>% 
            rownames_to_column("sample") %>% 
            as_tibble() %>% 
            pivot_longer(matches("PC"), names_to = "PC", values_to = "score") %>% 
        
            # combine with PC variance explained
            left_join(
                summary(.x) %>% 
                    .$importance %>% 
                    t %>% as.data.frame %>% 
                    rownames_to_column("PC") %>% 
                    as_tibble() %>%
                    dplyr::rename("proportion_of_variance" = "Proportion of Variance") %>% 
                    dplyr::select(PC, proportion_of_variance),
                by = join_by(PC)
                
            )
    )
    ) %>% 
    unnest(cols = c(sample_pcs)) %>% 
    dplyr::select(module, sample, PC, proportion_of_variance, score) %>% 
    
    # combine sample PC scores with sample covariate information
    left_join(df_covariates_mca %>% 
                  pivot_longer(2:ncol(.), 
                               names_to = "covariate", 
                               values_to = "covariate_val"),
              by = join_by(sample)
    )


# Identify module-dx relationships across all module PCs ---------------------------------------------------------------------

df_lm_dx_allPCs <- df_mods_pca_covariates %>% 
    group_by(module, PC, proportion_of_variance) %>% 
    filter(covariate %in% c("dx", paste0("MC", 1:8))) %>% 
    pivot_wider(id_cols = c(module, PC, sample, proportion_of_variance, score), names_from = covariate, values_from = covariate_val) %>% 
    
    # convert dx to factor variable for regression analyssi
    mutate(dx = case_when(
        dx == 0 ~"BD",
        dx == 1 ~ "Control",
        dx == 2 ~ "MDD",
        dx == 3 ~ "SCZ"
    ) %>% factor(levels = names(dx_colors))
    ) %>% 
    nest() %>% 
    
    # run linear models on each PC of each module to identify all dx associations
    mutate(
        
        # run linear models and extract results
        lm = map(
            .x = data,
            .f = ~ lm(score ~ dx + MC1 + MC2 + MC3 + MC4 + MC5 + MC6 + MC7 + MC8, data = .x)
        ),
        lm_sum = map(
            .x = lm,
            .f = ~ summary(.x)
        ),
        lm_res = map(
            .x = lm,
            .f = ~ tidy(.x)
        ),
        
        # take residuals without diagnosis for plotting
        data = map(
            .x = data,
            .f = ~ .x %>% 
                mutate(pc_resids = lm(score ~ MC1 + MC2 + MC3 + MC4 + MC5 + MC6 + MC7 + MC8) %>% residuals + median(score))
        )
        
    ) %>% 
    unnest(cols = c(lm_res)) %>% 
    clean_names %>% 
    dplyr::filter(term != "(Intercept)") %>% 
    group_by(module, term) %>% 
    mutate(p_adj = p.adjust(p_value, method = "fdr"),
           term = str_remove(term, "dx")
    ) %>% 
    #group_by(term) %>% 
    #mutate(p_adj2 = p.adjust(p_value, method = "fdr")) %>% 
    dplyr::select(-c(lm, lm_sum))


df_lm_dx_allPCs %>% 
    filter(p_adj < 0.05 & proportion_of_variance >= 0.02 & module != "geneM0") %>% 
    print(n = nrow(.))



# Extract module eigengene data only -----------------------------------------------


df_kme <- df_lm_dx_allPCs %>% 
    filter(pc == "PC1")




# Save results -----------------------------------------------


## Remove data to save
df_lm_dx_allPCs <- df_lm_dx_allPCs %>% dplyr::select(-data)

# Save
save(df_mods_pca_covariates, df_lm_dx_allPCs,
     file = paste0(objects_dir, "module_pca_res.RDS")
)
save(df_kme,
     file = paste0(objects_dir, "kME_analysis_res.RDS")
)


