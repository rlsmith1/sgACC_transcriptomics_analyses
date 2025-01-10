
########################################################################################

# Identify WGCNA modules that are associated with psychiatric diagnosis
# according to canonical univariate module eigengene analysis

########################################################################################

## Load data & functions
base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
figures_dir <- paste0(base_dir, "outputs/figures/supplement/")
#tables_dir <- paste0(base_dir, "outputs/tables/")
analysis_objects_dir <- paste0(base_dir, "outputs/objects/")

source(paste0(base_dir, "scripts/load_WGCNA_res.R"))


# Data format ---------------------------------------------------------------

## Link dx to group number (in the numeric covariates tibble)
dx_group <- c(
    "BD" = 0,
    "Control" = 1,
    "MDD" = 2,
    "SCZ" = 3
)

## Combine known covariates with drug MCA dimensions
df_covariates_mca <- df_covariates_numeric %>% 
    left_join(
        df_ind_loadings %>% 
            dplyr::rename_at(vars(contains("dim")), ~ str_replace(.x, "dim", "MC")),
        by = join_by(sample)
    ) %>% 
    
    # remove drug covariates (since we have the MCs)
    dplyr::select(-any_of(drug_covariates))


# Run PCA on each module -------------------------------------------------

df_mods_pca <- df_modules_filt %>% 
    
    # combine with sample expression information
    left_join(df_vsd_regress %>% 
                  pivot_longer(2:ncol(.), names_to = "ensembl_gene_id", values_to = "resids"),
              by = join_by(ensembl_gene_id)
    ) %>% 
    
    # create gene x sample expression matrix for each module
    pivot_wider(id_cols = c(ensembl_gene_id, module), 
                names_from = sample, 
                values_from = resids) %>% 
    group_by(module) %>% 
    nest() %>% 
    
    # run PCA on gene x sample expression matrix
    mutate(pca = map(.x = data,
                     .f = ~ .x %>% 
                         as.data.frame %>%
                         column_to_rownames("ensembl_gene_id") %>% 
                         prcomp()
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
                                        clean_names(),
                                    by = join_by(pc)
                                    
                                ) %>% 
                                dplyr::filter(proportion_of_variance > 0.01)
    )
    ) %>% 
    arrange(module) %>% 
    unnest(cols = c(sample_pcs)) %>% 
    dplyr::select(module, sample, pc, proportion_of_variance, pc_score) %>% 
    
    # combine sample PC scores with sample covariate information
    left_join(df_covariates_mca %>% 
                  pivot_longer(2:ncol(.), 
                               names_to = "covariate", 
                               values_to = "covariate_val"),
              by = join_by(sample)
    )


# Extract module eigengenes -----------------------------------------------

df_kme <- df_mods_pca %>% 
    dplyr::filter(pc == "pc1") %>% 
    pivot_wider(id_cols = c(module, sample, pc_score), 
                names_from = covariate, 
                values_from = covariate_val) %>% 
    dplyr::rename("kme" = "pc_score", "value" = "dx") %>% 
    left_join(enframe(dx_group, name = "dx"), by = join_by(value)) %>% 
    mutate(dx = factor(dx, levels = names(dx_colors))) %>% 
    dplyr::select(module, sample, kme, dx, everything())


# Identify significant modules using linear regression --------------------

## For each module, run a linear model on module eigengene with drug MCs as covariates
df_lm_dx <- df_kme %>% 
    group_by(module) %>% 
    nest() %>% 
    mutate(
        
        # run linear models and extract results
        lm = map(
            .x = data,
            .f = ~ lm(kme ~ dx + MC1 + MC2 + MC3 + MC4 + MC5 + MC6 + MC7 + MC8, data = .x)
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
                mutate(kme_resids = lm(kme ~ MC1 + MC2 + MC3 + MC4 + MC5 + MC6 + MC7 + MC8) %>% residuals + median(kme))
        )
        
    ) %>% 
    unnest(cols = c(lm_res)) %>% 
    clean_names %>% 
    dplyr::filter(term != "(Intercept)") %>% 
    group_by(term) %>% 
    mutate(p_adj = p.adjust(p_value, method = "fdr"),
           term = str_remove(term, "dx")
    ) %>% 
    dplyr::select(-c(lm, lm_sum))

## Identify significant associations
df_lm_dx %>% 
    dplyr::filter(p_value < 0.05)

## Identify modules significantly associated with SCZ specifically
sig_modules <- df_lm_dx %>% 
    dplyr::filter(p_value < 0.05 & term == "SCZ") %>% 
    pull(module)



# Save results for plotting -----------------------------------------------


# Figures
save(df_lm_dx,
     file = paste0(analysis_objects_dir, "kME_analysis_res.RDS")
)

# Downstream analyses
save(df_lm_dx,
     file = paste0(base_dir, "objects/kME_analysis_res.RDS")
)


