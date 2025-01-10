
################################################################################

# Identify gene-level modules that are associated with psychiatric diagnosis
# according to canonical univariate module eigengene analysis

################################################################################


# Setup -------------------------------------------------------------------

## Set directories
base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
figures_dir <- paste0(base_dir, "outputs/WGCNA_res/")

## Load data and functions
source(paste0(base_dir, "scripts/full_analysis_scripts/genes/load_WGCNA_res.R"))

## Link dx to group number (in the numeric covariates tibble)
group <- c(0, 1, 2, 3)
names(group) <- c("BD", "Control", "MDD", "SCZ")


# Data format: combine known covariates with drug MCA dimensions ---------

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
    left_join(enframe(group, name = "dx"), by = join_by(value)) %>% 
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
    dplyr::select(-c(data, lm, lm_sum))

## Identify significant associations
df_lm_dx %>% 
    dplyr::filter(p_value < 0.05)

## Identify modules significantly associated with SCZ specifically
sig_modules <- df_lm_dx %>% 
    dplyr::filter(p_value < 0.05 & term == "SCZ") %>% 
    pull(module)

## Save results & write as table
save(df_lm_dx, file = paste0(base_dir, "objects/kME_analysis_res.RDS"))
write_xlsx(df_lm_dx %>% dplyr::filter(p_value < 0.05),
           path = paste0(figures_dir, "sig_kME_analysis_res.xlsx")
)


# Plot significant modules & eigengene distributions ----------------------

df_lm_dx %>% 
    filter(p_value < 0.05 & term == "SCZ") %>% 
    unnest(cols = c(data)) %>% 
    mutate(label = case_when(
        dx == "SCZ" & p_value < 0.05 & p_adj > 0.05 ~ "*",
        dx == "SCZ" & p_value < 0.05 & p_adj < 0.05 ~ "**",
        TRUE ~ ""
    )
    ) %>% 
    
    ggplot(aes(x = kme_resids, y = dx, color = dx, fill = dx)) +
    geom_point(position = position_jitter(width = 0.1), size = 0.75) +
    geom_boxplot(alpha = 0.1, outlier.shape = NA) +
    geom_text(aes(label = label, x = -Inf), color = "black", size = 8,
              vjust = 0.75, hjust = -0.1, check_overlap = TRUE) +
    facet_wrap(vars(module), ncol = 1, scales = "free_x") +
    scale_fill_manual(values = dx_colors) +
    scale_color_manual(values = dx_colors) +
    labs(x = "kME (corrected for covariates)", y = "",
         caption = "* = p < 0.05; ** = FDR < 0.05",
         title = "Significant module eigengene-\ndiagnosis associations"
    ) +
    theme(legend.position = "none",
          axis.text.y = element_text(size = 10)
    )

map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(figures_dir, "WGCNA_eigengene_res", .x), width = 3.5, height = 6)
)


# Functional enrichment of WGCNA modules -------------------------------------------

### GENE ONTOLOGY ###

# Run Gene Ontology enrichment on each of the significantly associated modules

df_go_res <- map_dfr(
    .x = sig_modules,
    .f = ~ enrichGO(gene = df_modules_filt %>% filter(module == .x) %>% pull(ensembl_gene_id),
                    OrgDb = "org.Hs.eg.db",
                    universe = ensembl_ids,
                    keyType = "ENSEMBL",
                    ont = "ALL",
                    pAdjustMethod = "fdr",
                    pvalueCutoff = 1,
                    qvalueCutoff = 1,
                    readable = TRUE
    ) %>% 
        as_tibble() %>% 
        mutate(module = .x, .before = 1)
)


# Plot
l_module_go_res <- map(
    .x = sig_modules,
    .f = ~ f_plot_go_by_ontology(df_go_res %>% filter(module == .x), 
                                 pathway_text_width = 50,
                                 pathway_text_size = 10,
                                 plot_title = paste0(.x, " GO enrichment"),
                                 n_facet_rows = 1,
                                 legend_position = "none"
    )
)
wrap_plots(l_module_go_res, ncol = 1)

## Save plot
map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(figures_dir, "WGCNA_GO_enrichment", .x), width = 12, height = 9)
)

## Save table
write_xlsx(df_go_res %>% filter(pvalue < 0.05),
           path = paste0(figures_dir, "SCZ_module_GO_enrichments.xlsx")
)

