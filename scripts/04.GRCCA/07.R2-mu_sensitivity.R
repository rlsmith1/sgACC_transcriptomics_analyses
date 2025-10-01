
#------------------------------------------------------------------------------#
# GENE-LEVEL: Run GRCCA across mu values to demonstrate sensitivity of the main LV
# to the group penalty (mu ∈ {0.001 0.01 0.5 0.9 0.999}; quasi-logarithmic scale)
#------------------------------------------------------------------------------#


## Load data & functions
base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
figures_dir <- paste0(base_dir, "outputs/figures/response_to_reviewers/round2/")
#tables_dir <- paste0(base_dir, "outputs/tables/")
#analysis_objects_dir <- paste0(base_dir, "outputs/objects/")

source(paste0(base_dir, "scripts/load_WGCNA_res.R"))



# Load mu sensitivity results ---------------------------------------------


## Identify analysis directories for separate diagnosis runs
project_dir <- paste0("08Mar2024_GENES_qSVAgeSexRaceGC_sft3_minSize40_cutHeight0.98_ControlBDMDDSCZ_8MCA_regressBrainWeight_WITHGRAY/")
cca_dir <- paste0(base_dir, "RCCA_toolkit/GENES/", project_dir)


## Define variables altered for analysis
analysis_type <- "grcca" # rcca
analysis_framework <- "permutation" # holdout
VARx <- "0.1_1" # variance explained search grid (or set to specific value)
lambda <- "Nfeat"
mu <- c(0.001, 0.01, 0.5, 0.9, 0.999) # 0.1 is the original value
include_cmat <- "noCmat"


## Define analysis directory (variables altered in analysis with same X.mat and Y.mat)
analysis_dirs <- paste0(
    analysis_type, "_",
    analysis_framework, "_",
    "muSENS", "_", # turned on for mu sensitivity analysis; don't normally have this flag in
    "VARx", VARx, "_",
    "L2x", lambda, "_",
    "groupL2x", mu, "_",
    include_cmat
)


## Load GRCCA results for each analysis directory (across mu)
l_grcca_res <- 
    map(
        .x = analysis_dirs,
        .f = ~ f_load_cca_res(
            cca_directory = cca_dir,
            analysis_directory = .x,
            level = "gene",
            rename_covariates = TRUE
        )
    )
names(l_grcca_res) <- as.character(mu)



# Combine results into data frames for plotting ---------------------------


## Load original GRCCA results
load(paste0(base_dir, "objects/GRCCA_results.Rdata")) # df_lvs, df_results, df_y_res, df_x_res, df_rx_expr_cor


## Model results
df_results_muSENS <- map(
    .x = l_grcca_res,
    .f = ~ .x$model_results
) %>%
    list_rbind(names_to = "mu") %>%
    bind_rows(df_results %>% mutate(mu = "0.1"))


## Y results
df_y_res_muSENS <- map(
    .x = l_grcca_res,
    .f = ~ .x$y_res %>% 
        mutate_at(vars(c(weight, z_score, pearsons_r)), ~ -.x) # negate to align with original model
) %>%
    list_rbind(names_to = "mu") %>%
    bind_rows(df_y_res %>% mutate(mu = "0.1"))


## X results
df_x_res_muSENS <- map(
    .x = l_grcca_res,
    .f = ~ .x$x_res %>% 
        mutate_at(vars(c(weight, z_score, pearsons_r)), ~ -.x) # negate to align with original model
) %>%
    list_rbind(names_to = "mu") %>%
    left_join(df_ensembl_to_symbol) %>% # Add gene & module info
    left_join(df_modules_filt %>% dplyr::select(-c(mod_set, color)), by = join_by(ensembl_gene_id)) %>%
    bind_rows(df_x_res %>% mutate(mu = "0.1"))




# FIG XXA: mu sensitivity model results -----------------------------------


# Plot
df_results_muSENS %>% 
    
    # format
    group_by(mu) %>% 
    mutate(best_mod = ifelse(pval == min(pval), "black", "transparent")) %>%
    mutate(best_mod = ifelse(mu == 0.9 & set == 4, 0, best_mod)) %>%
    ungroup %>% 
    
    # layout
    ggplot(aes(x = factor(varx), y = -log10(pval))) +
    geom_hline(yintercept = -log10(0.005), lty = 2) +
    geom_segment(aes(y = 0, yend = -log10(pval), color = mu),
                 linewidth = 0.5, position = position_dodge(width = 0.6)) +
    geom_point(aes(size = correl, fill = mu, color = I(best_mod)),
               shape = 21, position = position_dodge(width = 0.6)) +
    
    # plot aesthetics
    scale_fill_nejm() +
    scale_color_nejm() +
    scale_size_continuous(range = c(0.5, 3), name = "X-Y LV cor") +
    labs(x = "X variance explained", y = "-log10(model P value)",
         title = "A | Mu sensitivity model comparison") +
    theme(
        legend.title = element_text(margin = margin(t = -5)),
        legend.margin = margin(l = -10)
    )


## Save plot
map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(figures_dir, "R3C2A.mu_sensitivity_model_results", .x), width = 6.5, height = 2.5)
)




# FIG XXB: Mu sensitivity covariate results -------------------------------


# Set covariate order for plotting
covariate_order <- c(
    "SCZ", "BD", "MDD",
    "Composite\ncase",
    paste0("Dim ", 1:8)
)


# Plot
df_y_res_muSENS %>%
    
    # format
    mutate(
        type = ifelse(str_detect(covariate, "Dim"), "tox", "main"),
        significant = ifelse(p_adj < 0.05 & abs(z_score) > 2, "black", "transparent"),
        covariate = factor(covariate, levels = rev(covariate_order))
    ) %>% 
    arrange(abs(pearsons_r)) %>% 
    
    # layout
    ggplot(aes(x = pearsons_r, y = covariate)) +
    geom_segment(aes(x = 0, xend = pearsons_r, color = mu),
                 position = position_dodge(width = 0.6), linewidth = 0.5) +
    geom_point(aes(size = -log10(p_adj + 10^-100), fill = mu, color = I(significant)), 
               position = position_dodge(width = 0.6), shape = 21) +
    geom_vline(xintercept = 0, color = "gray", linewidth = 0.5) +
    
    # Facet by covariate type
    facet_wrap(vars(type), scales = "free_y", ncol = 1) +
    force_panelsizes(rows = c(2, 8)) +
    
    # Plot aesthetics
    coord_flip() +
    scale_fill_nejm(guide = "none") +
    scale_color_nejm(guide = "none") +
    scale_size_continuous(range = c(0.5, 3), name = "-log10(FDR)") +
    #coord_cartesian(clip = "off") +
    labs(x = "Structure correlation (ry)", y = NULL,
         title = "B | Covariate results") +
    theme(
        strip.text = element_blank(),
        panel.spacing = unit(0.25, "cm"),
        legend.position = c(0.60, 0.20),
        legend.key.size = unit(0.25, "cm"),
        legend.title = element_text(margin = margin(t = -2))
        #legend.text = element_text(size = 5)
    )


## Save
map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(figures_dir, "R3C2B.mu_sensitivity_covariate_results", .x), width = 2.5, height = 3.7)
)




# FIG XXC: mu sensitivity gene results ------------------------------------


## Format to correlate x res at each mu sensitivity value with actual x res
df_x_res_muSENS_plot <- df_x_res_muSENS %>% 
    dplyr::select(mu, ensembl_gene_id, pearsons_r) %>% 
    filter(mu != "0.1") %>% 
    left_join(
        df_x_res %>% 
            dplyr::select(ensembl_gene_id, pearsons_r) %>% 
            dplyr::rename("actual_r" = "pearsons_r"),
        by = join_by(ensembl_gene_id)
    )

## Plot
df_x_res_muSENS_plot %>% 
    
    # layout
    ggplot(aes(x = pearsons_r, y = actual_r)) +
    geom_abline(color = "gray") +
    geom_point(size = 0.25) +
    geom_smooth(method = "lm", color = "maroon", linewidth = 0.5) +
    stat_cor(aes(label = after_stat(r.label)), size = 3) +
    facet_wrap(vars(mu), nrow = 1) +
    
    # aesthetics
    scale_y_continuous(limits = c(-1, 1), breaks = c(-1, 0, 1)) +
    scale_x_continuous(limits = c(-1, 1), breaks = c(-1, 0, 1)) +
    labs(x = "Mu sensitivity rx", y = "Original rx (mu = 0.1)",
         title = "C | Gene structure correlations across mu values") +
    theme_bw() +
    theme(
        plot.title = element_text(face = "bold", size = 11),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 9)
    )


## Save
map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(figures_dir, "R3C2C.mu_sensitivity_gene_results", .x), width = 6.5, height = 2)
)



# FIG XXD: Gene-level overlaps (after threshold for significance) ---------------------------------------------------


## Format data for upset plot
df_xres_muSENS_upset <- df_x_res_muSENS %>%
    
    # For each diagnosis, create logical column indicating if this gene is significant or not
    mutate(significant = ifelse(significant == 1, TRUE, FALSE)) %>% 
    arrange(mu) %>% 
    pivot_wider(id_cols = ensembl_gene_id, names_from = mu, values_from = significant) #%>% 
    
    # remove genes that weren't significant in any analysis
    filter(!(Actual == FALSE & BD == FALSE & MDD == FALSE & Composite == FALSE))

## Upset plot
upset(
    df_xres_muSENS_upset, 
    colnames(df_xres_muSENS_upset %>% dplyr::select(-ensembl_gene_id)),
    base_annotations = list(
        "Intersection size" = intersection_size(
            text = list(size = 3)
        )
    ),
    queries = list(
        upset_query(set = "0.001", color = pal_nejm()(1), fill = pal_nejm()(1)),
        upset_query(set = "0.01", color = pal_nejm()(2)[2], fill = pal_nejm()(2)[2]),
        upset_query(set = "0.1", color = pal_nejm()(3)[3], fill = pal_nejm()(3)[3]),
        upset_query(set = "0.5", color = pal_nejm()(4)[4], fill = pal_nejm()(4)[4]),
        upset_query(set = "0.9", color = pal_nejm()(5)[5], fill = pal_nejm()(5)[5]),
        upset_query(set = "0.999", color = pal_nejm()(6)[6], fill = pal_nejm()(6)[6])
    )
) +
    labs(x = NULL, title = "D | Shared significant genes across analyses") +
    theme(
        plot.title = element_text(face = "bold", size = 11),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 9)
    )


## Save
map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(figures_dir, "R3C3D.cross_disorder_gene_upset", .x), width = 6.5, height = 3)
)


## Stats: run hypergeometric tests across all combinations
l_sig_genes <- df_xdx_xres %>%
    filter(significant == 1) %>%
    group_split(Analysis) %>%
    map( ~ pull(.x, ensembl_gene_id)) %>%
    set_names(analysis_order)

df_hypergeometric_res <- expand_grid(
    analysis1 = analysis_order,
    analysis2 = analysis_order
) %>% 
    filter(analysis1 != analysis2) %>%
    group_by(grp = paste0(pmin(analysis1, analysis2), "_", pmax(analysis1, analysis2))) %>% 
    slice(1) %>%
    ungroup %>% 
    dplyr::select(-grp) %>% 
    mutate(
        hypergeometric_res = map2(
            .x = .$analysis1,
            .y = .$analysis2,
            .f = ~ f_hypergeometric(gene_list1 = l_sig_genes[[.x]], 
                                    gene_list2 = l_sig_genes[[.y]],
                                    gene_universe = ensembl_gene_ids
            )
        ),
        p_value = map(hypergeometric_res, ~ .x$p_value),
        odds_ratio = map(hypergeometric_res, ~ .x$odds_ratio)
    ) %>% 
    unnest(cols = c(p_value, odds_ratio))

df_hypergeometric_res %>% filter(p_value > 0.05)










