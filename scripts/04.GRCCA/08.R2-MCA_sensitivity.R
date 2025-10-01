#------------------------------------------------------------------------------#
# GENE-LEVEL: Run GRCCA using a different approach to imputing missing values in
# the toxicology data (see 01.preprocessing/05b.R2-MCAv2.R)
#------------------------------------------------------------------------------#


## Load data & functions
base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
figures_dir <- paste0(base_dir, "outputs/figures/response_to_reviewers/round2/")
analysis_objects_dir <- paste0(base_dir, "outputs/objects/")

#tables_dir <- paste0(base_dir, "outputs/tables/")

source(paste0(base_dir, "scripts/load_WGCNA_res.R"))



# Load mu sensitivity results ---------------------------------------------

## Identify analysis directories for separate diagnosis runs
project_dir <- paste0("2025.09.30_GENES_qSVAgeSexRaceGC_sft3_minSize40_cutHeight0.98_ControlBDMDDSCZ_8MCA_modeImpute/")
cca_dir <- paste0(base_dir, "RCCA_toolkit/GENES/", project_dir)

## Define variables altered for analysis
analysis_type <- "grcca"
VARx <- "0.1_1"
mu <- 0.1


## Load GRCCA results the GRCCA analysis run with the new MCA data
l_grcca_res_mcaImpute <- f_load_cca_res(
    cca_directory = cca_dir,
    level = "gene",
    analysis_type = analysis_type,
    mu = mu,
    VARx = VARx,
    include_Cmat = FALSE,
    rename_covariates = FALSE,
    all_effects = FALSE
)


## Load original GRCCA results
load(paste0(analysis_objects_dir, "GRCCA_results.Rdata"))




# FIG SXXA: Compare model results -----------------------------------------


## Combine new and original GRCCA results
df_figSXXA <- df_results %>% 
    mutate(Analysis = "Original", .before = 1) %>% 
    bind_rows(
        l_grcca_res_mcaImpute$model_results %>% 
            mutate(Analysis = "MCA sensitivity", .before = 1)
    ) %>% 
    
    # Format for plotting
    mutate(Analysis = factor(Analysis, levels = c("Original", "MCA sensitivity"))) %>%
    group_by(Analysis) %>% 
    mutate(best_mod = ifelse(pval == min(pval), "black", "transparent")) %>% 
    ungroup


# Plot
df_figSXXA %>% 
    ggplot(aes(x = factor(varx), y = -log10(pval))) +
    geom_hline(yintercept = -log10(0.005), lty = 2) +
    geom_segment(aes(y = 0, yend = -log10(pval), color = Analysis),
                 position = position_dodge(width = 0.5)) +
    geom_point(aes(size = correl, fill = Analysis, color = I(best_mod)),
               shape = 21, position = position_dodge(width = 0.5)) +
    
    # plot aesthetics
    #scale_fill_nejm() +
    #scale_color_nejm() +
    scale_size_continuous(range = c(0.5, 4), name = "X-Y LV cor") +
    labs(x = "X variance explained", y = "-log10(model P value)",
         title = "A | MCA sensitivity model comparison") +
    theme(
        legend.title = element_text(margin = margin(t = -5), size = 7),
        legend.text = element_text(size = 6),
        legend.position = c(0.1, 0.75),
        legend.key.size = unit(0.3, "cm")
    )


## Save plot
map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(figures_dir, "R3C1A.MCA_sensitivity_model_results", .x), width = 3.25, height = 2.5)
)



# FIG SXXB: Compare covariate results -------------------------------------


covariate_order <- c(
    "SCZ", "BD", "MDD",
    paste0("Dim ", 1:8)
)
df_figSXXB <- df_y_res %>% 
    mutate(Analysis = "Original", .before = 1) %>% 
    bind_rows(
        l_grcca_res_mcaImpute$y_res %>% 
            mutate(Analysis = "MCA sensitivity", .before = 1)
    ) %>% 
    mutate(Analysis = factor(Analysis, levels = c("Original", "MCA sensitivity"))) %>%
    
    # format for plotting
    mutate(
        significant = ifelse(p_adj < 0.05 & abs(z_score) > 2, "black", "transparent"),
        covariate = factor(covariate, levels = rev(covariate_order))
    ) %>% 
    
    # negate ry for Dim 1-3 in the MCA sensitivity analysis (to get everything in the same direction)
    mutate(
        pearsons_r = ifelse(
            Analysis == "MCA sensitivity" & covariate %in% paste0("Dim ", 1:3),
            -pearsons_r,
            pearsons_r)
    ) %>% 
    mutate(type = ifelse(str_detect(covariate, "Dim"), "tox", "main"))



# Plot
df_figSXXB %>% 
    ggplot(aes(x = pearsons_r, y = covariate)) +
    geom_segment(aes(x = 0, xend = pearsons_r, color = Analysis),
                 position = position_dodge(width = 0.75)) +
    geom_point(aes(size = -log10(p_adj + 10^-100), fill = Analysis, color = I(significant)), 
               position = position_dodge(width = 0.75), shape = 21) +
    geom_vline(xintercept = 0, color = "gray", linewidth = 0.5) +
    
    # Facet by covariate type
    facet_wrap(vars(type), scales = "free_y", ncol = 1) +
    force_panelsizes(rows = c(4, 8)) +
    
    # Plot aesthetics
    #scale_fill_nejm(guide = "none") +
    #scale_color_nejm(guide = "none") +
    scale_size_continuous(range = c(0.5, 4), name = "-log10(FDR)") +
    coord_cartesian(clip = "off") +
    labs(x = "Structure correlation (ry)", y = NULL,
         title = "B | Covariate results") +
    theme(
        strip.text = element_blank(),
        panel.spacing = unit(0.15, "cm"),
        legend.position = c(0.60, 0.20),
        legend.key.size = unit(0.25, "cm"),
        legend.title = element_text(margin = margin(t = -2))
        #legend.text = element_text(size = 5)
    )


## Save
map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(figures_dir, "R3C1B.MCA_sensitivity_covariate_results", .x), width = 3.15, height = 2.5)
)




# FIG SXXC: Compare gene correlations ----------------------------------------------------


## Combine new and original GRCCA results
df_figSXXC <- df_x_res %>% 
    mutate(Analysis = "Original", .before = 1) %>% 
    bind_rows(
        l_grcca_res_mcaImpute$x_res %>% 
            mutate(Analysis = "MCA sensitivity", .before = 1) %>% 
            left_join(df_ensembl_to_symbol) %>% 
            left_join(df_modules_filt %>% dplyr::select(-c(mod_set, color)), by = join_by(ensembl_gene_id))
    ) %>% 
    
    # Format for plotting
    mutate(Analysis = factor(Analysis, levels = c("Original", "MCA sensitivity")))


## Plot
df_figSXXC %>% 
    pivot_wider(id_cols = ensembl_gene_id, names_from = Analysis, values_from = pearsons_r) %>% 
    
    # Layout
    ggplot(aes(x = Original, y = `MCA sensitivity`)) +
    geom_vline(xintercept = 0, color = "gray") +
    geom_hline(yintercept = 0, color = "gray") +
    geom_point(aes(color = Original + `MCA sensitivity`), size = 0.25) +
    geom_smooth(method = "lm", color = "black", linewidth = 0.5) +
    stat_cor(aes(label = after_stat(r.label))) +
    
    # Aesthetics
    scale_color_gradientn(colors = rev(brewer.rdbu(100)), limits = c(-1.1, 1.1), guide = "none") +
    labs(title = "C | Gene rx correlation")


## Save
map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(figures_dir, "R3C1C.MCA_sensitivity_gene_results", .x), width = 3.15, height = 2.5)
)



# Run GSEA enrichments -----------------------------------------------------


### Get ranked gene list ###
my_grcca_ranked <- l_grcca_res_mcaImpute$x_res %>% 
    arrange(-pearsons_r) %>% 
    dplyr::select(ensembl_gene_id, pearsons_r) %>% 
    deframe


### Cell-type GSEA ###
df_grcca_fgsea_cell <- fgsea(pathways = cell_types, 
                             stats = my_grcca_ranked, 
                             eps = 0
)  %>% 
    as_tibble() %>% 
    dplyr::rename("cell_type" = "pathway") %>% 
    arrange(-abs(NES))


### GO GSEA ###
df_grcca_fgsea_go <- fgsea(pathways = go_pathways, 
                           stats = my_grcca_ranked, 
                           eps = 0
) %>% 
    as_tibble() %>% 
    dplyr::rename("go_id" = "pathway") %>% 
    arrange(-abs(NES)) %>% 
    left_join(df_go_terms) %>% # add names of GO terms
    dplyr::select(go_id, term, ontology, everything())


### Benchmark lists ###
df_grcca_fgsea_benchmark <- fgsea(pathways = benchmarking_lists, 
                                  stats = my_grcca_ranked, 
                                  eps = 0
) %>% 
    as_tibble() %>% 
    dplyr::rename("benchmark_list" = "pathway") %>% 
    arrange(-abs(NES))



# FIG SXXD: Risk gene GSEA ------------------------------------------------


## Prep
df_figSXXD <- df_grcca_fgsea_benchmark %>% 
    mutate(class = ifelse(str_detect(benchmark_list, "SCZ"), "SCZ", "other") %>% factor(levels = c("SCZ", "other"))) %>% 
    mutate(benchmark_list = str_remove_all(benchmark_list, "SCZ|\\)|\\(") %>% 
               str_trim %>% 
               str_replace(", ", ",\n")
    ) %>% 
    mutate(label_color = ifelse(padj < 0.04, "white", "black")) %>% 
    mutate(box_color = ifelse(padj < 0.04, "black", "transparent")) %>% 
    mutate(label = round(NES, 2) %>% format(digits = 3),
           label = ifelse(padj < 0.04, paste0(label, "*"), label)
    )
#mutate(sig = ifelse(padj < 0.05, scientific(padj, digits = 3), round(padj, 2) %>% format(digits = 3)))
# mutate(sig = case_when(padj < 0.001 ~ "***", padj < 0.01 ~ "**", padj < 0.05 ~ "*", TRUE ~ ""))


## Plot
df_figSXXD %>% 
    ggplot(aes(x = benchmark_list, y = 1)) +
    geom_tile(aes(fill = NES, color = I(box_color)), width = 0.95, height = 0.95, linewidth = 0.75) +
    geom_text(aes(label = label, color = I(label_color)), size = 3) +
    facet_wrap(vars(class), scales = "free", nrow = 2) +
    scale_fill_gradientn(colors = rev(brewer.rdbu(100)[1:50]), limits = c(0.99, 1.50)) +
    #scale_fill_gradientn(colors = rev(brewer.rdbu(100)), limits = c(-1.58, 1.58)) +
    guides(fill = guide_colorbar(title.hjust = 0.5, title.position = "left")) +
    labs(x = NULL, y = NULL,
         title = "D | Risk gene enrichment") +
    theme(axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          axis.text.x = element_text(margin = margin(t = -3)),
          strip.background = element_rect(fill = "white", color = "black"),
          legend.margin = margin(r = -15),
          legend.position = "left",
          legend.justification = "top",
          legend.title = element_text(angle = 90),
          legend.key.height = unit(0.75, "cm"),
          legend.key.width = unit(0.15, "cm")
    )


## Save
map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(figures_dir, "R3C1D.MCA_sensitivity_risk_gene_enrichment", .x), width = 3.15, height = 2.25)
)







