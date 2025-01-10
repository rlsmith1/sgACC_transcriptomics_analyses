
################################################################################

# GENE-LEVEL: Characterize RCCA results (gene weights & enrichments)

################################################################################

## Load data & functions
base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
source(paste0(base_dir, "scripts/full_analysis_scripts/genes/04.RCCA/A.load_RCCA_results.R"))
source(paste0(base_dir, "scripts/full_analysis_scripts/functions/module_hypergeometric_overlap.R")) # f_module_hypergeometric
figures_dir <- paste0(base_dir, "outputs/RCCA_res/")


# X-Y LV correlation ------------------------------------------------------

cor.test(df_lvs$lvx, df_lvs$lvy, method = "spearman")

# plot correlation
df_lvs %>% 
    ggplot(aes(x = lvx, y = lvy)) +
    geom_point(aes(color = dx), alpha = 0.7, size = 3) +
    geom_smooth(method = lm, se = FALSE, color = "black") +
    stat_cor() +
    scale_color_manual(values = dx_colors) +
    guides(color = guide_legend(title = "diagnosis")) +
    labs(x = "X latent variable (x matrix • x weights)",
         y = "Y latent variable (y matrix • y weights)",
         title = "Gene x-y latent variable correlation") +
    theme(legend.position = c(0.8, 0.15))

# plot X LV
df_lvs %>% 
    ggplot(aes(x = lvx)) +
    geom_density(aes(fill = dx), alpha = 0.7) +
    scale_fill_manual(values = dx_colors)


 
# Plot Y weights ---------------------------------------------------------------

df_y_res %>% 
    mutate(significant = ifelse(abs(z_score) >= 2 & p_adj < 0.05, "*", "")) %>% 
    
    ggplot(aes(x = reorder(covariate, weight), y = weight)) +
    geom_col(aes(fill = covariate, alpha = abs(weight)), color = "black") +
    geom_errorbar(aes(ymin = weight - sd, ymax = weight + sd), width = 0.2) +
    geom_text(aes(label = significant, y = 0), size = 6, vjust = 1.5) +
    geom_hline(aes(yintercept = 0), color = "black") +
    scale_fill_manual(values = c(dx_colors, rep("gray", 6))) +
    scale_alpha_continuous(range = c(0.1, 1)) +
    scale_x_discrete(labels = function(x) {str_replace(x, "_", " ")}) +
    ylim(c(-0.10, 0.10)) +
    labs(x = "", y = "weight", 
         title = "Covariate weights",
         #title = bquote(bold("A")~" | Covariate (y) coefficient weights"),
         caption = "* = abs(z-score) >= 2 & p-adj < 0.05; \nError bar shows +/- 1 sd") +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 0.9))
map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(figures_dir, "yres_weights", .x), width = 4, height = 3)
)
  

# Y structure correlations ---------------------------------------

df_y_res %>% 
    mutate(significant = ifelse(abs(z_score) >= 2 & p_adj < 0.05, "*", "")) %>% 
    
    ggplot(aes(x = pearsons_r, y = reorder(covariate, pearsons_r))) +
    geom_point(aes(size = -log10(p_adj), fill = pearsons_r), shape = 21) +
    geom_vline(xintercept = 0, lty = 2) +
    geom_text(aes(label = significant), size = 6, vjust = 0.75, hjust = -1.5, color = "black") +
    scale_fill_gradientn(colors = gene_weight_color_scale, limits = c(-1, 1)) +
    scale_size_continuous(range = c(2, 6)) +
    guides(fill = guide_colorbar(title = "Pearson's r"),
           size = guide_legend(title = "-log10(p-adj)")) +
    xlim(c(-1, 1)) +
    labs(y = "", x = "Pearson's r", 
         title = "Covariate structure correlations") +
    #title = bquote(bold("B")~" | Y data correlation with Y latent variable")) +
    theme(legend.position = c(0.15, 0.6), 
          legend.box = "vertical",
          legend.key.height = unit(0.35, "cm"),
          legend.key.width = unit(0.30, "cm")
    )
map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(figures_dir, "yres_struc_cor", .x), width = 4, height = 3)
)



# Correlate X structure correlations with expression values --------

### ILLUSTRATE STRUCTURE CORRELATION AND EXPRESSION IN SCZ
  
# normalized & regressed expression of GRCCA Pearson's r
df_vsd_regress %>% 
    pivot_longer(contains("ENSG"), names_to = "ensembl_gene_id", values_to = "expression_value") %>% 
    mutate(dx = case_when(
        str_detect(sample, "control") ~ "Control",
        str_detect(sample, "bipolar") ~ "BD",
        str_detect(sample, "mdd") ~ "MDD",
        str_detect(sample, "schizo") ~ "SCZ"
    ) %>% factor(levels = names(dx_colors))) %>% 
    group_by(ensembl_gene_id, dx) %>% 
    summarise(mean_expr = mean(expression_value)) %>% 
    left_join(df_x_res) %>% 
    
    # plot
    ggplot(aes(x = mean_expr, y = pearsons_r)) +
    geom_hline(yintercept = 0, lty = 2, color = "black") +
    geom_vline(xintercept = 0, lty = 2, color = "black") +
    geom_point(aes(color = pearsons_r, alpha = abs(pearsons_r))) +
    geom_smooth(method = "lm", color = "black", linewidth = 0.75) +
    stat_cor(aes(label = after_stat(r.label)), label.sep = "\n", label.y.npc = "bottom", vjust = 1) +
    facet_wrap(vars(dx), nrow = 1) +
    scale_color_gradientn(colors = gene_weight_color_scale) +
    labs(x = "Mean normalized & corrected expression", y = "Structure correlation (r)",
         title = "Gene structure correlations relationship with mean expression"
         #title = bquote(bold("C")~" | Correlations of GRCCA structure correlations with mean expression value across all genes")
    ) +
    theme(legend.position = "none")
map(
      .x = c(".png", ".pdf"),
      .f = ~ ggsave(paste0(figures_dir, "xres_struc_cor_mean_expr", .x), width = 12, height = 4)
  )
  
  

# Significant X res ---------------------------------------------------------------


### VECTOR OF ALL SIGNIFICANT GENES ###

df_x_res %>% 
    dplyr::filter(p_adj < 0.05 & abs(z_score) >= 2) %>% 
    
    ggplot(aes(x = pearsons_r, y = reorder(ensembl_gene_id, pearsons_r))) +
    geom_col(aes(fill = pearsons_r)) +
    geom_vline(xintercept = 0, lty = 2) +
    scale_fill_gradientn(colors = gene_weight_color_scale, limits = c(-1, 1)) +
    guides(fill = guide_colorbar(title = "Pearson's r")) +
    #scale_color_manual(values = c("*" = "black", "-" = "transparent")) +
    xlim(c(-1, 1)) +
    labs(y = "Gene", x = "Pearson's r", 
         title = "Significant gene structure correlations",
         #title = bquote(bold("D")~" | Significant gene structure correlations"),
         caption = "n = 1194"
    ) +
    theme(axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(),
          legend.position = c(0.25, 0.75), 
          legend.box = "horizontal",
          legend.key.size = unit(0.35, "cm")
    )
map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(figures_dir, "xres_struc_cor_significant", .x), # or grcca/withGray
                  width = 4, height = 3.15)
)

### STRUCTURE CORRELATIONS OF TOP N GENES ###

df_x_res %>% 
    dplyr::filter(p_adj < 0.05 & abs(z_score) >= 2) %>%
    top_n(n = 25, wt = abs(pearsons_r)) %>%
    mutate(direction = ifelse(pearsons_r < 0, "2neg", "1pos")) %>% 
    
    ggplot(aes(x = pearsons_r, y = reorder(gene_symbol, pearsons_r))) +
    geom_point(aes(size = -log10(p_adj), fill = pearsons_r), shape = 21) +
    geom_vline(xintercept = 0, lty = 2) +
    facet_wrap(vars(direction), ncol = 1, scales = "free_y") +
    force_panelsizes(rows = c(15, 10)) +
    scale_fill_gradientn(colors = gene_weight_color_scale, 
                         limits = c(-1, 1), guide = "none") +
    scale_size_continuous(range = c(2, 6)) +
    guides(size = guide_legend(title = "-log10(p-adj)")) +
    xlim(c(-1, 1)) +
    labs(y = "", x = "Pearson's r", 
         title ="Top 25 by abs(Pearson's r)"
    ) +
    theme(strip.text = element_blank(), 
          panel.spacing = unit(3, "lines"),
          legend.position = c(0.75, 0.20), 
          legend.box = "horizontal",
          legend.key.size = unit(0.35, "cm")
    )
map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(figures_dir, "xres_struc_cor_top25", .x), # or grcca/withGray
                  width = 4, height = 6)
)

  
### PLOT EXPRESSION LEVELS FOR SPECIFIC GENES ###

top_pos_gene <- df_x_res %>% dplyr::filter(pearsons_r > 0) %>% arrange(-abs(pearsons_r)) %>% head(1) %>% pull(ensembl_gene_id) # EZH1
top_neg_gene <- df_x_res %>% dplyr::filter(pearsons_r < 0) %>% arrange(-abs(pearsons_r)) %>% head(1) %>% pull(ensembl_gene_id) # CCDC80
df_vsd_regress %>% 
    dplyr::select(sample, all_of(c(top_neg_gene, top_pos_gene))) %>% 
    pivot_longer(contains("ENSG"), names_to = "ensembl_gene_id", values_to = "expr") %>% 
    left_join(df_ensembl_to_symbol) %>% 
    mutate(dx = case_when(
        str_detect(sample, "control") ~ "Control",
        str_detect(sample, "bipolar") ~ "BD",
        str_detect(sample, "mdd") ~ "MDD",
        str_detect(sample, "schizo") ~ "SCZ"
    ) %>% factor(levels = names(dx_colors))) %>% 
    mutate(gene_symbol = factor(gene_symbol, levels = c("EZH1", "CCDC80"))) %>% 
    
    # plot
    ggplot(aes(x = dx, y = expr)) +
    geom_hline(yintercept = 0, lty = 2, color = "black") +
    geom_point(aes(fill = dx), shape = 21, position = position_jitter(width = 0.2)) +
    geom_boxplot(aes(color = dx), fill = "transparent", outlier.shape = NA) +
    scale_fill_manual(values = dx_colors, guide = "none") +
    scale_color_manual(values = dx_colors, guide = "none") +
    facet_wrap(vars(gene_symbol), scales = "free_y", ncol = 1) +
    labs(x = "", y = "Normalized & corrected expression value",
         title = "Expression values for top positive and negative genes") +
    theme(panel.spacing = unit(5, "lines"),)
map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(figures_dir, "top_genes_expr_values", .x), # or grcca/withGray
                  width = 4.5, height = 6.5)
)

## DE stats
  df_vsd_regress %>% 
      dplyr::select(sample, all_of(c(top_neg_gene, top_pos_gene))) %>% 
      pivot_longer(contains("ENSG"), names_to = "ensembl_gene_id", values_to = "expr") %>% 
      left_join(df_ensembl_to_symbol) %>% 
      mutate(dx = case_when(
          str_detect(sample, "control") ~ "Control",
          str_detect(sample, "bipolar") ~ "BD",
          str_detect(sample, "mdd") ~ "MDD",
          str_detect(sample, "schizo") ~ "SCZ"
      ) %>% factor(levels = names(dx_colors))) %>% 
      dplyr::filter(gene_symbol == "CCDC80") %>% 
      aov(expr ~ dx, data = .) %>% TukeyHSD
  
  

# Module hypergeometric overrepresentation ----------------------------------------------------------------

sig_rcca_genes <- df_x_res %>% 
      dplyr::filter(p_adj < 0.05 & abs(z_score) >= 2) %>% 
      pull(ensembl_gene_id)
length(sig_rcca_genes) # n = 1194
  
## Hypergeometric test for significant overlap with modules
df_rcca_module_hypergeometric <- map_dfr(
    .x = names(module_colors),
    .f = ~ f_module_hypergeometric(
        gene_list = sig_rcca_genes, 
        gene_list_name = "", 
        wgcna_module = .x, 
        gene_universe = ensembl_ids
    ) %>% 
        dplyr::select(-benchmark_gene_list)
) %>% 
    mutate(p_adj = p.adjust(p_value, method = "fdr")) %>% 
    arrange(p_adj) %>% 
    mutate(module = factor(module, levels = names(module_colors))) %>% 
    
    # specify modules to label
    mutate(module_label = ifelse(p_adj < 0.01, paste0("n = ", overlap_n), ""))

## Plot
f_plot_module_hypergeometric(df_rcca_module_hypergeometric,
                             y_axis = "p-adj",
                             plot_title = "Module overrepresentation of RCCA genes",
                             module_text_size = 10,
                             subset_x_labs = TRUE,
                             include_module_labels = TRUE,
                             legend_position = "none",
                             p_value_threshold = 0.01,
                             include_threshold_label = FALSE
)

# save to project dir
map(
  .x = c(".png", ".pdf"),
  .f = ~ ggsave(paste0(figures_dir, "module_hypergeometric_res", .x), # or grcca/withGray
                width = 4, height = 3)
)



# Direction of significant genes in top modules ------------------

## Identify significant modules
sig_modules <- df_rcca_module_hypergeometric %>% dplyr::filter(p_adj < 0.01) %>% pull(module) %>% sort()

## Identify panel sizes (number of significant genes in each module)
panel_widths <- df_x_res %>% 
    dplyr::filter(module %in% sig_modules & ensembl_gene_id %in% sig_rcca_genes) %>% 
    count(module) %>% 
    pull(n)

## Plot significant genes in significant modules
df_x_res %>% 
    dplyr::filter(module %in% sig_modules & ensembl_gene_id %in% sig_rcca_genes) %>% 
    mutate(weight = ifelse(significant == 1, weight, 0)) %>% 

    ggplot(aes(x = reorder_within(ensembl_gene_id, within = module, by = pearsons_r), y = pearsons_r)) +
    geom_col(aes(fill = module, alpha = abs(pearsons_r))) +
    facet_wrap(vars(module), nrow = 1, scales = "free_x") +
    force_panelsizes(cols = panel_widths) +
    scale_fill_manual(values = module_colors) +
    scale_alpha_continuous(range = c(0.3, 1)) +
    labs(title = "Sig module gene structure correlations",
         x = "Gene", y = "Pearson's r") +
    guides(fill = guide_legend(nrow = 2)) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          strip.text = element_blank(),
          legend.position = "none"
    )
map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(figures_dir, "xres_sig_mods_gene_directions", .x), # or grcca/withGray
                  width = 4, height = 2.5)
)


# Significant module GO enrichments ----------------------------------------------------------------

## Run Gene Ontology enrichment on each of the significantly associated modules ##

sig_modules <- df_rcca_module_hypergeometric %>% dplyr::filter(p_adj < 0.01) %>% arrange(module) %>% pull(module) %>% as.character()
df_go_res <- map_dfr(
    .x = sig_modules,
    .f = ~ enrichGO(gene = df_modules_filt %>% dplyr::filter(module == .x) %>% pull(ensembl_gene_id),
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
) %>% 
    mutate(module = factor(module, levels = names(module_colors)))

# Plot
l_module_go_res <- map(
    .x = sig_modules,
    .f = ~ f_plot_go_by_ontology(df_go_res %>% dplyr::filter(module == .x), 
                                 pathway_text_width = 50,
                                 pathway_text_size = 10,
                                 plot_title = paste0(.x, " GO enrichment"),
                                 n_facet_rows = 1,
                                 legend_position = "none"
    )
)
wrap_plots(l_module_go_res, ncol = 1)

# save
map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(figures_dir, "xres_sig_mods_GO", .x), # or grcca/withGray
                  width = 12, height = 9)
)



# Gene set enrichment analyses (full RCCA gene list) ---------------------------------------------------
  
## Generate gene vector to run GSEA on 
ranked_gene_list <- df_x_res %>% 
    arrange(-pearsons_r) %>% 
    dplyr::select(ensembl_gene_id, pearsons_r) %>% 
    distinct() %>% 
    deframe

# Run GSEA for each ontology
set.seed(20240801)
df_gsea_res <- map_dfr(
    .x = c("BP", "CC", "MF"),
    .f = ~ gseGO(ranked_gene_list,
                 ont = .x,
                 keyType = "ENSEMBL",
                 OrgDb = "org.Hs.eg.db",
                 pvalueCutoff = 1,
                 eps = 1e-300
    ) %>% 
        as_tibble() %>% 
        clean_names %>% 
        mutate(ontology = .x, .before = 1)
)

## Plot
f_plot_gsea_by_ontology(df_gsea_res, 
                        pathway_text_size = 3,
                        pathway_text_width = 30,
                        pathway_text_overlaps = 30,
                        plot_title = "RCCA GSEA enrichments",
                        legend_position = c(0.9, 0.1)
)

## Save
map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(figures_dir, "RCCA_GSEA_enrichment", .x), width = 6, height = 10)
)

  # # export table
  # write_xlsx(list("Biological process" = as.data.frame(df_gsea_res %>% filter(ontology == "BP") %>% dplyr::select(-ontology)),
  #                 "Cellular component" = as.data.frame(df_gsea_res %>% filter(ontology == "CC") %>% dplyr::select(-ontology)),
  #                 "Molecular function" = as.data.frame(df_gsea_res %>% filter(ontology == "MF") %>% dplyr::select(-ontology))
  # ),
  # paste0(base_dir, "outputs/tables/for_manuscript/TableS4_GENES_GSEA_res.xlsx")
  # )
 
 