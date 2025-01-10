
###################################################################

# Updated DGE analyis including all covariates used in WGCNA/GRCCA

###################################################################

## Set directories & load data
base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
figures_dir <- paste0(base_dir, "outputs/DGE_res/")
source(paste0(base_dir, "scripts/full_analysis_scripts/genes/setup.R"))


# Convert count data to matrix for DESeq2 --------------------------------------

m_vsd_regress <- df_vsd_regress %>% 
    dplyr::filter(str_detect(sample, "schizo|control")) %>% # subset for only schizophrenia-control comparisons
    column_to_rownames("sample") %>% 
    t() + 10 # add pseudo-count of 10 so there are no negative values in the matrix
m_vsd_regress <- round(m_vsd_regress*10^5) # multiply by 10000 and round to convert to integers but maintain precision


# Generate covariates matrix ----------------------------------------------

## Load drug MCA results
load(paste0(base_dir, "objects/drug_MCA_results.Rdata")) # df_ind_loadings, df_var_loadings

## Combine drug MCs with known covariates
df_covariates_mca <- df_covariates %>% 
    left_join(
        df_ind_loadings %>% 
            dplyr::rename_at(vars(contains("dim")), ~ str_replace(.x, "dim", "MC")),
        by = join_by(sample)
    )

## Convert to matrix
n_mcs <- 8
paste0("MC", 1:n_mcs)
m_covariates <- df_covariates_mca %>% 
    dplyr::filter(dx %in% c("SCZ", "Control")) %>% # SCZ vs controls only
    mutate(dx = factor(dx, levels = c("Control", "SCZ"))) %>% 
    dplyr::select(sample, dx, all_of(paste0("MC", 1:n_mcs))) %>% 
    column_to_rownames("sample")


# Run DESeq ---------------------------------------------------------------

## Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = m_vsd_regress, 
                              colData = m_covariates, 
                              design = ~ dx + 
                                  MC1 + MC2 + MC3 + MC4 + MC5 + MC6 + MC7 + MC8
                              )

## Run DESeq
de_res <- DESeq(dds)  

## Create tibble out of results
df_de_res <- results(de_res, name = "dx_SCZ_vs_Control") %>% 
    as.data.frame() %>% 
    rownames_to_column("ensembl_gene_id") %>% 
    as_tibble() %>% 
    clean_names %>% 
    left_join(df_ensembl_to_symbol) %>% 
    dplyr::select(ensembl_gene_id, gene_symbol, everything()) %>% 
    arrange(-abs(log2fold_change)) %>% 
    mutate(base_mean = (base_mean/10^5) - 10) %>% 
    arrange(-abs(stat)) %>% 
    mutate(padj = p.adjust(pvalue, method = "BH"))
save(df_de_res, file = paste0(base_dir, "objects/DE_results.RDS"))

## Export as table
write_xlsx(df_de_res, path = paste0(figures_dir, "DEGs.xlsx"))


# Plot DEGs (volcano) ---------------------------------------------------------

## Count
df_de_res %>% dplyr::filter(padj < 0.05) # n = 5
df_de_res %>% dplyr::filter(pvalue < 0.05) # n = 1690

## Format df for plot
df_de_volcano <- df_de_res %>% 
    mutate(gene_symbol = ifelse(is.na(gene_symbol), ensembl_gene_id, gene_symbol)) %>% 
    mutate(fill = ifelse(pvalue > 0.05, NA_real_, log2fold_change)) %>% 
    mutate(significant = ifelse(padj < 0.05, "yes", "no")) %>% 
    mutate(label = ifelse(-log10(pvalue) > 4 | (pvalue < 0.05 & abs(log2fold_change) > 0.1), gene_symbol, NA))
summary(df_de_volcano$log2fold_change)

## Volcano plot of results
df_de_volcano %>% 
    ggplot(aes(x = log2fold_change, y = -log10(pvalue))) +
    geom_point(aes(fill = fill, color = significant), shape = 21, size = 1) +
    geom_hline(aes(yintercept = -log10(0.05)), lty = 2, color = "black") +
    geom_text_repel(aes(label = label), min.segment.length = 0, size = 2.5) +
    scale_fill_gradientn(colors = gene_weight_color_scale, 
                         values = rescale(c(-0.17, 0, 0.17)),
                         limits = c(-0.17, 0.17),
                         na.value = "lightgray", guide = "none") +
    scale_color_manual(values = c("yes" = "black", "no" = "transparent"), na.value = "transparent") +
    guides(color = guide_legend(title = "FDR < 0.05", title.vjust = -1.5)) +
    labs(x = "log2(fold change)", y = "-log10(p-value)",
         title = "DEG volcano plot") +
    theme(legend.position = c(0.01, 0.91),
          legend.spacing.x = unit(0.1, "mm"),
          legend.background = element_blank(),
          legend.box.background = element_rect(colour = "black")
    )
map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(figures_dir, "DEGs_volcano", .x), width = 4, height = 3)
)


# Compare with Nirmala's results ------------------------------------------

## Combine with Akula et al 2021 results
df_de_akula_res <- df_de_res %>% 
    mutate(analysis = "current", .before = 1) %>% 
    
    bind_rows(df_akula_res %>% 
                  dplyr::filter(level == "gene" & comparison == "schizo_ctrl") %>% 
                  dplyr::select(-c(level, comparison)) %>% 
                  mutate(analysis = "akula") %>% 
                  dplyr::rename("ensembl_gene_id" = "id", "gene_symbol" = "symbol")
    )

## Venn diagram (overlap with Nirmala's results)
my_degs <- df_de_res %>% dplyr::filter(pvalue < 0.05) %>% pull(ensembl_gene_id)
akula_degs <- df_akula_res %>% 
    dplyr::filter(level == "gene" & comparison == "schizo_ctrl" & pvalue < 0.05) %>% 
    pull(id)
hypergeometric_overlap <- f_hypergeometric(my_degs, akula_degs, gene_universe = ensembl_ids)
length(hypergeometric_overlap$overlap_genes) # n = 454; p << 0.0001

## Create df for plot
df_correlation_plot <- df_de_akula_res %>% 
    
    # indicate if gene is significant in the current analysis
    mutate(significant = ifelse(ensembl_gene_id %in% my_degs, "yes", "no")) %>% 
    
    # pivot wider to plot
    pivot_wider(id_cols = c(ensembl_gene_id, gene_symbol, significant), names_from = analysis, values_from = log2fold_change) %>% 
    dplyr::filter(!is.na(akula)) %>% 
    
    # normalize values to compare across analyses
    mutate(akula = scale(akula), current = scale(current)) %>% 
    
    # add gene labels
    mutate(label = ifelse(abs(current) > 1 & abs(akula) > 1.1, gene_symbol, "")) %>% 
    
    # arrange to plot
    arrange(-abs(current))

## Plot relationship between L2FC in each analysis
df_correlation_plot %>% 
    ggplot(aes(x = akula, y = current)) +
    geom_point(aes(fill = akula, color = significant), shape = 21, size = 1) +
    geom_hline(aes(yintercept = 0), lty = 2, color = "grey") +
    geom_vline(aes(xintercept = 0), lty = 2, color = "grey") +
    geom_smooth(method = "lm", color = "black", linewidth = 0.75) +
    geom_abline(lty = 2) +
    #geom_text(x = 2.1, y = 2.6, label = "y = x", check_overlap = TRUE) +
    geom_text_repel(aes(label = label), min.segment.length = 0, max.overlaps = 15, size = 2.5) +
    stat_cor() +
    #stat_regline_equation(label.y.npc = "top", vjust = 3) +
    scale_fill_gradientn(colors = gene_weight_color_scale, guide = "none") +
    scale_color_manual(values = c("transparent", "black")) +
    guides(color = guide_legend(title = "FDR < 0.05", title.vjust = -1.5)) +
    labs(x = "Akula 2021 log2(fold change)", y = "Current analysis log2(fold change)", 
         title = "Correlation between DGE analysis results") +
    theme(legend.position = c(0.05, 0.7),
          legend.spacing.x = unit(0.1, "mm"),
          legend.background = element_blank(),
          legend.box.background = element_rect(colour = "black")
    )
map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(figures_dir, "correlate_with_Akula_DEGs", .x), width = 4, height = 3)
)



# Functional enrichments --------------------------------------------------

### GENE ONTOLOGY ###

# Run Gene Ontology enrichment on DEGs with pvalue < 0.05 (n = 1373)
df_go_res <- enrichGO(gene = my_degs,
                      OrgDb = "org.Hs.eg.db",
                      universe = ensembl_ids,
                      keyType = "ENSEMBL",
                      ont = "ALL",
                      pAdjustMethod = "fdr",
                      pvalueCutoff = 1,
                      qvalueCutoff = 1,
                      readable = TRUE
) %>% 
    as_tibble()

# Plot
f_plot_go_by_ontology(df_go_res, 
                      pathway_text_width = 45,
                      pathway_text_size = 9,
                      plot_title = "DEG GO pathway results",
                      n_facet_rows = 2,
                      legend_position = c(0.7, 0.25)
)
map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(figures_dir, "DEGs_GO_enrichment", .x), width = 7, height = 5)
)

## Export as table
write_xlsx(df_go_res, path = paste0(figures_dir, "DEG_GO_pathway_enrichments.xlsx"))

### GSEA ###

# Run GSEA for each ontology
my_degs_ranked <- df_de_res %>% 
    arrange(-log2fold_change) %>% 
    dplyr::select(ensembl_gene_id, log2fold_change) %>% 
    deframe
set.seed(20240730)
df_gsea_res <- map_dfr(
    .x = c("BP", "CC", "MF"),
    .f = ~ gseGO(my_degs_ranked,
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

# Plot
f_plot_gsea_by_ontology(df_gsea_res, 
                        pathway_text_size = 2.5,
                        pathway_text_width = 25,
                        pathway_text_overlaps = 10,
                        pathway_box_padding = 0.25,
                        plot_title = "DEGs GSEA enrichments",
                        legend_position = "none"
)
map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(figures_dir, "DEGs_GSEA_enrichment", .x), width = 4.25, height = 5.5)
)

## Export as table
write_xlsx(df_gsea_res, path = paste0(figures_dir, "DEG_GSEA_enrichments.xlsx"))
