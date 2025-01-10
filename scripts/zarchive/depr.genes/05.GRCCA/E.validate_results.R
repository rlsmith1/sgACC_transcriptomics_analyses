
################################################################################

# GENE-LEVEL: Validate GRCCA results using external data

################################################################################

base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
source(paste0(base_dir, "scripts/full_analysis_scripts/genes/GRCCA/04.3.load_GRCCA_results.R"))
figures_dir <- paste0(base_dir, "outputs/benchmarking_res/")

# load data ---------------------------------------------------------------


## GWAS DATA
df_scz_gwas <- readxl::read_excel(paste0(base_dir, "data/trubetskoy_2023_scz_gwas.xlsx"),
                                  sheet = 2
)

## DE data
load(paste0(base_dir, "objects/akula_et_al_DE_results.RDS")) # df_akula_res

## Consensus DE data
df_consensus_degs <- readxl::read_excel(paste0(base_dir, "data/merikangas_2022_consensus_DEGs.xlsx")) %>% 
    left_join(df_ensembl_to_symbol)

## Rare variants
df_exome_urvs <- readxl::read_excel(paste0(base_dir, "data/singh_2022_exome_URVs.xlsx")) %>% 
    left_join(df_ensembl_to_symbol)


# GWAS hypergeometric test -----------------------------------------------------

# pull genes in analysis
gene_universe <- colnames(X.mat)
length(gene_universe) # 18677

# pull GWAS hits
gwas_genes <- df_scz_gwas %>% 
    pull(gene_ensembl) %>% 
    unique
length(gwas_genes) # 656

# find GWAS hits that were part of our analysis
gwas_genes_in_universe <- intersect(gene_universe, gwas_genes) %>% unique
length(gwas_genes_in_universe) # 442

# pull significant GRCCA genes
grcca_genes <- df_x_res %>%
    filter(significant == 1) %>%
    pull(ensembl_gene_id) %>%
    unique()
length(grcca_genes) # 2026

# find genes that are both GWAS hits and GRCCA hits
overlap_genes <- intersect(grcca_genes, gwas_genes_in_universe)
length(overlap_genes) # 56

phyper_gwas <- 1 - phyper(q = length(overlap_genes), 
                          m = length(gwas_genes_in_universe), 
                          n = length(gene_universe) - length(gwas_genes_in_universe), 
                          k = length(grcca_genes)
) 
phyper_gwas # p = 0.095


# GWAS enrichment test ---------------------------------------------------------

ranked_gene_list <- df_x_res %>% 
    filter(significant == 1) %>% 
    arrange(-pearsons_r) %>% 
    dplyr::select(ensembl_gene_id, pearsons_r) %>% 
    distinct() %>% 
    deframe
m_gwas <- df_x_res %>% 
    dplyr::select(ensembl_gene_id) %>% 
    mutate(gwas = ifelse(ensembl_gene_id %in% gwas_genes, "hit", "no_hit"), .before = 1)

set.seed(20240306)
sig_gwas_res <- GSEA(ranked_gene_list, TERM2GENE = m_gwas, pvalueCutoff = 1) %>% as_tibble() %>% clean_names
p_enrich <- sig_gwas_res %>% pull(pvalue)
nes_enrich <- sig_gwas_res %>% pull(nes)


# Consensus DEG overlap --------------------------------------------

df_consensus_degs %>% 
    count(direction_of_effect)

## Hypergeometric

# pull genes in analysis
gene_universe <- df_x_res %>% pull(gene_symbol) %>% unique
length(gene_universe) # 18640

# pull DEG hits
consensus_degs <- df_consensus_degs %>% 
    pull(gene_symbol) %>% 
    unique
length(consensus_degs) # 160

# find GWAS hits that were part of our analysis
consensus_degs_in_universe <- intersect(gene_universe, consensus_degs) %>% unique
length(consensus_degs_in_universe) # 141

# pull significant GRCCA genes
grcca_genes <- df_x_res %>%
    filter(significant == 1) %>%
    pull(gene_symbol) %>%
    unique()
length(grcca_genes) # 2025

# find genes that are both GWAS hits and GRCCA hits
overlap_genes <- intersect(grcca_genes, consensus_degs_in_universe)
length(overlap_genes) # 13

phyper_consensus_degs <- 1 - phyper(q = length(overlap_genes), 
                          m = length(consensus_degs_in_universe), 
                          n = length(gene_universe) - length(consensus_degs_in_universe), 
                          k = length(grcca_genes)
) 
phyper_consensus_degs # p = 0.68

## GSEA
ranked_gene_list <- df_x_res %>% 
    filter(significant == 1) %>% 
    arrange(-pearsons_r) %>% 
    dplyr::select(gene_symbol, pearsons_r) %>% 
    distinct() %>% 
    deframe
m_consensus_degs <- df_x_res %>% 
    dplyr::select(gene_symbol) %>% 
    mutate(deg = ifelse(gene_symbol %in% consensus_degs, "hit", "no_hit"), .before = 1)

set.seed(20240710)
sig_consensus_degs_res <- GSEA(ranked_gene_list, TERM2GENE = m_consensus_degs, pvalueCutoff = 1) %>% as_tibble() %>% clean_names
sig_consensus_degs_res

## Plot position of consensus DEGs in distribution
df_x_res_degs <- df_x_res %>% 
    filter(gene_symbol %in% consensus_degs)
df_x_res %>% 
    ggplot(aes(x = z_score)) +
    geom_density() +
    geom_point(data = df_x_res_degs, shape = 21, size = 3,
               aes(x = z_score, y = 0, fill = z_score)) +
    geom_vline(xintercept = c(-2, 2), lty = 2) +
    geom_text_repel(data = df_x_res_degs, min.segment.length = 0.1, box.padding = 1, size = 3,
                    aes(x = z_score, y = 0, label = gene_symbol)) +
    scale_fill_gradientn(colors = rev(brewer.pal(n = 11, "RdBu")), guide = "none") +
    labs(x = "Gene weight z-score", y = "", title = "Consensus DEG position in GRCCA results")
map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(figures_dir, "GRCCA_genes_consensus_DEG_dist", .x), width = 4, height = 3)
)


# Rare variant overlap ---------------------------------------------

## Hypergeometric

# pull genes in analysis
gene_universe <- df_x_res %>% pull(ensembl_gene_id) %>% unique
length(gene_universe) # 18677

# pull DEG hits
exome_urvs <- df_exome_urvs %>% 
    pull(ensembl_gene_id) %>% 
    unique
length(exome_urvs) # 10

# find GWAS hits that were part of our analysis
exome_urvs_in_universe <- intersect(gene_universe, exome_urvs) %>% unique
length(exome_urvs_in_universe) # 9

# pull significant GRCCA genes
grcca_genes <- df_x_res %>%
    filter(significant == 1) %>%
    pull(ensembl_gene_id) %>%
    unique()
length(grcca_genes) # 2026

# find genes that are both GWAS hits and GRCCA hits
overlap_genes <- intersect(grcca_genes, exome_urvs_in_universe)
length(overlap_genes) # 2

phyper_exome_urvs <- 1 - phyper(q = length(overlap_genes), 
                                    m = length(exome_urvs_in_universe), 
                                    n = length(gene_universe) - length(exome_urvs_in_universe), 
                                    k = length(grcca_genes)
) 
phyper_exome_urvs # p = 0.68

## GSEA
ranked_gene_list <- df_x_res %>% 
    filter(significant == 1) %>% 
    arrange(-pearsons_r) %>% 
    dplyr::select(ensembl_gene_id, pearsons_r) %>% 
    distinct() %>% 
    deframe
m_exome_urvs <- df_x_res %>% 
    dplyr::select(ensembl_gene_id) %>% 
    mutate(urv = ifelse(ensembl_gene_id %in% exome_urvs, "hit", "no_hit"), .before = 1)

set.seed(20240710)
sig_exome_urvs_res <- GSEA(ranked_gene_list, TERM2GENE = m_exome_urvs, pvalueCutoff = 1) %>% as_tibble() %>% clean_names
sig_exome_urvs_res

## Plot position of URV genes in distribution
df_x_res_urvs <- df_x_res %>% 
    filter(ensembl_gene_id %in% exome_urvs)
df_x_res %>% 
    ggplot(aes(x = z_score)) +
    geom_density() +
    geom_point(data = df_x_res_urvs, shape = 21, size = 3,
               aes(x = z_score, y = 0, fill = z_score)) +
    geom_vline(xintercept = c(-2, 2), lty = 2) +
    geom_text_repel(data = df_x_res_urvs, min.segment.length = 0.1, box.padding = 1, size = 3,
                    aes(x = z_score, y = 0, label = gene_symbol)) +
    scale_fill_gradientn(colors = rev(brewer.pal(n = 11, "RdBu")), guide = "none") +
    labs(x = "Gene weight z-score", y = "", title = "Rare variant position in GRCCA results")
map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(figures_dir, "GRCCA_genes_rare_variant_dist", .x), width = 4, height = 3)
)



# Fig 4A | SCZ GWAS genes x GRCCA hits ------------------------------------


### PLOT GWAS x GRCCA hits

df_x_res %>% 
    filter(ensembl_gene_id %in% overlap_genes) %>% 
    arrange(pearsons_r) %>% 
    mutate(gene_symbol = factor(gene_symbol, levels = unique(.$gene_symbol))) %>% 
    
    ggplot(aes(x = pearsons_r, y = gene_symbol, fill = module)) +
    geom_segment(aes(x = 0, xend = pearsons_r, color = module)) +
    geom_point(aes(size = -log10(p_adj)), shape = 21, color = "black") +
    geom_vline(xintercept = 0, lty = 2) +
    scale_fill_manual(values = module_colors) +
    scale_color_manual(values = module_colors) +
    scale_alpha_continuous(range = c(0.5, 1)) +
    labs(y = "", x = "Structure correlation",
         title = paste0("Structure correlations of GRCCA x GWAS hits (n = ", length(overlap_genes),")"),
         #title = bquote(bold("E")~" | Significant GRCCA genes in SCZ GWAS"),
         caption = paste0("p enrichment = ", round(p_enrich, 3), "\n",
                          "NES = ", round(nes_enrich, 3))
    ) +
    guides(fill = guide_legend(ncol = 2)) +
    theme(legend.position = c(0.22, 0.70),
          legend.direction = "vertical",
          legend.box = "vertical"
    )

# save to project dir
map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(base_dir, "outputs/figures/for_manuscript/genes/grcca/withGray/4A.overlap_SCZ_GWAS_gene_list", .x), # or grcca/withGray
                  width = 6.5, height = 8.5)
)

## ALT PLOT: HEATMAP?
m_overlap_genes_expr <- df_vsd_regress %>% 
    dplyr::select(sample, all_of(overlap_genes)) %>% 
    column_to_rownames("sample") %>% 
    as.matrix %>% 
    scale
sample_order <- m_overlap_genes_expr[hclust(dist(m_overlap_genes_expr))$order,] %>% rownames
gene_order <- m_overlap_genes_expr[hclust(dist(m_overlap_genes_expr))$order,] %>% colnames
library(viridis)
df_vsd_regress %>% 
    dplyr::select(sample, all_of(overlap_genes)) %>% 
    pivot_longer(contains("ENSG"), names_to = "ensembl_gene_id", values_to = "expr") %>% 
    left_join(df_ensembl_to_symbol) %>% 
    mutate(dx = str_remove(sample, ".*_"),
           dx = case_when(
               dx == "control" ~ "Control",
               dx == "bipolar" ~ "BD",
               dx == "mdd" ~ "MDD",
               dx == "schizo" ~ "SCZ"
           ) %>% factor(levels = names(dx_colors))
    ) %>% 
    arrange(dx) %>% 
    mutate(sample = factor(sample, levels = unique(.$sample))) %>% 
    mutate(ensembl_gene_id = factor(ensembl_gene_id, levels = gene_order)) %>% 
    arrange(ensembl_gene_id) %>% 
    mutate(gene_symbol = factor(gene_symbol, levels = unique(.$gene_symbol))) %>% 
    
    ggplot(aes(x = sample, y = gene_symbol)) +
    geom_tile(aes(fill = expr)) +
    scale_fill_viridis() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))


# Fig 4B | GWAS enrichment across significant thresholds ---------------------------------------------


### CALCULATE ENRICHMENT P-VALUE ACROSS THRESHOLDS
z_thresholds <- seq(1.5, 5, by = 0.01)
p_thresholds <- seq(0.05, 0.0001, by = -0.01)
df_gwas_enrichment <- tibble()
for(z in z_thresholds) {
    for(p in p_thresholds) {
        ranked_gene_list <- df_x_res %>% 
            filter(p_adj <= p & abs(z_score) > z) %>% 
            arrange(-pearsons_r) %>% 
            dplyr::select(ensembl_gene_id, pearsons_r) %>% 
            distinct() %>% 
            deframe
        
        m_gwas <- df_x_res %>% 
            dplyr::select(ensembl_gene_id) %>% 
            mutate(gwas = ifelse(ensembl_gene_id %in% gwas_genes, "hit", "no_hit"), .before = 1)
        
        set.seed(20240306)
        sig_gwas_res <- GSEA(ranked_gene_list, TERM2GENE = m_gwas, pvalueCutoff = 1) %>% as_tibble() %>% clean_names
        p_enrichment <- sig_gwas_res %>% pull(p_adjust)
        
        df_tmp <- tibble(
            z_score = z,
            p_adj = p,
            p_enrichment = p_enrichment
        )
        df_gwas_enrichment <- df_gwas_enrichment %>% bind_rows(df_tmp)
        
    }
}

### PLOT

df_gwas_enrichment %>% #filter(z_score <= 3) %>% 
    ggplot(aes(x = z_score, y = -log10(p_enrichment))) +
    geom_point(aes(fill = -log10(p_enrichment)), shape = 21) +
    geom_hline(yintercept = -log10(0.05), lty = 2, color = "red") +
    geom_vline(xintercept = c(1.96, 2.58)) +
    annotate(geom = "text", x = 1.96, y = 0.25, label = "p = 0.05", hjust = -0.1, size = 3.5) +
    annotate(geom = "text", x = 2.58, y = 0, label = "p = 0.01", hjust = -0.1, size = 3.5) +
    #annotate(geom = "text", x = 4.5, y = 1.4, label = "Enrichment p = 0.05", hjust = 0.85, size = 3) +
    facet_wrap(vars(p_adj), ncol = 2) +
    scale_fill_gradient(low = "white", high = "midnightblue") +
    guides(fill = guide_colorbar(title.position = "top", title.hjust = 0, title = "-log10(Enrichment p-value)")) +
    labs(x = "Abs(z-score) threshold", y = "-log10(Enrichment p-value)",
         caption = "Faceted by correlation p-value threshold \nRed line: enrichment p = 0.05 (anything above is significant)",
         title = "GRCCA genes enrichment for SCZ GWAS across significance thresholds") +
    theme(legend.position = c(0.75, 0.15),
          legend.box = "vertical",
          legend.direction = "horizontal")


map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(base_dir, "outputs/figures/for_manuscript/genes/grcca/withGray/4B.enrichment_SCZ_GWAS_thresholds", .x), # or grcca/withGray
                  width = 6.5, height = 7.5)
)




# DEGs hypergeometric test ------------------------------------------------

# pull genes in analysis
gene_universe <- colnames(X.mat)
length(gene_universe) # 18677

# pull DEGs
SCZ_DEGs <- df_akula_res %>% 
    filter(level == "gene" & padj < 0.05 & str_detect(comparison, "schizo")) %>% 
    pull(id) %>% 
    unique
length(SCZ_DEGs) # 95

# find GWAS hits that were part of our analysis
SCZ_DEGs_in_universe <- intersect(gene_universe, SCZ_DEGs) %>% unique
length(SCZ_DEGs_in_universe) # 69

# pull significant GRCCA genes
grcca_genes <- df_x_res %>%
    filter(significant == 1) %>%
    pull(ensembl_gene_id) %>%
    unique()
length(grcca_genes) # 2026

# find genes that are both GWAS hits and GRCCA hits
overlap_genes <- intersect(grcca_genes, SCZ_DEGs_in_universe)
length(overlap_genes) # 11

phyper_DEGs <- 1 - phyper(q = length(overlap_genes), 
                          m = length(SCZ_DEGs_in_universe), 
                          n = length(gene_universe) - length(SCZ_DEGs_in_universe), 
                          k = length(grcca_genes)
) 
phyper_DEGs # p = 0.07



# Fig 4C | Module assignments of GRCCAxDEGs -------------------------------------------------------


df_modules_filt %>% 
    filter(ensembl_gene_id %in% overlap_genes) %>% 
    dplyr::count(module, color) %>% 
    
    ggplot(aes(x = n, y = fct_reorder(module, n), fill = I(color))) +
    geom_col(color = "black") +
    geom_text(aes(label = n), hjust = -1) +
    labs(y = "", x = "Number of genes that are differentially expressed",
         title = "Module assignment of significant GRCCA x DEGs \n(Akula et al 2021)") +
    theme(legend.position = "none")

# save to project dir
map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(base_dir, "outputs/figures/for_manuscript/genes/grcca/withGray/4C.DEGxGRCCA_modules", .x), # or grcca/withGray
                  width = 5, height = 3)
)


# DEG OVERLAP WITH GWAS
DEG_GWAS_overlap_genes <- intersect(SCZ_DEGs, gwas_genes_in_universe)
length(DEG_GWAS_overlap_genes) # 2 _withGray; n =  _WITHGRAY

phyper_gwas <- 1 - phyper(q = length(DEG_GWAS_overlap_genes), 
                          m = length(gwas_genes_in_universe), 
                          n = length(gene_universe) - length(gwas_genes_in_universe), 
                          k = length(SCZ_DEGs)
) 
phyper_gwas # p = 0.39 _withGray; p = 0.015 _WITHGRAY


# DEG ENRICHMENT
sig_deg_list <- df_akula_res %>% 
    filter(padj <= 0.05 & str_detect(comparison, "schizo") & level == "gene") %>% 
    arrange(-log2fold_change) %>% 
    dplyr::select(id, log2fold_change) %>% 
    distinct() %>% 
    deframe
m_gwas <- df_x_res %>% 
    dplyr::select(ensembl_gene_id) %>% 
    mutate(gwas = ifelse(ensembl_gene_id %in% gwas_genes, "hit", "no_hit"), .before = 1)

set.seed(20240306)
sig_gwas_res <- GSEA(sig_deg_list, TERM2GENE = m_gwas, pvalueCutoff = 1) %>% as_tibble() %>% clean_names
sig_gwas_res




# Fig 4D | Overrepresented modules in SCZ DEG list ----------------------------------------------------------------


# DEG HYPERGEOMETRIC

# pull DEGs
SCZ_DEGs <- df_akula_res %>% 
    filter(level == "gene" & padj < 0.05 & str_detect(comparison, "schizo")) %>% 
    pull(id) %>% 
    unique
length(SCZ_DEGs) # 95

# run hypergeometric test for each module
df_hypergeometric_degs <- df_modules_filt %>% 
    filter(ensembl_gene_id %in% SCZ_DEGs) %>% 
    dplyr::count(module, color) %>% 
    dplyr::rename("q" = "n") %>% 
    left_join(
        df_modules_filt %>% 
            dplyr::count(module, color) %>% 
            dplyr::rename("k" = "n")
    ) %>% 
    mutate(m = length(SCZ_DEGs),
           n = length(gene_universe) - length(SCZ_DEGs),
           p_value = phyper(q, m, n, k, lower.tail = FALSE),
           p_adj = p.adjust(p_value, method = "fdr")
    ) %>% 
    arrange(p_value)

# PLOT
df_hypergeometric_degs %>%
    bind_rows(df_modules_filt %>% 
                  dplyr::count(module) %>% 
                  filter(!(module %in% c(df_hypergeometric_degs %>% pull(module)))) %>% 
                  dplyr::rename("k" = "n")
    ) %>% 
    mutate(q = ifelse(is.na(q), 0, q),
           p_adj = ifelse(is.na(p_adj), 1, p_adj)
    ) %>%
    
    ggplot() +
    geom_segment(aes(x = module, xend = module, y = 0, yend = -log10(p_adj))) +
    #geom_col(width = 0.05, fill = "black") +
    geom_point(aes(x = module, y = -log10(p_adj), fill = I(color), size = (q/k)*100), shape = 21) +
    geom_hline(aes(yintercept = -log10(0.01)), linewidth = 0.25) +
    annotate(geom = "text", label = "FDR < 0.01", x = 23, y = -log10(0.01) + 0.25) +
    guides(size = guide_legend(title = "Percent of module that was \n differentially expressed")) +
    labs(x = "", y = "-log10(hypergeometric FDR)",
         title = "Modules with overrepresentation of DEGs") +
    theme(legend.position = c(0.85, 0.75),
          axis.text.x = element_text(angle = 45, hjust = 0.9, size = 7)
          )

# save to project dir
map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(base_dir, "outputs/figures/for_manuscript/genes/grcca/withGray/4D.DEG_modules_hypergeometric", .x), # or grcca/withGray
                  width = 6.5, height = 5)
)


