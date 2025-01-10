
################################################################################

# GENE-LEVEL: Compare RCCA vs GRCCA results in the context of the SCZ literature

################################################################################

base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
source(paste0(base_dir, "scripts/full_analysis_scripts/genes/load_WGCNA_res.R"))
figures_dir <- paste0(base_dir, "outputs/GRCCA_res/")
tables_dir <- paste0(base_dir, "outputs/tables/02Aug_update/")


# Load RCCA & GRCCA data ---------------------------------------------------------

df_rcca_x_res <- read_xlsx(paste0(tables_dir, "gene_RCCA_results.xlsx"), sheet = 3)
df_grcca_x_res <- read_xlsx(paste0(tables_dir, "gene_GRCCA_results.xlsx"), sheet = 3)

## Consensus DE data
df_consensus_degs <- readxl::read_excel(paste0(base_dir, "data/merikangas_2022_consensus_DEGs.xlsx")) %>% 
    left_join(df_ensembl_to_symbol)



# Combine RCCA and GRCCA z-score and structure correlations ---------------

df_rcca_grcca_res <- df_grcca_x_res %>% 
    dplyr::select(ensembl_gene_id, gene_symbol, z_score, pearsons_r, p_adj) %>% 
    dplyr::rename_at(3:5, ~ paste0("grcca_", .)) %>% 
    left_join(
        df_rcca_x_res %>% 
            dplyr::select(ensembl_gene_id, gene_symbol, z_score, pearsons_r, p_adj) %>% 
            dplyr::rename_at(3:5, ~ paste0("rcca_", .)),
        by = join_by(ensembl_gene_id, gene_symbol)
        
    )


# Plot correlation between gene weight vectors ----------------------------

## Correlate z-scores
df_rcca_grcca_res %>% 
    ggplot(aes(x = rcca_z_score, y = grcca_z_score)) +
    geom_point(aes(color = rcca_z_score), size = 0.75) +
    geom_abline(lty = 2, color = "black") +
    geom_smooth(method = "lm", color = "black", se = FALSE) +
    stat_regline_equation(label.y.npc = "top", size = 3) +
    stat_cor(label.y.npc = "top", size = 3, vjust = 3) +
    scale_color_gradientn(colors = brewer.rdbu(100) %>% rev, guide = "none") +
    labs(x = "RCCA z-score", y = "GRCCA z-score",
         title = "Relationship between RCCA and GRCCA \ngene weight z-scores")
map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(figures_dir, "RCCA_GRCCA_zscore_cor", .x), width = 3.5, height = 3)
)

## Correlate structure correlations
df_rcca_grcca_res %>% 
    ggplot(aes(x = rcca_pearsons_r, y = grcca_pearsons_r)) +
    geom_point(aes(color = rcca_pearsons_r), size = 0.75) +
    geom_abline(lty = 2, color = "black") +
    geom_smooth(method = "lm", color = "black", se = FALSE) +
    stat_regline_equation(label.y.npc = "top", size = 3) +
    stat_cor(label.y.npc = "top", size = 3, vjust = 3) +
    scale_color_gradientn(colors = brewer.rdbu(100) %>% rev, guide = "none") +
    labs(x = "RCCA structure correlation", y = "GRCCA structure correlation",
         title = "Relationship between RCCA and GRCCA \ngene weight structure correlations")
map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(figures_dir, "RCCA_GRCCA_struct_cor", .x), width = 3.5, height = 3)
)



# Overlap of significant RCCA and GRCCA genes (venn diagram) --------------

sig_rcca_genes <- df_rcca_x_res %>% dplyr::filter(significant == 1) %>% pull(ensembl_gene_id)
sig_grcca_genes <- df_grcca_x_res %>% dplyr::filter(significant == 1) %>% pull(ensembl_gene_id)

sig_overlap <- intersect(sig_rcca_genes, sig_grcca_genes)
length(sig_overlap)
length(sig_rcca_genes) - length(sig_overlap)
length(sig_grcca_genes) - length(sig_overlap)


# Functional enrichment of RCCA and GRCCA-specific genes ------------------

### RCCA sig genes ###

## Run GO
df_rcca_go <- enrichGO(gene = sig_rcca_genes[!(sig_rcca_genes %in% sig_overlap)],
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

## Plot
f_plot_go_by_ontology(df_rcca_go, 
                      pathway_text_width = 35,
                      pathway_text_size = 10,
                      plot_title = "RCCA-only genes enrichment",
                      n_facet_rows = 1,
                      legend_position = "none"
)

## Save
map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(figures_dir, "RCCA_only_GO_enrichment", .x), # or grcca/withGray
                  width = 10, height = 3.5)
)

### GRCCA sig genes ###

## Run GO
df_grcca_go <- enrichGO(gene = sig_grcca_genes[!(sig_grcca_genes %in% sig_overlap)],
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

## Plot
f_plot_go_by_ontology(df_grcca_go, 
                      pathway_text_width = 35,
                      pathway_text_size = 10,
                      plot_title = "GRCCA-only genes enrichment",
                      n_facet_rows = 1,
                      legend_position = "none"
)

## Save
map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(figures_dir, "GRCCA_only_GO_enrichment", .x), # or grcca/withGray
                  width = 10, height = 3.5)
)


# Load benchmarking data ---------------------------------------------------------------

## Current analysis DEGs
load(paste0(base_dir, "objects/DE_results.RDS")) # df_de_res
my_degs <- df_de_res %>% dplyr::filter(pvalue < 0.05) %>% pull(ensembl_gene_id)

## Akula DEGs
load(paste0(base_dir, "objects/akula_et_al_DE_results.RDS")) # df_akula_res
akula_degs <- df_akula_res %>% dplyr::filter(level == "gene" & comparison == "schizo_ctrl" & pvalue < 0.05) %>% pull(id)

## Consensus DE data
df_consensus_degs <- readxl::read_excel(paste0(base_dir, "data/merikangas_2022_consensus_DEGs.xlsx")) %>% 
    left_join(df_ensembl_to_symbol)
consensus_degs <- df_consensus_degs %>% dplyr::filter(!is.na(ensembl_gene_id)) %>% pull(ensembl_gene_id) %>% unique

## GWAS DATA
df_scz_gwas <- readxl::read_excel(paste0(base_dir, "data/trubetskoy_2023_scz_gwas.xlsx"),
                                  sheet = 3
)
common_variants <- df_scz_gwas %>% dplyr::filter(str_detect(gene_ensembl, "^ENSG0")) %>% pull(gene_ensembl) %>% unique

# sheet 2 = 95% credible set
# sheet 3 = 95% credible set k <= 3.5
# sheet 4 = 95% credible set <= 5 SNPs
# sheet 5 = 95% credible set <= 5 PIP

## Rare variants
df_exome_urvs <- readxl::read_excel(paste0(base_dir, "data/singh_2022_exome_URVs.xlsx"))
rare_variants <- df_exome_urvs %>% pull(ensembl_gene_id) %>% unique

### Combine into list ###
l_benchmark_lists <- list(
    "my_degs" = my_degs,
    "akula_degs" = akula_degs,
    "consensus_degs" = consensus_degs,
    "common_variants" = common_variants,
    "rare_variants" = rare_variants
)


# test hypergeometric overlap of GRCCA and RCCA hits across thresholds with each benchmark list ------------------


### Loop through and run hypergeometric tests across thresholds ###

thresholds <- seq(0.01, 0.25, by = 0.01)
df_rcca_grcca_hypergeometric_benchmark <- tibble()
for (list in 1:length(l_benchmark_lists)) {
    
    print(paste0("Calculating overlap with ", names(l_benchmark_lists[list])))
    
    for (thresh in thresholds) {
        
        n_genes <- thresh*length(ensembl_ids)
        
        # RCCA hypergeometric
        df_rcca_hypergeometric_res <- f_hypergeometric(gene_list1 = df_rcca_grcca_res %>% 
                                                           top_n(n = n_genes, wt = abs(rcca_pearsons_r)) %>% 
                                                           #top_n(n = n_genes, wt = abs(rcca_z_score)) %>% 
                                                           pull(ensembl_gene_id), 
                                                       gene_list2 = l_benchmark_lists[[list]],
                                                       gene_universe = ensembl_ids
        ) %>% as_tibble() %>% 
            count(p_value) %>% 
            mutate(model = "RCCA", .before = 1)
        
        # GRCCA hypergeometric
        df_grcca_hypergeometric_res <- f_hypergeometric(gene_list1 = df_rcca_grcca_res %>% 
                                                            top_n(n = n_genes, wt = abs(grcca_pearsons_r)) %>% 
                                                            #top_n(n = n_genes, wt = abs(grcca_z_score)) %>% 
                                                            pull(ensembl_gene_id), 
                                                        gene_list2 = l_benchmark_lists[[list]],
                                                        gene_universe = ensembl_ids
        ) %>% as_tibble() %>% 
            count(p_value) %>% 
            mutate(model = "GRCCA", .before = 1)
        
        # Combine results
        df_tmp <- bind_rows(df_rcca_hypergeometric_res, df_grcca_hypergeometric_res) %>% 
            dplyr::rename("overlap_n" = "n") %>% 
            mutate(
                benchmark_list = names(l_benchmark_lists[list]),
                threshold = thresh, .before = 1
            ) %>% 
            mutate(total_n = length(l_benchmark_lists[[list]]))
        df_rcca_grcca_hypergeometric_benchmark <- df_rcca_grcca_hypergeometric_benchmark %>% bind_rows(df_tmp)
        
    }
    
}

df_rcca_grcca_hypergeometric_benchmark <- df_rcca_grcca_hypergeometric_benchmark %>% 
    mutate(p_adj = p.adjust(p_value, method = "fdr"))

## Plot results
df_rcca_grcca_hypergeometric_benchmark %>%
    ggplot(aes(x = threshold, y = -log10(p_adj))) +
    geom_point(aes(color = model)) +
    geom_line(aes(group = model, color = model)) +
    geom_hline(aes(yintercept = -log10(0.05)), lty = 2) +
    facet_wrap(vars(benchmark_list), nrow = 1, scales = "free_y") +
    labs(title = "RCCA vs GRCCA hypergeometric overlaps across structure correlation thresholds") +
    theme(legend.position = "top",
          legend.justification = "right")
map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(figures_dir, "RCCA_GRCCA_struct_cor_hypergeometric_thresholds", .x), width = 9, height = 3)
)



# Hypergeometric test of benchmark lists in RCCA/GRCCA-specific gene lists --------------------

df_rcca_grcca_only_hypergeometric_benchmark <- tibble()
for (list in 1:length(l_benchmark_lists)) {
    
    # RCCA hypergeometric
    df_rcca_only_hypergeometric_res <- f_hypergeometric(gene_list1 = sig_rcca_genes[!(sig_rcca_genes %in% sig_overlap)], 
                                                        gene_list2 = l_benchmark_lists[[list]],
                                                        gene_universe = ensembl_ids
    ) %>% as_tibble() %>% 
        count(p_value) %>% 
        mutate(model = "RCCA", .before = 1)
    
    # GRCCA hypergeometric
    df_grcca_only_hypergeometric_res <- f_hypergeometric(gene_list1 = sig_grcca_genes[!(sig_grcca_genes %in% sig_overlap)], 
                                                         gene_list2 = l_benchmark_lists[[list]],
                                                         gene_universe = ensembl_ids
    ) %>% as_tibble() %>% 
        count(p_value) %>% 
        mutate(model = "GRCCA", .before = 1)
    
    # Combine results
    df_tmp <- bind_rows(df_rcca_only_hypergeometric_res, df_grcca_only_hypergeometric_res) %>% 
        dplyr::rename("overlap_n" = "n") %>% 
        mutate(
            benchmark_list = names(l_benchmark_lists[list]),
        ) %>% 
        mutate(total_n = length(l_benchmark_lists[[list]]))
    df_rcca_grcca_only_hypergeometric_benchmark <- df_rcca_grcca_only_hypergeometric_benchmark %>% bind_rows(df_tmp)
    
}

df_rcca_grcca_only_hypergeometric_benchmark # relatively even


# Were any genes pp-regulated in one analysis & down-regulated in the other? ---------------


df_rcca_grcca_res %>% 
    #dplyr::filter(grcca_p_adj < 0.05 & rcca_p_adj < 0.05) %>% 
    dplyr::filter(
        (grcca_pearsons_r < 0 & rcca_pearsons_r > 0) |
            (rcca_pearsons_r > 0 & grcca_pearsons_r < 0)
    )

df_rcca_grcca_res %>% 
    dplyr::filter(grcca_p_adj < 0.05 & rcca_p_adj < 0.05) %>% 
    dplyr::filter(
        (grcca_z_score < 0 & rcca_z_score > 0) |
            (grcca_z_score > 0 & rcca_z_score < 0)
    ) %>% print(n = nrow(.))



# Calculate orthogonal residuals -----------------------------------------

## Write function to calculate the shortest distance from a point to a line
# m = slope (1 for identity line)
# b = intercept (0 for identity line)
# x0 = x coordinate of point
# y0 = y coordinate of point
f_distance_from_identity <- function(m = 1, b = 0, x0, y0) { abs(m * x0 - y0 + b) / sqrt(m^2 + 1) }

## Calculate distance between each structure correlation point and the identity line
df_rcca_grcca_distance <- df_rcca_grcca_res %>% 
    mutate(distance = f_distance_from_identity(m = 1, b = 0, x = rcca_pearsons_r, y = grcca_pearsons_r)) %>% 
    arrange(-distance)

df_rcca_grcca_distance %>% 
    filter(grcca_p_adj < 0.05 | rcca_p_adj < 0.05 )
df_consensus_degs %>% filter(gene_symbol == "CX3CR1")


