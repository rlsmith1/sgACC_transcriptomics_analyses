
################################################################################

# GENE-LEVEL: Validate GRCCA results using external data

################################################################################

base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
source(paste0(base_dir, "scripts/full_analysis_scripts/genes/05.GRCCA/C.load_GRCCA_results.R"))
figures_dir <- paste0(base_dir, "outputs/GRCCA_res/")

## Source plot data
source(paste0(base_dir, "scripts/full_analysis_scripts/functions/plot_functions.R"))


# Identify GRCCA significant genes -------------------------------

## Significant genes
sig_grcca_genes <- df_x_res %>% 
    dplyr::filter(abs(z_score) > 2 & p_adj < 0.05) %>% 
    pull(ensembl_gene_id)

## Ranked gene vector
ranked_grcca_genes <- df_x_res %>% 
    arrange(-pearsons_r) %>% 
    dplyr::select(ensembl_gene_id, pearsons_r) %>% 
    deframe


# Load benchmarking data ---------------------------------------------------------------

## Current analysis DEGs
load(paste0(base_dir, "objects/DE_results.RDS")) # df_de_res

## Akula DEGs
load(paste0(base_dir, "objects/akula_et_al_DE_results.RDS")) # df_akula_res

## Consensus DE data
df_consensus_degs <- readxl::read_excel(paste0(base_dir, "data/merikangas_2022_consensus_DEGs.xlsx")) %>% 
    left_join(df_ensembl_to_symbol)

## GWAS DATA
df_scz_gwas <- readxl::read_excel(paste0(base_dir, "data/trubetskoy_2023_scz_gwas.xlsx"),
                                  sheet = 3
)

# sheet 2 = 95% credible set
# sheet 3 = 95% credible set k <= 3.5
# sheet 4 = 95% credible set <= 5 SNPs
# sheet 5 = 95% credible set <= 5 PIP

## Rare variants
df_exome_urvs <- readxl::read_excel(paste0(base_dir, "data/singh_2022_exome_URVs.xlsx")) %>% 
    left_join(df_ensembl_to_symbol)


# Current analysis DEGs ---------------------------------------------------

### HYPERGEOMETRIC (overall enrichment) ###

# Pull significant DEGs
my_degs <- df_de_res %>% dplyr::filter(pvalue < 0.05) %>% pull(ensembl_gene_id) %>% unique
length(my_degs) # n = 1747
my_degs %in% ensembl_ids %>% sum # n = 1747 in our universe

# Run hypergeometric test for gene symbols
my_degs_hypergeometric <- f_hypergeometric(gene_list1 = sig_grcca_genes, 
                                           gene_list2 = my_degs,
                                           gene_universe = ensembl_ids
) 
n_overlap <- my_degs_hypergeometric$overlap_genes %>% length
phyper <- my_degs_hypergeometric$p_value


### GSEA (pole enrichment) ###

# Create success vs failure matrix for DEGs
m_gsea_my_degs <- tibble(ensembl_gene_id = ensembl_ids) %>% 
    mutate(hit = ifelse(ensembl_gene_id %in% my_degs, "success", "failure"), .before = 1)

# Run GSEA
set.seed(20240801)
my_degs_gsea <- GSEA(ranked_grcca_genes, TERM2GENE = m_gsea_my_degs, pvalueCutoff = 1) %>% as_tibble() %>% clean_names
my_degs_gsea
p_enrich <- my_degs_gsea %>% pull(pvalue)
nes_enrich <- my_degs_gsea %>% pull(nes)

# Plot position of consensus DEGs in grcca gene z-score distribution
df_x_res_my_degs <- df_x_res %>% 
    dplyr::filter(ensembl_gene_id %in% my_degs)
f_plot_gene_position(df = df_x_res, 
                     df_subset = df_x_res_my_degs %>% mutate(significant = ifelse(row_number() <= 10, 1, 0)), 
                     x_axis_variable = "structure correlation",
                     plot_title = "Current analysis DEG position in grcca results",
                     max_text_overlaps = 15
) +
    labs(caption = paste0(
        "phyper = ", round(phyper, 2), "; n sig overlap = ", n_overlap,
        "\np GSEA = ", round(p_enrich, 2), "; NES = ", round(nes_enrich, 2)
    )
    )

map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(figures_dir, "GRCCA_current_DEGs_dist", .x), width = 4, height = 3)
)


# Akula DEGs overlap -------------------------------------------------------

### HYPERGEOMETRIC (overall enrichment) ###

# Pull Akula DEG IDs
akula_degs <- df_akula_res %>% dplyr::filter(level == "gene" & comparison == "schizo_ctrl" & pvalue < 0.05) %>% pull(id)
length(akula_degs) # n = 1373
akula_degs %in% ensembl_ids %>% sum # n = 1263 in our universe

# Run hypergeometric test for gene symbols
akula_degs_hypergeometric <- f_hypergeometric(gene_list1 = sig_grcca_genes, 
                                              gene_list2 = akula_degs,
                                              gene_universe = ensembl_ids
) 
n_overlap <- akula_degs_hypergeometric$overlap_genes %>% length
phyper <- akula_degs_hypergeometric$p_value


### GSEA (pole enrichment) ###

# Create success vs failure matrix for Akula DEGs
m_gsea_akula_degs <- tibble(ensembl_gene_id = ensembl_ids) %>% 
    mutate(hit = ifelse(ensembl_gene_id %in% akula_degs, "success", "failure"), .before = 1)

# Run GSEA
set.seed(20240801)
akula_degs_gsea <- GSEA(ranked_grcca_genes, TERM2GENE = m_gsea_akula_degs, pvalueCutoff = 1) %>% as_tibble() %>% clean_names
akula_degs_gsea
p_enrich <- akula_degs_gsea %>% pull(pvalue)
nes_enrich <- akula_degs_gsea %>% pull(nes)

# Plot position of Akula DEGs in grcca gene z-score distribution
df_x_res_akula_degs <- df_x_res %>% 
    dplyr::filter(ensembl_gene_id %in% akula_degs)
f_plot_gene_position(df = df_x_res, 
                     df_subset = df_x_res_akula_degs %>% mutate(significant = ifelse(row_number() <= 10, 1, 0)),
                     x_axis_variable = "structure correlation",
                     plot_title = "Akula DEG position in grcca results",
                     max_text_overlaps = 20) +
    labs(caption = paste0(
        "phyper = ", round(phyper, 2), "; n sig overlap = ", n_overlap,
        "\np GSEA = ", round(p_enrich, 2), "; NES = ", round(nes_enrich, 2)
    )
    )

map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(figures_dir, "GRCCA_Akula_DEG_dist", .x), width = 4, height = 3)
)


# Consensus DEGs -----------------------------------------------------------

### HYPERGEOMETRIC (overall enrichment) ###

# Pull consensus DEG IDs
consensus_degs <- df_consensus_degs %>% dplyr::filter(!is.na(ensembl_gene_id)) %>% pull(ensembl_gene_id)
length(consensus_degs) # n = 161 (154 in our universe)
consensus_degs %in% ensembl_ids %>% sum # n = 154 in our universe

# Run hypergeometric test for gene symbols
consensus_degs_hypergeometric <- f_hypergeometric(gene_list1 = sig_grcca_genes, 
                                                  gene_list2 = consensus_degs,
                                                  gene_universe = ensembl_ids
) 
n_overlap <- consensus_degs_hypergeometric$overlap_genes %>% length
phyper <- consensus_degs_hypergeometric$p_value


### GSEA (pole enrichment) ###

# Create success vs failure matrix for consensus DEGs
m_gsea_consensus_degs <- tibble(ensembl_gene_id = ensembl_ids) %>% 
    mutate(hit = ifelse(ensembl_gene_id %in% consensus_degs, "success", "failure"), .before = 1)

# Run GSEA
set.seed(20240801)
consensus_degs_gsea <- GSEA(ranked_grcca_genes, TERM2GENE = m_gsea_consensus_degs, pvalueCutoff = 1) %>% as_tibble() %>% clean_names
consensus_degs_gsea
p_enrich <- consensus_degs_gsea %>% pull(pvalue)
nes_enrich <- consensus_degs_gsea %>% pull(nes)

# Plot position of consensus DEGs in grcca gene z-score distribution
df_x_res_consensus_degs <- df_x_res %>% 
    dplyr::filter(ensembl_gene_id %in% consensus_degs)
f_plot_gene_position(df = df_x_res, 
                     df_subset = df_x_res_consensus_degs %>% mutate(significant = ifelse(row_number() <= 10, 1, 0)),
                     x_axis_variable = "structure correlation",
                     plot_title = "Consensus DEG position in GRCCA results",
                     max_text_overlaps = 20) +
    labs(caption = paste0(
        "phyper = ", round(phyper, 2), "; n sig overlap = ", n_overlap,
        "\np GSEA = ", round(p_enrich, 2), "; NES = ", round(nes_enrich, 2)
    )
    )

map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(figures_dir, "GRCCA_consensus_DEG_dist", .x), width = 4, height = 3)
)


# Common variants ---------------------------------------------------------


# sig_grcca_genes <- df_x_res %>% 
#     dplyr::filter(abs(z_score) > 2 & p_adj < 0.05) %>% 
#     head(1194) %>% 
#     pull(ensembl_gene_id)

### HYPERGEOMETRIC (overall enrichment) ###

# Pull common variant IDs
common_variants <- df_scz_gwas %>% dplyr::filter(!is.na(gene_ensembl)) %>% pull(gene_ensembl) %>% unique
length(common_variants) # n = 629
common_variants %in% ensembl_ids %>% sum # n = 426 in our universe

# Run hypergeometric test for gene symbols
common_variants_hypergeometric <- f_hypergeometric(gene_list1 = sig_grcca_genes, 
                                                   gene_list2 = common_variants,
                                                   gene_universe = ensembl_ids
) 
n_overlap <- common_variants_hypergeometric$overlap_genes %>% length
phyper <- common_variants_hypergeometric$p_value


### GSEA (pole enrichment) ###

# Create success vs failure matrix for consensus DEGs
m_gsea_common_variants <- tibble(ensembl_gene_id = ensembl_ids) %>% 
    mutate(hit = ifelse(ensembl_gene_id %in% common_variants, "success", "failure"), .before = 1)

# Run GSEA
set.seed(20240801)
common_variants_gsea <- GSEA(ranked_grcca_genes, TERM2GENE = m_gsea_common_variants, pvalueCutoff = 1) %>% as_tibble() %>% clean_names
common_variants_gsea
p_enrich <- common_variants_gsea %>% pull(pvalue)
nes_enrich <- common_variants_gsea %>% pull(nes)

# Plot position of consensus DEGs in grcca gene z-score distribution
df_x_res_common_variants <- df_x_res %>% 
    dplyr::filter(ensembl_gene_id %in% common_variants)
f_plot_gene_position(df = df_x_res, 
                     df_subset = df_x_res_common_variants %>% mutate(significant = ifelse(row_number() <= 10, 1, 0)),
                     x_axis_variable = "structure correlation",
                     plot_title = "Common variant position in GRCCA results",
                     max_text_overlaps = 27) +
    labs(caption = paste0(
        "phyper = ", round(phyper, 2), "; n sig overlap = ", n_overlap,
        "\np GSEA = ", round(p_enrich, 2), "; NES = ", round(nes_enrich, 2)
    )
    )

map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(figures_dir, "GRCCA_common_variant_dist", .x), width = 4, height = 3)
)


# Rare variants ---------------------------------------------------------


sig_grcca_genes <- df_x_res %>%
    dplyr::filter(abs(z_score) > 2 & p_adj < 0.05) %>%
    head(0.1*length(ensembl_ids) %>% ceiling) %>%
    pull(ensembl_gene_id)


### HYPERGEOMETRIC (overall enrichment) ###

# Pull rare variant IDs
rare_variants <- df_exome_urvs %>% dplyr::filter(!is.na(ensembl_gene_id)) %>% pull(ensembl_gene_id) %>% unique
length(rare_variants) # n = 9
rare_variants %in% ensembl_ids %>% sum # n = 9 in our universe

# Run hypergeometric test for gene symbols
rare_variants_hypergeometric <- f_hypergeometric(gene_list1 = sig_grcca_genes, 
                                                 gene_list2 = rare_variants,
                                                 gene_universe = ensembl_ids
) 
n_overlap <- rare_variants_hypergeometric$overlap_genes %>% length
phyper <- rare_variants_hypergeometric$p_value


### GSEA (pole enrichment) ###

# Create success vs failure matrix for consensus DEGs
m_gsea_rare_variants <- tibble(ensembl_gene_id = ensembl_ids) %>% 
    mutate(hit = ifelse(ensembl_gene_id %in% rare_variants, "success", "failure"), .before = 1)

# Run GSEA
set.seed(20240801)
rare_variants_gsea <- GSEA(ranked_grcca_genes, TERM2GENE = m_gsea_rare_variants, pvalueCutoff = 1) %>% as_tibble() %>% clean_names
rare_variants_gsea
p_enrich <- rare_variants_gsea %>% pull(pvalue)
nes_enrich <- rare_variants_gsea %>% pull(nes)

# Plot position of consensus DEGs in grcca gene z-score distribution
df_x_res_rare_variants <- df_x_res %>% 
    dplyr::filter(ensembl_gene_id %in% rare_variants)
f_plot_gene_position(df = df_x_res, 
                     df_subset = df_x_res_rare_variants %>% mutate(significant = ifelse(row_number() <= 10, 1, 0)), 
                     x_axis_variable = "structure correlation",
                     plot_title = "Rare variant position in GRCCA results",
                     max_text_overlaps = 27) +
    labs(caption = paste0(
        "phyper = ", round(phyper, 2), "; n sig overlap = ", n_overlap,
        "\np GSEA = ", round(p_enrich, 2), "; NES = ", round(nes_enrich, 2)
    )
    )

map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(figures_dir, "GRCCA_rare_variant_dist", .x), width = 4, height = 3)
)
