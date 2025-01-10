
###################################################################

# GENE LEVEL: Describe Akula et al 2021 functional enrichment and #
# overlap with published SCZ gene lists.                          #

###################################################################

## Set directories & load data
base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
figures_dir <- paste0(base_dir, "outputs/DGE_res/")
source(paste0(base_dir, "scripts/full_analysis_scripts/genes/setup.R"))


# Load data ---------------------------------------------------------------

## Current analysis DEG results
load(paste0(base_dir, "objects/DE_results.RDS")) # df_de_res

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


# Identify SCZ-control DEGs -------------------------------

## DEGs
my_degs <- df_de_res %>% 
    dplyr::filter(pvalue < 0.05) %>% 
    pull(ensembl_gene_id)

### DEGs (ranked) ###
my_degs_ranked <- df_de_res %>% 
    arrange(-log2fold_change) %>% 
    dplyr::select(ensembl_gene_id, log2fold_change) %>% 
    deframe


# Consensus DEGs -----------------------------------------------------------

### HYPERGEOMETRIC (overall enrichment) ###

# Pull consensus DEG symbols
consensus_degs <- df_consensus_degs %>% pull(ensembl_gene_id)
length(consensus_degs) # n = 161
consensus_degs %in% ensembl_ids %>% sum # n = 154 in our universe

# Run hypergeometric test for gene symbols
consensus_degs_hypergeometric <- f_hypergeometric(gene_list1 = my_degs, 
                                                  gene_list2 = consensus_degs,
                                                  gene_universe = ensembl_ids
) 
n_overlap <- consensus_degs_hypergeometric$overlap_genes %>% length
phyper <- consensus_degs_hypergeometric$p_value
# p = 0.37 for DEGs pvalue < 0.05


### GSEA (pole enrichment) ###

# Create success vs failure matrix for DEGs
m_gsea_consensus_degs <- tibble(ensembl_gene_id = ensembl_ids) %>% 
    mutate(hit = ifelse(ensembl_gene_id %in% consensus_degs, "success", "failure"), .before = 1)

# Run GSEA
set.seed(20240730)
consensus_degs_gsea <- GSEA(my_degs_ranked, TERM2GENE = m_gsea_consensus_degs, pvalueCutoff = 1) %>% as_tibble() %>% clean_names
consensus_degs_gsea
p_enrich <- consensus_degs_gsea %>% pull(pvalue)
nes_enrich <- consensus_degs_gsea %>% pull(nes)

# Plot position of consensus DEGs in L2FC distribution
df_de_res_consensus_degs <- df_de_res %>% 
    dplyr::filter(ensembl_gene_id %in% consensus_degs)
df_de_res %>% 
    ggplot(aes(x = log2fold_change)) +
    geom_density() +
    geom_point(data = df_de_res_consensus_degs, shape = 21, size = 2.5,
               aes(x = log2fold_change, y = 0, fill = log2fold_change, color = padj < 0.05)) +
    geom_vline(xintercept = 0, lty = 2, color = "darkgrey") +
    #geom_vline(xintercept = c(-2, 2), lty = 2) +
    geom_text_repel(data = df_de_res_consensus_degs, 
                    min.segment.length = 0.1, box.padding = 0.5, size = 3, max.overlaps = 50,
                    aes(x = log2fold_change, y = 0, label = gene_symbol)) +
    scale_fill_gradientn(colors = gene_weight_color_scale, guide = "none") +
    scale_color_manual(values = c("transparent", "yellow"), guide = "none") +
    labs(x = "Gene L2FC", y = "", 
         caption = paste0(
           "phyper = ", round(phyper, 2), "; n sig overlap = ", n_overlap,
           "\np GSEA = ", round(p_enrich, 2), "; NES = ", round(nes_enrich, 2)
         ),
         title = "SCZ consensus DEG position \nin DE results")
map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(figures_dir, "DEGs_consensus_DEG_dist", .x), width = 3.5, height = 3)
)


# Common variants ---------------------------------------------------------

### HYPERGEOMETRIC (overall enrichment) ###

# Pull common variant IDs
common_variants <- df_scz_gwas %>% pull(gene_ensembl) %>% unique
length(common_variants) # n = 629
common_variants %in% ensembl_ids %>% sum # n = 426 in our universe

# Run hypergeometric test for gene symbols
common_variant_hypergeometric <- f_hypergeometric(gene_list1 = my_degs, 
                                                  gene_list2 = common_variants,
                                                  gene_universe = ensembl_ids
) 
n_overlap <- common_variant_hypergeometric$overlap_genes %>% length
phyper <- common_variant_hypergeometric$p_value


### GSEA (pole enrichment) ###

# Create success vs failure matrix for DEGs
m_gsea_common_variants <- tibble(ensembl_gene_id = ensembl_ids) %>% 
    mutate(hit = ifelse(ensembl_gene_id %in% common_variants, "success", "failure"), .before = 1)

# Run GSEA
set.seed(20240730)
common_variants_gsea <- GSEA(my_degs_ranked, TERM2GENE = m_gsea_common_variants, pvalueCutoff = 1) %>% as_tibble() %>% clean_names
common_variants_gsea
p_enrich <- common_variants_gsea %>% pull(pvalue)
nes_enrich <- common_variants_gsea %>% pull(nes)

# Plot position of consensus DEGs in L2FC distribution
df_de_res_common_variants <- df_de_res %>% 
    dplyr::filter(ensembl_gene_id %in% common_variants)
df_de_res %>%
    ggplot(aes(x = log2fold_change)) +
    geom_density() +
    geom_point(data = df_de_res_common_variants %>% arrange(-padj), 
               shape = 21, size = 2.5,
               aes(x = log2fold_change, y = 0, fill = log2fold_change, color = padj < 0.05)) +
    geom_vline(xintercept = 0, lty = 2, color = "darkgrey") +
    #geom_vline(xintercept = c(-2, 2), lty = 2) +
    geom_text_repel(data = df_de_res_common_variants %>% dplyr::filter(padj < 0.05),
                    min.segment.length = 0.1, box.padding = 0.5, size = 3, max.overlaps = 50,
                    aes(x = log2fold_change, y = 0, label = gene_symbol)) +
    scale_fill_gradientn(colors = gene_weight_color_scale, guide = "none") +
    scale_color_manual(values = c("transparent", "yellow"), guide = "none") +
    labs(x = "Gene L2FC", y = "", caption = paste0(
        "phyper = ", round(phyper, 2), "; n sig overlap = ", n_overlap,
        "\np GSEA = ", round(p_enrich, 2), "; NES = ", round(nes_enrich, 2)
    ),
    title = "SCZ common variant position in \nDE results")
map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(figures_dir, "DEGs_common_variant_dist", .x), width = 3.5, height = 3)
)


# Rare variants ---------------------------------------------------------

### HYPERGEOMETRIC (overall enrichment) ###

# Pull rare variant IDs
rare_variants <- df_exome_urvs %>% pull(ensembl_gene_id) %>% unique
length(rare_variants) # n = 10
rare_variants %in% ensembl_ids %>% sum # n = 9 in our universe

# Run hypergeometric test for gene symbols
rare_variant_hypergeometric <- f_hypergeometric(gene_list1 = my_degs, 
                                                  gene_list2 = rare_variants,
                                                  gene_universe = ensembl_ids
) 
n_overlap <- rare_variant_hypergeometric$overlap_genes %>% length
phyper <- rare_variant_hypergeometric$p_value


### GSEA (pole enrichment) ###

# Create success vs failure matrix for DEGs
m_gsea_rare_variants <- tibble(ensembl_gene_id = ensembl_ids) %>% 
    mutate(hit = ifelse(ensembl_gene_id %in% rare_variants, "success", "failure"), .before = 1)

# Run GSEA
set.seed(20240730)
rare_variants_gsea <- GSEA(my_degs_ranked, TERM2GENE = m_gsea_rare_variants, pvalueCutoff = 1) %>% as_tibble() %>% clean_names
rare_variants_gsea
p_enrich <- rare_variants_gsea %>% pull(pvalue)
nes_enrich <- rare_variants_gsea %>% pull(nes)

# Plot position of consensus DEGs in L2FC distribution
df_de_res_rare_variants <- df_de_res %>% 
    dplyr::filter(ensembl_gene_id %in% rare_variants)
df_de_res %>%
    ggplot(aes(x = log2fold_change)) +
    geom_density() +
    geom_point(data = df_de_res_rare_variants %>% arrange(-padj), 
               shape = 21, size = 2.5,
               aes(x = log2fold_change, y = 0, fill = log2fold_change, color = padj < 0.05)) +
    geom_vline(xintercept = 0, lty = 2, color = "darkgrey") +
    #geom_vline(xintercept = c(-2, 2), lty = 2) +
    geom_text_repel(data = df_de_res_rare_variants,
                    min.segment.length = 0.1, box.padding = 0.5, size = 3, max.overlaps = 20,
                    aes(x = log2fold_change, y = 0, label = gene_symbol)) +
    scale_fill_gradientn(colors = gene_weight_color_scale, guide = "none") +
    scale_color_manual(values = c("transparent", "yellow"), guide = "none") +
    labs(x = "Gene L2FC", y = "", caption = paste0(
        "phyper = ", round(phyper, 2), "; n sig overlap = ", n_overlap,
        "\np GSEA = ", round(p_enrich, 2), "; NES = ", round(nes_enrich, 2)
    ),
    title = "SCZ rare variant position in \nDE results")
map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(figures_dir, "DEGs_rare_variant_dist", .x), width = 3.5, height = 3)
)

