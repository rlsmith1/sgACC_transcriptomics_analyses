
################################################################################

# TRANSCRIPT-LEVEL: Validate GRCCA results using external data

################################################################################

base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
source(paste0(base_dir, "scripts/full_analysis_scripts/transcripts/overall_subset/GRCCA/04.3.load_GRCCA_results.R"))
figures_dir <- paste0(base_dir, "outputs/benchmarking_res/")


# load data ---------------------------------------------------------------

## GWAS DATA
df_scz_gwas <- readxl::read_excel(paste0(base_dir, "data/trubetskoy_2023_scz_gwas.xlsx"),
                                  sheet = 2
) %>% 
    dplyr::select(gene_symbol, gene_ensembl) %>% 
    distinct() %>% 
    left_join(df_gene_to_transcript)

## DE data
load(paste0(base_dir, "objects/akula_et_al_DE_results.RDS")) # df_akula_res

## Consensus DE data
df_consensus_degs <- readxl::read_excel(paste0(base_dir, "data/merikangas_2022_consensus_DEGs.xlsx")) %>% 
    left_join(df_gene_to_transcript)

## Rare variants
df_exome_urvs <- readxl::read_excel(paste0(base_dir, "data/singh_2022_exome_URVs.xlsx")) %>% 
    left_join(df_gene_to_transcript)


# Akula DET overlap -------------------------------------------------------

## Hypergeometric

# pull genes with transcripts that were in analysis
transcript_universe <- df_x_res %>% 
    filter(!is.na(ensembl_transcript_id)) %>% 
    pull(ensembl_transcript_id) %>% 
    unique
length(transcript_universe) # 18031

# pull Akula hits
akula_transcripts <- df_akula_res %>% 
    filter(level == "transcript" & comparison == "schizo_ctrl" & padj < 0.05) %>% 
    pull(id) %>% 
    unique
length(akula_transcripts) # 470

# find Akula hits that were part of our analysis
akula_transcripts_in_universe <- intersect(transcript_universe, akula_transcripts) %>% unique
length(akula_transcripts_in_universe) # 341

# pull significant transcripts
grcca_transcripts <- df_x_res %>%
    filter(significant == 1 & !is.na(ensembl_transcript_id)) %>%
    pull(ensembl_transcript_id) %>%
    unique()
length(grcca_transcripts) # 257

# find genes that are both GWAS hits and GRCCA hits
overlap_transcripts <- intersect(grcca_transcripts, akula_transcripts_in_universe)
length(overlap_transcripts) # 20

phyper_akula <- 1 - phyper(q = length(overlap_transcripts), 
                          m = length(akula_transcripts_in_universe), 
                          n = length(transcript_universe) - length(akula_transcripts_in_universe), 
                          k = length(grcca_transcripts)
) 
phyper_akula # p = 0



# GWAS overlap -----------------------------------------------------

## Hypergeometric

# pull genes in transcript-level analysis
gene_universe <- df_x_res %>% 
    filter(!is.na(ensembl_gene_id)) %>% 
    pull(ensembl_gene_id) %>% 
    unique
length(gene_universe) # 18031

# pull GWAS hits
gwas_genes <- df_scz_gwas %>% 
    pull(gene_ensembl) %>% 
    unique
length(gwas_genes) # 656

# find GWAS hits that were part of our analysis
gwas_genes_in_universe <- intersect(gene_universe, gwas_genes) %>% unique
length(gwas_genes_in_universe) # 427

# pull significant genes that map to significant GRCCA transcripts
grcca_genes <- df_x_res %>%
    filter(significant == 1 & !is.na(ensembl_gene_id)) %>%
    pull(ensembl_gene_id) %>%
    unique()
length(grcca_genes) # 234

# find genes that are both GWAS hits and GRCCA hits
overlap_genes <- intersect(grcca_genes, gwas_genes_in_universe)
length(overlap_genes) # 4

phyper_gwas <- 1 - phyper(q = length(overlap_genes), 
                          m = length(gwas_genes_in_universe), 
                          n = length(gene_universe) - length(gwas_genes_in_universe), 
                          k = length(grcca_genes)
) 
phyper_gwas # p = 0.095


## GSEA
ranked_gene_list <- df_x_res %>% 
    filter(significant == 1 & !is.na(ensembl_gene_id)) %>% 
    group_by(ensembl_gene_id) %>% 
    slice_max(abs(pearsons_r)) %>% 
    #mutate(pearsons_r = abs(pearsons_r)) %>% 
    arrange(-pearsons_r) %>% 
    dplyr::select(ensembl_gene_id, pearsons_r) %>% 
    distinct() %>% 
    deframe
m_gwas <- df_x_res %>% 
    dplyr::select(ensembl_gene_id) %>% 
    mutate(gwas = ifelse(ensembl_gene_id %in% gwas_genes, "hit", "no_hit"), .before = 1)

set.seed(20240306)
sig_gwas_res <- GSEA(ranked_gene_list, TERM2GENE = m_gwas, scoreType = "std", pvalueCutoff = 1) %>% as_tibble() %>% clean_names
sig_gwas_res

## Plot position of common variants in distribution
df_x_res_gwas <- df_x_res %>% 
    filter(ensembl_gene_id %in% gwas_genes)
df_x_res %>% 
    ggplot(aes(x = z_score)) +
    geom_density() +
    geom_point(data = df_x_res_gwas %>% arrange(significant), shape = 21, size = 3,
               aes(x = z_score, y = 0, fill = z_score, color = as.factor(significant))) +
    geom_vline(xintercept = c(-2, 2), lty = 2) +
    geom_text_repel(data = df_x_res_gwas %>% filter(significant == 1), min.segment.length = 0.1, box.padding = 1, size = 3,
                    aes(x = z_score, y = 0, label = transcript_symbol)) +
    scale_fill_gradientn(colors = rev(brewer.pal(n = 11, "RdBu")), guide = "none") +
    scale_color_manual(values = c("transparent", "yellow"), guide = "none") +
    labs(x = "Transcript weight z-score", y = "", title = "SCZ common variant position in GRCCA results")
map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(figures_dir, "GRCCA_transcripts_GWAS_dist", .x), width = 4, height = 3)
)



# Consensus DEG overlap --------------------------------------------

df_consensus_degs %>% 
    count(direction_of_effect)

## Hypergeometric

# pull genes in transcript-level analysis
gene_universe <- df_x_res %>% 
    filter(!is.na(ensembl_gene_id)) %>% 
    pull(gene_symbol) %>% 
    unique
length(gene_universe) # 17976

# pull consensus DEG hits
consensus_degs <- df_consensus_degs %>% 
    pull(gene_symbol) %>% 
    unique
length(consensus_degs) # 160

# find consensus DEGs that were part of our analysis
consensus_degs_in_universe <- intersect(gene_universe, consensus_degs) %>% unique
length(consensus_degs_in_universe) # 132

# pull significant genes that map to significant GRCCA transcripts
grcca_genes <- df_x_res %>%
    filter(significant == 1 & !is.na(ensembl_gene_id)) %>%
    pull(gene_symbol) %>%
    unique()
length(grcca_genes) # 234

# find genes that are both GWAS hits and GRCCA hits
overlap_genes <- intersect(grcca_genes, consensus_degs_in_universe)
length(overlap_genes) # 2

phyper_consensus_degs <- 1 - phyper(q = length(overlap_genes), 
                          m = length(consensus_degs_in_universe), 
                          n = length(gene_universe) - length(consensus_degs_in_universe), 
                          k = length(grcca_genes)
) 
phyper_consensus_degs # p = 0.2

## GSEA
ranked_gene_list <- df_x_res %>% 
    filter(significant == 1 & !is.na(ensembl_gene_id)) %>% 
    group_by(gene_symbol) %>% 
    slice_max(abs(pearsons_r)) %>% 
    #mutate(pearsons_r = abs(pearsons_r)) %>% 
    arrange(-pearsons_r) %>% 
    dplyr::select(gene_symbol, pearsons_r) %>% 
    distinct() %>% 
    deframe
m_consensus_degs <- df_x_res %>% 
    filter(!is.na(ensembl_gene_id)) %>% 
    dplyr::select(gene_symbol) %>% 
    mutate(deg = ifelse(gene_symbol %in% consensus_degs, "hit", "no_hit"), .before = 1)

set.seed(20240306)
sig_consensus_degs_res <- GSEA(ranked_gene_list, TERM2GENE = m_consensus_degs, scoreType = "std", pvalueCutoff = 1) %>% as_tibble() %>% clean_names
sig_consensus_degs_res

## Plot position of consensus DEGs in distribution
df_x_res_degs <- df_x_res %>% 
    filter(gene_symbol %in% consensus_degs)
df_x_res_degs %>% filter(gene_symbol %in% overlap_genes)
df_x_res %>% 
    ggplot(aes(x = z_score)) +
    geom_density() +
    geom_point(data = df_x_res_degs %>% arrange(significant), shape = 21, size = 3,
               aes(x = z_score, y = 0, fill = z_score, color = as.factor(significant))) +
    geom_vline(xintercept = c(-2, 2), lty = 2) +
    geom_text_repel(data = df_x_res_degs %>% filter(significant == 1), min.segment.length = 0.1, box.padding = 1, size = 3,
                    aes(x = z_score, y = 0, label = transcript_symbol)) +
    scale_fill_gradientn(colors = rev(brewer.pal(n = 11, "RdBu")), guide = "none") +
    scale_color_manual(values = c("transparent", "yellow"), guide = "none") +
    labs(x = "Transcript weight z-score", y = "", title = "Consensus DEG position in GRCCA results")
map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(figures_dir, "GRCCA_transcripts_consensus_DEG_dist", .x), width = 4, height = 3)
)


# Rare variant overlap ---------------------------------------------

## Hypergeometric

# pull genes in transcript-level analysis
gene_universe <- df_x_res %>% 
    filter(!is.na(ensembl_gene_id)) %>% 
    pull(ensembl_gene_id) %>% 
    unique
length(gene_universe) # 18031

# pull consensus DEG hits
exome_urvs <- df_exome_urvs %>% 
    pull(ensembl_gene_id) %>% 
    unique
length(exome_urvs) # 10

# find URVs that were part of our analysis
exome_urvs_in_universe <- intersect(gene_universe, exome_urvs) %>% unique
length(exome_urvs_in_universe) # 8

# pull significant genes that map to significant GRCCA transcripts
grcca_genes <- df_x_res %>%
    filter(significant == 1 & !is.na(ensembl_gene_id)) %>%
    pull(ensembl_gene_id) %>%
    unique()
length(grcca_genes) # 234

# find genes that are both GWAS hits and GRCCA hits
overlap_genes <- intersect(grcca_genes, exome_urvs_in_universe)
length(overlap_genes) # 0

phyper_exome_urvs <- 1 - phyper(q = length(overlap_genes), 
                                    m = length(exome_urvs_in_universe), 
                                    n = length(gene_universe) - length(exome_urvs_in_universe), 
                                    k = length(grcca_genes)
) 
phyper_exome_urvs # p = 0.1

## GSEA
ranked_gene_list <- df_x_res %>% 
    filter(significant == 1 & !is.na(ensembl_gene_id)) %>% 
    group_by(ensembl_gene_id) %>% 
    slice_max(abs(pearsons_r)) %>% 
    #mutate(pearsons_r = abs(pearsons_r)) %>% 
    arrange(-pearsons_r) %>% 
    dplyr::select(ensembl_gene_id, pearsons_r) %>% 
    distinct() %>% 
    deframe
m_exome_urvs <- df_x_res %>% 
    filter(!is.na(ensembl_gene_id)) %>% 
    dplyr::select(ensembl_gene_id) %>% 
    mutate(urv = ifelse(ensembl_gene_id %in% exome_urvs, "hit", "no_hit"), .before = 1)

set.seed(20240306)
sig_exome_urvs_res <- GSEA(ranked_gene_list, TERM2GENE = m_exome_urvs, scoreType = "std", pvalueCutoff = 1) %>% as_tibble() %>% clean_names
sig_exome_urvs_res

## Plot position of consensus DEGs in distribution
df_x_res_urvs <- df_x_res %>% 
    filter(ensembl_gene_id %in% exome_urvs)
df_x_res %>% 
    ggplot(aes(x = z_score)) +
    geom_density() +
    geom_point(data = df_x_res_urvs %>% arrange(significant), shape = 21, size = 3,
               aes(x = z_score, y = 0, fill = z_score, color = as.factor(significant))) +
    geom_vline(xintercept = c(-2, 2), lty = 2) +
    geom_text_repel(data = df_x_res_urvs %>% filter(abs(z_score) > 2), min.segment.length = 0.1, box.padding = 1, size = 3,
                    aes(x = z_score, y = 0, label = transcript_symbol)) +
    scale_fill_gradientn(colors = rev(brewer.pal(n = 11, "RdBu")), guide = "none") +
    scale_color_manual(values = c("transparent", "yellow"), guide = "none") +
    labs(x = "Transcript weight z-score", y = "", title = "Rare variant position in GRCCA results")
map(
    .x = c(".png", ".pdf"),
    .f = ~ ggsave(paste0(figures_dir, "GRCCA_transcripts_exome_urvs_dist", .x), width = 4, height = 3)
)


