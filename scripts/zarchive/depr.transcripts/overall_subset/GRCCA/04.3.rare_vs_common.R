
########################################################################################

# OVERALL SUBSET transcript level 
# characterize role of rare and common transcripts in GRCCA results

########################################################################################

# libraries ---------------------------------------------------------------

library(tidyverse)
library(patchwork)
library(janitor)
library(RColorBrewer)
library(biomaRt)
library(tidymodels)  
library(ggrepel)
library(ggpubr)
library(clusterProfiler)
library(msigdbr)
library(org.Hs.eg.db)
library(ggh4x)
library(xlsx)
library(readxl)


#library(ggchicklet)


# set theme for plots -----------------------------------------------------

theme_set(theme_bw() +
            theme(plot.title = element_text(size = 11),
                  axis.title = element_text(size = 10),
                  axis.text = element_text(size = 10),
                  strip.text = element_text(size = 10),
                  legend.title = element_text(size = 8),
                  legend.text = element_text(size = 8)
            )
)

# DX COLORS
dx_colors <- c("#0072B2", "#E69F00", "#009E73", "#9966FF")
names(dx_colors) <- c("Control", "BD", "MDD", "SCZ")


# data --------------------------------------------------------------------

base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
prefix <- "08Mar2024_TRANSCRIPTS_qSVAgeSexRaceGC_CVq1"

# LOAD OBJECTS 
load(paste0(base_dir, "objects/14Nov2023_covariates.Rdata")) # df_covariates, df_covariates_numeric (generated in clean_covariates.R)
load(paste0(base_dir, "objects/", prefix, ".RDS")) # df_vsd_regress

# TRANSCRIPT SYMBOL TO ENSEMBL TRANSCRIPT ID MAPPING
load(paste0(base_dir, "objects/hsapiens_genome_v110.RDS")) # df_hsapiens_genome

df_gene_to_transcript <- df_hsapiens_genome %>% 
  dplyr::select(transcript_id, transcript_name, gene_id, gene_name) %>% 
  distinct() %>% 
  dplyr::rename("gene_symbol" = "gene_name", "transcript_symbol" = "transcript_name", 
                "ensembl_gene_id" = "gene_id", "ensembl_transcript_id" = "transcript_id"
  ) %>% 
  filter(!is.na(ensembl_transcript_id)) %>% 
  mutate(transcript_symbol = ifelse(is.na(transcript_symbol), ensembl_transcript_id, transcript_symbol),
         gene_symbol = ifelse(is.na(gene_symbol), ensembl_gene_id, gene_symbol)
  )

df_ensembl_to_symbol <- df_gene_to_transcript %>% 
  dplyr::select(ensembl_transcript_id, transcript_symbol) %>% 
  distinct()

# LOAD TRANSCRIPT GRCCA RESULTS
df_x_res <- read_xlsx(paste0(base_dir, "outputs/tables/for_manuscript/TRANSCRIPTS_GRCCA_res_withGray.xlsx"), sheet = 1) %>% 
  mutate(module = factor(module, levels = 0:max(as.numeric(module))))
df_y_res <- read_excel(paste0(base_dir, "outputs/tables/for_manuscript/TRANSCRIPTS_GRCCA_res_withGray.xlsx"), sheet = 2)

# LOAD GENE X RES
df_x_res_genes <- read_excel(paste0(base_dir, "outputs/tables/for_manuscript/GENES_GRCCA_res_withGray.xlsx"), sheet = 1) %>% 
    mutate(module = factor(module, levels = 0:max(as.numeric(module))))

# TRANSCRIPT RAW COUNTS
df_transcript_raw_counts <- read.csv(paste0(base_dir, "data/transcript_count_matrix.csv")) %>% 
  as_tibble() %>% 
  dplyr::rename("ensembl_transcript_id" = "X") %>% 
  rename_all(~str_remove(.x, "_MergedBam_stringtieOutput")) %>% 
  clean_names() %>% 
  rename_all(~str_remove(.x, "x")) 



# 5A | are implicated transcripts rare or common? ------------------------------


# identify all transcripts included in the GRCCA analysis
transcript_universe <- df_x_res %>% 
  pull(ensembl_transcript_id) 
length(transcript_universe) # 54302

# transcript identified as significant in GRCCA analysis
grcca_transcripts <- df_x_res %>%
  filter(p_adj < 0.05 & abs(z_score) >= 2) %>%
  pull(ensembl_transcript_id) 
length(grcca_transcripts) # 257

# IDENTIFY RARE TRANSCRIPTS (At least 10 counts and less than 100 counts in at least 80% of samples = rare)

rare_transcripts <- df_transcript_raw_counts %>% 
  filter(ensembl_transcript_id %in% transcript_universe) %>% 
  mutate_if(is.numeric, ~ifelse(.x >= 10 & .x < 100, 1, 0)) %>% 
  mutate(thresh = rowSums(dplyr::select(., -ensembl_transcript_id))/(ncol(.) - 1), .before = 2) %>% 
  filter(thresh >= 0.80) %>% 
  pull(ensembl_transcript_id) 
length(rare_transcripts) # 23211

# HYPERGEOMETRIC TEST
rare_grcca_intersect <- intersect(rare_transcripts, grcca_transcripts) # 2802

phyper(q = length(rare_grcca_intersect), 
       m = length(grcca_transcripts), 
       n = length(transcript_universe) - length(grcca_transcripts), 
       k = length(rare_transcripts),
       lower.tail = FALSE
) # rare transcripts are not overrepresented in these results


# PLOT COMMON VS RARE
df_x_res %>%
  mutate(prevalence = ifelse(ensembl_transcript_id %in% rare_transcripts, "rare", "common")) %>% 
  mutate(significant = ifelse(significant == 1, "significant", "not significant")) %>% 
  
  ggplot(aes(x = significant, fill = prevalence)) +
  geom_histogram(stat = "count") +
  guides(fill = guide_legend(title = "")) +
    facet_wrap(vars(significant), scales = "free") +
  labs(x = "", y = "Number of transcripts",
       title = bquote(bold("A") ~ " | Prevalence of significant transcripts")) +
  theme(legend.position = c(0.8, 0.8))

map(
  .x = c(".png", ".pdf"),
  .f = ~ ggsave(paste0("~/Documents/PhD/manuscripts/BiolPsych2024/figures/main_text/Fig4A", .x),
                width = 4, height = 3)
)



# supplement? | module prevalence of rare and common transcripts ------------------------

# hypergeometric test (module 8 looks significant maybe)

# q = number of rare transcripts in the module
# m = total rare transcripts in universe
# n = transcripts in universe that were *not* rare (total transcripts - rare transcripts)
# k = total transcripts in module

df_x_res %>%
  dplyr::select(ensembl_transcript_id, module) %>% 
  count(module) %>% 
  dplyr::rename("k" = "n") %>% 
  left_join(
    df_x_res %>%
      dplyr::select(ensembl_transcript_id, module) %>% 
      filter(ensembl_transcript_id %in% rare_transcripts) %>% 
      count(module) %>% 
      dplyr::rename("q" = "n"),
    by = join_by(module)
  ) %>% 
  mutate(
    m = length(rare_transcripts),
    n = length(transcript_universe) - m,
    p_value = phyper(q = q, 
                     m = m, 
                     n = n, 
                     k = k,
                     lower.tail = FALSE
    ),
    p_adj = p.adjust(p_value, method = "fdr")
  ) %>% 
  arrange(p_value) # only the gray module is enriched for rare transcripts (makes sense)


# 5B-E | rare vs common GSEA ---------------------------------------------------------------

# plotting function
facet_plot <- function(data, fill_scale_high = "black", legend_position = "none") {
  data %>% 
    ggplot(aes(x = reorder(str_wrap(description, 20), pvalue), y = -log10(p_adjust))) +
    geom_hline(aes(yintercept = -log10(0.05)), lty = 2, color = "black") +
    geom_point(aes(fill = -log10(pvalue), size = set_size), shape = 21) +
    facet_wrap(vars(prevalence)) +
    scale_size_continuous(range = c(3, 6), limits = c(10, 460)) +
    scale_fill_gradient2(low = colorRampPalette(c("white", fill_scale_high))(10)[3], high = fill_scale_high) +
    guides(size = guide_legend(title = "n genes"),
           fill = FALSE) +
    labs(x = "", y = "-log10(FDR)") +
    ylim(c(0, 10)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
          legend.position = legend_position)
}

# DEFINE SIGNIFICANT GENE LISTS ASSOCIATED WITH RARE AND COMMON TRANSCRIPTS

rare_gene_list <- df_x_res %>%
    filter(p_adj < 0.05 & abs(z_score) >= 2 & ensembl_transcript_id %in% rare_transcripts) %>% 
    arrange(-pearsons_r) %>%
    dplyr::select(ensembl_gene_id, pearsons_r) %>% 
    group_by(ensembl_gene_id) %>% 
    summarise(pearsons_r = median(pearsons_r)) %>% 
    arrange(-pearsons_r) %>%
    deframe 
length(rare_gene_list) # n = 84

common_gene_list <- df_x_res %>%
    left_join(df_gene_to_transcript) %>% 
    filter(p_adj < 0.05 & abs(z_score) >= 2 & !(ensembl_transcript_id %in% rare_transcripts)) %>% 
    arrange(-pearsons_r) %>%
    dplyr::select(ensembl_gene_id, pearsons_r) %>% 
    group_by(ensembl_gene_id) %>% 
    summarise(pearsons_r = median(pearsons_r)) %>% 
    arrange(-pearsons_r) %>%
    deframe 
length(common_gene_list) # n = 159

# B | GENE ONTOLOGY GSEA (BIOLOGICAL PROCESS, CELLULAR COMPONENT, MOLECULAR FUNCTION)
set.seed(20240319)
df_combos <- expand_grid(
    prevalence = c("rare", "common"),
    ontology = c("BP", "CC", "MF")
)

df_go_res_rareVcommon <- map2_dfr(
    .x = df_combos$prevalence,
    .y = df_combos$ontology,
    .f = ~ gseGO(geneList = if(.x == "rare") {rare_gene_list} else{common_gene_list},
                 ont = .y,
                 keyType = "ENSEMBL",
                 OrgDb = "org.Hs.eg.db",
                 pvalueCutoff = 1,
                 eps = 1e-300
    ) %>% 
        as_tibble() %>% 
        clean_names %>% 
        mutate(prevalence = .x, ontology = .y, .before = 1)
)

# plot all
df_go_res_rareVcommon %>% 
    group_by(prevalence, ontology) %>% 
    top_n(n = 10, wt = -pvalue) %>%
    mutate(`survives FDR` = ifelse(p_adjust < 0.05, "yes", "no")) %>% 
    
    ggplot(aes(y = reorder(str_wrap(description, width = 35), -pvalue), x = -log10(pvalue))) +
    geom_vline(aes(xintercept = -log10(0.05)), lty = 2) +
    geom_point(aes(size = set_size, fill = ontology, color = `survives FDR`), shape = 21) +
    facet_grid(ontology ~ prevalence, scales = "free_y") +
    scale_fill_manual(values = c("#F8766D", "#00BA38", "#619CFF"), guide = "none") +
    scale_color_manual(values = c("yes" = "black", "no" = "transparent")) +
    scale_size_continuous(range = c(3, 6), limits = c(10, 500)) +
    guides(size = guide_legend(title = "n genes")) +
    labs(y = "", x = "-log10(p-value)", title = bquote(bold("E")~" | Gene set ontology enrichments")) +
    theme(legend.position = c(0.75, 0.15),
          legend.box = "horizontal")

# CELL TYPE  
rare_cell_res <- GSEA(rare_gene_list, TERM2GENE = m_cell, pvalueCutoff = 1) %>% as_tibble()
common_cell_res <- GSEA(common_gene_list, TERM2GENE = m_cell, pvalueCutoff = 1) %>% as_tibble()

df_cell_gsea <- rare_cell_res %>% 
    mutate(prevalence = "rare", .before = 1) %>% 
    bind_rows(common_cell_res %>% 
                  mutate(prevalence = "common", .before = 1)
    ) %>% 
    clean_names %>% 
    group_by(prevalence)

p_cell_gsea <- map2(
    .x = c("rare", "common"),
    .y = c("#00BFC4", "#F8766D"),
    .f = ~ facet_plot(
        df_cell_gsea %>% filter(prevalence == .x & pvalue < 0.05), fill_scale_high = .y, legend_position = "none"
    )
) %>% 
    wrap_plots(widths = c(2, 4)) +
    plot_annotation(title = bquote(bold("E") ~ " | Cell type"))




#### ****************EXTRA*****************

# POSITIONAL
rare_pge_res <- GSEA(rare_gene_list, TERM2GENE = m_position, pvalueCutoff = 1) %>% as_tibble()
common_pge_res <- GSEA(common_gene_list, TERM2GENE = m_position, pvalueCutoff = 1) %>% as_tibble()

df_position_gsea <- rare_pge_res %>% 
  mutate(prevalence = "rare", .before = 1) %>% 
  bind_rows(common_pge_res %>% 
              mutate(prevalence = "common", .before = 1)
  ) %>% 
  clean_names %>% 
  group_by(prevalence)

p_position_gsea <- map2(
  .x = c("rare", "common"),
  .y = c("#00BFC4", "#F8766D"),
  .f = ~ facet_plot(
    df_position_gsea %>% filter(prevalence == .x & pvalue < 0.05), fill_scale_high = .y, legend_position = "none"
  )
) %>% 
  wrap_plots(widths = c(5, 15)) +
  plot_annotation(title = bquote(bold("D") ~ " | Chromosome position"))

# BIOTYPE
rare_biotype_res <- GSEA(rare_gene_list, TERM2GENE = m_biotype, pvalueCutoff = 1) %>% as_tibble()
common_biotype_res <- GSEA(common_gene_list, TERM2GENE = m_biotype, pvalueCutoff = 1) %>% as_tibble()

df_biotype_gsea <- rare_biotype_res %>% 
  mutate(prevalence = "rare", .before = 1) %>% 
  bind_rows(common_biotype_res %>% 
              mutate(prevalence = "common", .before = 1)
  ) %>% 
  clean_names %>% 
  group_by(prevalence)

p_biotype_gsea <- map2(
  .x = c("rare", "common"),
  .y = c("#00BFC4", "#F8766D"),
  .f = ~ facet_plot(
    df_biotype_gsea %>% filter(prevalence == .x & pvalue < 0.05) %>% mutate(description = str_replace_all(description, "_", " ")), 
    fill_scale_high = .y, legend_position = "none"
  )
) %>% 
  wrap_plots(widths = c(2, 3)) +
  plot_annotation(title = bquote(bold("F") ~ " | Biotype"))


# COMBINE AND SAVE
p_functional_gsea /
  p_position_gsea /
  (p_biotype_gsea | p_cell_gsea)


map(
  .x = c(".png", ".pdf"),
  .f = ~ ggsave(paste0("~/Documents/PhD/manuscripts/BiolPsych2024/figures/main_text/Fig4BCDE", .x),
                width = 12, height = 10)
)



# GO enrichment (too few genes for reliable enrichments) ------------------

# DEFINE SIGNIFICANT GENE LISTS ASSOCIATED WITH RARE AND COMMON TRANSCRIPTS

rare_gene_list <- df_x_res %>%
    filter(p_adj < 0.05 & abs(z_score) >= 2 & ensembl_transcript_id %in% rare_transcripts) %>% 
    arrange(-pearsons_r) %>%
    dplyr::select(ensembl_gene_id, pearsons_r) %>% 
    group_by(ensembl_gene_id) %>% 
    summarise(pearsons_r = median(pearsons_r)) %>% 
    arrange(-pearsons_r) %>%
    deframe 
length(rare_gene_list) # n = 84

common_gene_list <- df_x_res %>%
    left_join(df_gene_to_transcript) %>% 
    filter(p_adj < 0.05 & abs(z_score) >= 2 & !(ensembl_transcript_id %in% rare_transcripts)) %>% 
    arrange(-pearsons_r) %>%
    dplyr::select(ensembl_gene_id, pearsons_r) %>% 
    group_by(ensembl_gene_id) %>% 
    summarise(pearsons_r = median(pearsons_r)) %>% 
    arrange(-pearsons_r) %>%
    deframe 
length(common_gene_list) # n = 159

# RUN GO ON EACH LIST
df_go_res <- map_dfr(
    .x = c("rare", "common"),
    .f = ~ enrichGO(gene = if(.x == "rare") {names(rare_gene_list)} else{names(common_gene_list)},
                    OrgDb = "org.Hs.eg.db",
                    universe = df_x_res %>% pull(ensembl_gene_id) %>% unique,
                    keyType = "ENSEMBL",
                    ont = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 1,
                    qvalueCutoff = 1,
                    readable = TRUE
    ) %>% 
        as_tibble() %>% 
        clean_names %>% 
        mutate(prevalence = .x, .before = 1)
)

# PLOT RESULTS

# function
f_plot_go_res <- function(go_res) {
    go_res %>% 
        mutate(`survives FDR` = ifelse(p_adjust < 0.05, "yes", "no")) %>% 
        mutate(description = str_wrap(description, width = 35)) %>%
        filter(pvalue < 0.05) %>% 
        group_by(ontology) %>% top_n(10, wt = -pvalue) %>% 
        mutate(gene_ratio = eval(parse(text = gene_ratio))) %>% 
        
        # plot
        ggplot(aes(x = -log10(pvalue), y = reorder(description, -pvalue))) +
        geom_point(aes(size = gene_ratio, fill = ontology, color = `survives FDR`),
                   shape = 21, stroke = 1) +
        geom_vline(aes(xintercept = -log10(0.05)), lty = 2, color = "black") +
        scale_color_manual(values = c("yes" = "black", "no" = "transparent")) +
        scale_size_continuous(range = c(2, 4)) +
        scale_y_reordered() +
        facet_wrap(vars(ontology), scales = "free_y") +
        labs(y = "", x = "log10(p-value)") +
        guides(size = guide_legend(title = "Gene ratio")) +
        theme(legend.position = "right",
              legend.direction = "horizontal",
              legend.box = "vertical",
              legend.justification = "center",
              legend.key.size = unit(0.3, "cm")
        )
}

# plot
l_go_res_plots <- map(
    .x = c("rare", "common"),
    .f = ~ f_plot_go_res(df_go_res %>% filter(prevalence == .x))
)
wrap_plots(l_go_res_plots)


# dig into individual transcripts from GSEA -------------------------------------

rare_astrocyte_transcripts <- rare_cell_res %>% 
  filter(ID == "Astro") %>% 
  pull(core_enrichment) %>% 
  str_split(pattern = "/") %>% .[[1]]

df_gene_to_transcript %>% 
  left_join(df_gene_to_transcript) %>% 
  filter(ensembl_gene_id %in% rare_astrocyte_transcripts) %>% 
  count(gene_symbol) %>% 
  print(n = nrow(.))

rare_pge_res

# they are genes :(

rare_17q_transcripts <- rare_pge_res %>% 
  filter(ID == "chr17q21") %>% 
  pull(core_enrichment) %>% 
  str_split(pattern = "/") %>% .[[1]] # region in Eva's paper, associated with SA

df_gene_to_transcript %>% 
  left_join(df_gene_to_transcript %>% dplyr::select(contains("gene")) %>% distinct()) %>% 
  filter(ensembl_gene_id %in% rare_17q_transcripts) %>% 
  count(gene_symbol) %>% 
  print(n = nrow(.)) # CRHR1 also in Eva's paper

rare_14q_transcripts <- rare_pge_res %>% 
  filter(ID == "chr14q32") %>% 
  pull(core_enrichment) %>% 
  str_split(pattern = "/") %>% .[[1]] # region in Eva's paper
# genes contributing to covariation between CT and schizophrenia were enriched at chromosome 14q32

df_gene_to_transcript %>% 
  left_join(df_gene_to_transcript %>% dplyr::select(contains("gene")) %>% distinct()) %>% 
  filter(ensembl_gene_id %in% rare_14q_transcripts) %>% 
  count(gene_symbol) %>% 
  print(n = nrow(.))

# Eva's paper
df_x_res %>% mutate(rare = ifelse(ensembl_transcript_id %in% rare_transcripts, 1, 0)) %>% filter(str_detect(transcript_symbol, "ARHGAP27"))
df_x_res %>% mutate(rare = ifelse(ensembl_transcript_id %in% rare_transcripts, 1, 0)) %>% filter(str_detect(transcript_symbol, "CRHR1"))
df_x_res %>% mutate(rare = ifelse(ensembl_transcript_id %in% rare_transcripts, 1, 0)) %>% filter(str_detect(transcript_symbol, "DBF4B"))
df_x_res %>% mutate(rare = ifelse(ensembl_transcript_id %in% rare_transcripts, 1, 0)) %>% filter(str_detect(transcript_symbol, "FMNL1"))
df_x_res %>% mutate(rare = ifelse(ensembl_transcript_id %in% rare_transcripts, 1, 0)) %>% filter(str_detect(transcript_symbol, "MAPT"))
df_x_res %>% mutate(rare = ifelse(ensembl_transcript_id %in% rare_transcripts, 1, 0)) %>% filter(str_detect(transcript_symbol, "SPPL2C"))
df_x_res %>% mutate(rare = ifelse(ensembl_transcript_id %in% rare_transcripts, 1, 0)) %>% filter(str_detect(transcript_symbol, "MEIOC"))

df_x_res %>% mutate(rare = ifelse(ensembl_transcript_id %in% rare_transcripts, 1, 0)) %>% filter(str_detect(transcript_symbol, "SRRM2"))

common_6q_transcripts <- common_pge_res %>% 
  filter(ID == "chr6q25") %>% 
  pull(core_enrichment) %>% 
  str_split(pattern = "/") %>% .[[1]] # region in Eva's paper
# genes contributing to covariation between CT and schizophrenia were enriched at chromosome 14q32

df_gene_to_transcript %>% 
  left_join(df_gene_to_transcript %>% dplyr::select(contains("gene")) %>% distinct()) %>% 
  filter(ensembl_gene_id %in% common_6q_transcripts) %>% 
  count(gene_symbol) %>% 
  print(n = nrow(.))

# overlapping genes between schizophrenia and MRI metrics
df_s7 <- read.xlsx(paste0(base_dir, "data/stauffer_2023_supplementary_tables.xlsx"), sheetName = "S7") %>% 
  as_tibble() %>% 
  row_to_names(1) # overlapping genes between schizophrenia and MRI metrics, concentrated in 17q21 and 3p21

df_x_res %>% 
  mutate(rare = ifelse(ensembl_transcript_id %in% rare_transcripts, 1, 0)) %>% 
  filter(str_detect(transcript_symbol, paste0(df_s7 %>% filter(CHR == 3) %>% pull(`Gene name`), collapse = "|"))) %>% 
  filter(significant == 1)

# Genes with highest contribution to covariance between brain structure and schizophrenia
df_s14 <- read.xlsx(paste0(base_dir, "data/stauffer_2023_supplementary_tables.xlsx"), sheetName = "S14") %>% 
  as_tibble() %>% 
  row_to_names(2)

df_s14 <- df_s14[1:2] %>% 
  bind_rows(df_s14[4:5]) %>% 
  bind_rows(df_s14[7:8]) %>% 
  clean_names

df_ct_scz_transcripts <- df_x_res %>% 
  mutate(rare = ifelse(ensembl_transcript_id %in% rare_transcripts, 1, 0)) %>% 
  filter(str_detect(transcript_symbol, paste0(df_s14 %>% filter(mri_metric == "CT") %>% pull(gene_symbol), collapse = "|"))) %>% 
  filter(significant == 1)
df_ct_scz_transcripts %>% count(rare)

df_sa_scz_transcripts <- df_x_res %>% 
  mutate(rare = ifelse(ensembl_transcript_id %in% rare_transcripts, 1, 0)) %>% 
  filter(str_detect(transcript_symbol, paste0(df_s14 %>% filter(mri_metric == "SA") %>% pull(gene_symbol), collapse = "|"))) %>% 
  filter(significant == 1)
df_sa_scz_transcripts %>% count(rare)

df_ndi_scz_transcripts <- df_x_res %>% 
  mutate(rare = ifelse(ensembl_transcript_id %in% rare_transcripts, 1, 0)) %>% 
  filter(str_detect(transcript_symbol, paste0(df_s14 %>% filter(mri_metric == "NDI") %>% pull(gene_symbol), collapse = "|"))) %>% 
  filter(significant == 1)
df_ndi_scz_transcripts %>% count(rare)

