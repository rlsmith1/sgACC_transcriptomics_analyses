
## Load data & functions
base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
figures_dir <- paste0(base_dir, "outputs/updated_figures/figures/Fig5_GRCCAtranscriptsRes/")
tables_dir <- paste0(base_dir, "outputs/updated_figures/tables/")
analysis_objects_dir <- paste0(base_dir, "outputs/updated_figures/objects/Fig5_GRCCAtranscriptsRes/")

source(paste0(base_dir, "scripts/full_analysis_scripts/scripts/load_WTCNA_res.R"))



# load GRCCA results ------------------------------------------------------

## Identify CCA directory
dx_groups_to_analyze <- c("Control", "BD", "MDD", "SCZ") # select diagnostic groups of interest
n_mca_dim <- 8
type <- "OVERALLsubs/"
project_dir <- paste0(prefix, # date & preprocessing info
                      "_sft", soft_power, "_minSize", minimum_size, "_cutHeight", tree_cut_height, "_", # WTCNA parameters
                      paste0(dx_groups_to_analyze, collapse = ""), "_", # diagnostic groups included
                      n_mca_dim, "MCA_regressDim2_", # covariates included
                      #"broadRareSCZ/"
                      "broadRareSCZBDMDDASD_sigGRCCAfdr0_05z2/" # define subset of transcripts included (based on gene-level results)
)
cca_dir <- paste0(base_dir, "RCCA_toolkit/", type, project_dir)

## Read in GRCCA analysis results
analysis_type <- "grcca"
VARx <- "0.1_1"
l_grcca_res <- f_load_cca_res(
    cca_directory = cca_dir,
    level = "transcript",
    analysis_type = analysis_type,
    mu = 0.1,
    lambda = "Nfeat", #"GeneLevel",
    VARx = VARx,
    include_Cmat = TRUE,
    rename_covariates = TRUE # FALSE if newer analysis where disorders are abbreviated, TRUE for old analysis with old labels
)
l_grcca_res

## Assign results to tibbles (negate if SCZ is negative to make SCZ positive)
df_results_transcripts <- l_grcca_res$model_results
df_lvs_transcripts <- l_grcca_res$latent_variables
df_y_res_transcripts <- l_grcca_res$y_res #%>% 
    #mutate_at(vars(c(weight, z_score, pearsons_r)), ~ -.x)
df_x_res_transcripts <- l_grcca_res$x_res %>% 
    #mutate_at(vars(c(weight, z_score, pearsons_r)), ~ -.x) %>% 
    left_join(df_transcript_to_gene, by = join_by(ensembl_transcript_id)) %>% 
    left_join(df_modules_filt %>% dplyr::select(-c(mod_set, color)), by = join_by(ensembl_transcript_id)) %>% 
    dplyr::select(ensembl_transcript_id, transcript_symbol, module, everything()) %>% 
    mutate(transcript_symbol = ifelse(is.na(transcript_symbol), ensembl_transcript_id, transcript_symbol)) %>% 
    distinct()


# Plot y res --------------------------------------------------------------


df_y_res_transcripts %>% 
    pivot_longer(c(z_score, pearsons_r), names_to = "metric", values_to = "value") %>% 
    mutate(metric = ifelse(metric == "z_score", "weight z-score", "structure correlation")) %>% 
    mutate(p_adj = ifelse(p_adj == 0, p_adj + 0.00001, p_adj)) %>% 
    mutate(significance = case_when(
        metric == "weight z-score" & abs(value) > 2 ~ "*",
        metric == "structure correlation" & p_adj < 0.05 & p_adj > 0.01 ~ "*",
        metric == "structure correlation" & p_adj < 0.01 & p_adj > 0.001 ~ "**",
        metric == "structure correlation" & p_adj < 0.001 ~ "***",
        TRUE ~ ""
    )
    ) %>% 
    
    ggplot(aes(x = value, y = reorder_within(covariate, value, metric))) +
    geom_point(aes(size = -log10(p_adj), fill = value), shape = 21) +
    geom_vline(xintercept = 0, lty = 2) +
    geom_text(aes(label = significance), size = 4, vjust = -0.15, color = "black") +
    facet_wrap(vars(metric), scales = "free", nrow = 1) +
    
    scale_fill_gradientn(colors = gene_weight_color_scale, limits = c(-3, 3)) +
    scale_size_continuous(range = c(2, 6)) +
    scale_y_reordered() +
    guides(fill = guide_colorbar(title = NULL),
           size = guide_legend(title = "-log10(FDR)")) +
    coord_cartesian(clip = "off") +
    labs(y = NULL, x = NULL, 
         title = "Covariate results") +
    theme(legend.position = "bottom", 
          legend.justification = "center",
          legend.box = "horizontal",
          #legend.justification = "center",
          #legend.key.size = unit(0.2, 'cm')
          legend.key.height = unit(0.20, "cm"),
          legend.key.width = unit(1.0, "cm")
    )

## Plot x-res



# Hypergeometric test with risk genes -------------------------------------


## Identify all genes included in transcript-level analysis
gene_universe <- df_x_res_transcripts %>% 
    #filter(pearsons_r > 0) %>% 
    pull(ensembl_gene_id) %>% 
    unique
length(gene_universe) # 1719

## Map significant transcripts to parent genes
sig_grcca_genes <- df_x_res_transcripts %>% 
    #filter(pearsons_r > 0) %>% 
    filter(p_adj < 0.05 & abs(z_score) >= 2) %>% 
    pull(ensembl_gene_id) %>% 
    unique
length(sig_grcca_genes) # 209


map(
    .x = benchmarking_lists,
    .f = ~ f_hypergeometric(
        gene_list1 = sig_grcca_genes,
        gene_list2 = .x,
        gene_universe = gene_universe
    )
)


# GSEA res --------------------------------------------------------


## Generate vector
transcript_grcca_ranked <- df_x_res_transcripts %>% 
    group_by(ensembl_gene_id) %>% 
    summarise(pearsons_r = median(pearsons_r)) %>% 
    arrange(-pearsons_r) %>% 
    dplyr::select(ensembl_gene_id, pearsons_r) %>% 
    deframe

## Benchmark
fgsea(pathways = benchmarking_lists, 
      stats = transcript_grcca_ranked, 
      eps = 0
) %>% 
    as_tibble() %>% 
    dplyr::rename("benchmark_list" = "pathway") %>% 
    arrange(-abs(NES))

## Cell type
fgsea(pathways = cell_types, 
      stats = transcript_grcca_ranked, 
      eps = 0
)  %>% 
    as_tibble() %>% 
    dplyr::rename("cell_type" = "pathway") %>% 
    arrange(-abs(NES))

## GO
fgsea(pathways = go_pathways, 
      stats = transcript_grcca_ranked, 
      eps = 0
) %>% 
    as_tibble() %>% 
    dplyr::rename("go_id" = "pathway") %>% 
    arrange(-abs(NES)) %>% 
    left_join(df_go_terms) %>% # add names of GO terms
    dplyr::select(go_id, term, ontology, everything())



# Compare with gene-level GRCCA -------------------------------------------

load(paste0(base_dir, "outputs/objects//Fig3_GRCCAres/GRCCA_results.Rdata")) # df_x_res
df_x_res_genes <- df_x_res
rm(list = "df_x_res")

## Identify genes included in this transcript analysis
grcca_genes_to_include <- df_x_res_genes %>% # significant GRCCA genes only
    filter(abs(z_score) >= 1.96) %>%
    pull(ensembl_gene_id) %>%
    unique
#risk_genes_to_include <- c( benchmarking_lists[["SCZ (common, broad)"]], benchmarking_lists[["SCZ (rare)"]] ) # SCZ risk genes (common and rare)
risk_genes_to_include <- unlist(benchmarking_lists) %>% unique # all psychiatric risk genes
genes_to_include <- unique(c(grcca_genes_to_include, risk_genes_to_include))

## Check correspondence between gene and transcript-level results
df_x_res_transcripts %>% 
    filter(!is.na(ensembl_gene_id)) %>% 
    group_by(ensembl_gene_id) %>% 
    summarise(pearsons_r = max(pearsons_r)) %>% 
    dplyr::rename("transcript_res" = "pearsons_r") %>% 
    
    left_join(df_x_res_genes %>% 
                  dplyr::select(ensembl_gene_id, gene_symbol, pearsons_r) %>% 
                  dplyr::rename("gene_res" = "pearsons_r")
    ) %>% 
    mutate(
        risk_gene = ifelse(ensembl_gene_id %in% benchmarking_lists[["SCZ (common, broad)"]], "yes", "no"),
        label = ifelse(risk_gene == "yes", gene_symbol, "")
    ) %>% 
    arrange(risk_gene) %>% 
    
    ggplot(aes(x = gene_res, y = transcript_res)) +
    geom_vline(xintercept = 0, color = "gray") +
    geom_hline(yintercept = 0, color = "gray") +
    geom_abline(lty = 2) +
    geom_point(aes(fill = gene_res + transcript_res, color = risk_gene, size = risk_gene), shape = 21) +
    geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 10) +
    stat_cor() +
    geom_smooth(method = "lm", color = "black", linewidth = 0.75) +
    xlim(c(-0.6, 0.6)) +
    ylim(c(-0.6, 0.6)) +
    scale_fill_gradientn(colors = gene_weight_color_scale, guide = "none") +
    scale_color_manual(values = c("yes" = "black", "no" = "transparent")) +
    scale_size_manual(values = c("yes" = 1.5, "no" = 0.5)) +
    labs(x = "Gene GRCCA structure correlations", y = "Transcript GRCCA structure correlations",
         title = "Correspondence between gene & \ntranscript GRCCA results")

## Identify genes with more than one transcript in the transcript-level analysis
more_than_one_transcript <- df_x_res_transcripts %>% 
    dplyr::select(contains("ensembl")) %>% 
    count(ensembl_gene_id, sort = TRUE) %>% 
    filter(n > 1) %>% 
    pull(ensembl_gene_id)

## Check directionality of significant risk genes (from GRCCA results)
sig_risk_genes <- df_x_res_genes %>% 
    filter(ensembl_gene_id %in% intersect(genes_to_include, more_than_one_transcript) & p_adj < 0.05 & abs(z_score) > 2) %>%
    top_n(10, abs(pearsons_r)) %>% 
    pull(ensembl_gene_id)
length(sig_risk_genes)
sig_risk_gene_order <- df_x_res_genes %>% 
    filter(ensembl_gene_id %in% sig_risk_genes) %>% 
    arrange(-pearsons_r) %>% 
    pull(gene_symbol)

df_x_res_transcripts %>% 
    filter(ensembl_gene_id %in% sig_risk_genes) %>% 
    #filter(ensembl_gene_id %in% intersect(sig_risk_genes, more_than_one_transcript)) %>%
    mutate(
        gene_symbol = factor(gene_symbol, levels = sig_risk_gene_order),
        significant = ifelse(p_adj < 0.05, "yes", "no")
    ) %>% 
    
    ggplot(aes(x = pearsons_r, y = reorder(transcript_symbol, pearsons_r))) +
    geom_vline(xintercept = 0, color = "gray") +
    geom_point(aes(fill = pearsons_r, size = -log10(p_adj), color = significant), shape = 21) +
    facet_wrap(vars(gene_symbol), scales = "free_y") +
    scale_fill_gradientn(colors = gene_weight_color_scale) +
    scale_color_manual(values = c("yes" = "black", "no" = "gray")) +
    coord_cartesian(clip = "off")


## Check directionality of non-significant risk genes (from GRCCA results)
non_sig_risk_genes <- df_x_res_genes %>% 
    filter(ensembl_gene_id %in% intersect(genes_to_include, more_than_one_transcript) & p_adj > 0.05 & abs(z_score) < 0.1) %>%
    top_n(10, -abs(pearsons_r)) %>% 
    pull(ensembl_gene_id)
length(non_sig_risk_genes)
non_sig_risk_gene_order <- df_x_res_genes %>% 
    filter(ensembl_gene_id %in% non_sig_risk_genes) %>% 
    arrange(-pearsons_r) %>% 
    pull(gene_symbol) %>% 
    unique

df_x_res_transcripts %>% 
    filter(ensembl_gene_id %in% non_sig_risk_genes) %>%
    mutate(
        gene_symbol = factor(gene_symbol, levels = non_sig_risk_gene_order),
        significant = ifelse(p_adj < 0.05, "yes", "no")
    ) %>% 
    
    ggplot(aes(x = pearsons_r, y = reorder(transcript_symbol, pearsons_r))) +
    geom_vline(xintercept = 0, color = "gray") +
    geom_point(aes(fill = pearsons_r, size = -log10(p_adj), color = significant), shape = 21) +
    facet_wrap(vars(gene_symbol), scales = "free_y") +
    scale_fill_gradientn(colors = gene_weight_color_scale) +
    scale_color_manual(values = c("yes" = "black", "no" = "gray")) +
    coord_cartesian(clip = "off")


## Weight transcripts of same gene by relative expression?
example_transcripts <- df_transcript_to_gene %>% 
    filter(gene_symbol == "SLC25A12") %>% 
    pull(ensembl_transcript_id) %>% 
    unique
df_example_transcripts <- df_vsd_regress_filt %>% 
    dplyr::select(sample, all_of(example_transcripts))
df_example_transcripts %>% 
    pivot_longer(contains("ENST"), names_to = "ensembl_transcript_id", values_to = "value") %>% 
    left_join(df_transcript_to_gene) %>% 
    left_join(df_covariates) %>% 
    distinct() %>% 
    
    ggplot(aes(x = value, y = transcript_symbol)) +
    #geom_point() +
    geom_boxplot() +
    scale_color_manual(values = dx_colors)

