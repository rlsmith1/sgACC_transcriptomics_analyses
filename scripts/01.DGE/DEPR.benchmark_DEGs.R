
###################################################################

# Quantify DEG overlap with published SCZ gene lists                        

###################################################################

## Set directories & load data
base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
figures_dir <- paste0(base_dir, "outputs/updated_figures/figures/Fig2_DGEres/")
tables_dir <- paste0(base_dir, "outputs/updated_figures/tables/")
analysis_objects_dir <- paste0(base_dir, "outputs/updated_figures/objects/Fig2_DGEres/")

source(paste0(base_dir, "scripts/full_analysis_scripts/scripts/load_genes.R"))


# Identify SCZ-control DEGs -------------------------------

## DEGs
my_degs <- df_de_res %>% 
    dplyr::filter(pvalue < 0.05) %>% 
    pull(ensembl_gene_id)
length(my_degs) # 1690

### DEGs (ranked) ###
my_degs_ranked <- df_de_res %>% 
    arrange(-log2fold_change) %>% 
    dplyr::select(ensembl_gene_id, log2fold_change) %>% 
    deframe

## Identify max L2FC (for plot color scales)
scale_max <- df_de_res %>% pull(log2fold_change) %>% abs %>% max



# Run hypergeometric test (testing for overall enrichment) on DEGs and each benchmark list -----------------

df_hypergeometric_res <- tibble(
    benchmark_list = names(benchmarking_lists)
) %>% 
    mutate(
        benchmark_list_n = map(
          .x = benchmark_list,
          .f = ~ length(benchmarking_lists[[.x]])
        ),
        benchmark_list_n_in_universe = map(
            .x = benchmark_list,
            .f = ~ intersect(benchmarking_lists[[.x]], ensembl_gene_ids) %>% length
        ),
        hypergeometric_test = map(
            .x = benchmark_list,
            .f = ~ f_hypergeometric(gene_list1 = my_degs, 
                                    gene_list2 = benchmarking_lists[[.x]],
                                    gene_universe = ensembl_gene_ids
            )
        ),
        overlap_n = map(
            .x = hypergeometric_test,
            .f = ~ length(.x$overlap_genes)
        ),
        hyper_p_val = map(
            .x = hypergeometric_test,
            .f = ~ .x$p_value
        )
    ) %>% 
    unnest(cols = c(benchmark_list_n, benchmark_list_n_in_universe, overlap_n, hyper_p_val))



# Run GSEA (testing for pole enrichment) on DEGs and each benchmark list -----------------


df_gsea_res <- tibble(
    benchmark_list = names(benchmarking_lists)
) %>% 
    mutate(
        gsea_mat = map(
            .x = benchmark_list,
            .f = ~ tibble(ensembl_gene_id = ensembl_gene_ids) %>% 
                mutate(hit = ifelse(ensembl_gene_id %in% benchmarking_lists[[.x]], "success", "failure"), .before = 1)
        ),
        gsea_res = map(
            .x = gsea_mat,
            .f = ~ GSEA(my_degs_ranked, TERM2GENE = .x, pvalueCutoff = 1) %>% 
                as_tibble() %>% 
                clean_names
        ),
        gsea_p_val = map(
            .x = gsea_res,
            .f = ~ .x %>% pull(pvalue)
        ),
        gsea_nes = map(
            .x = gsea_res,
            .f = ~ .x %>% pull(nes)
        )
    ) %>%
    unnest(cols = c(gsea_p_val, gsea_nes))



# Combine benchmarking results and save --------------------------------------------


df_benchmarking_res <- df_hypergeometric_res %>% 
    left_join(df_gsea_res)

save(benchmarking_lists, df_benchmarking_res,
     file = paste0(analysis_objects_dir, "DE_benchmark_results.Rdata")
)
