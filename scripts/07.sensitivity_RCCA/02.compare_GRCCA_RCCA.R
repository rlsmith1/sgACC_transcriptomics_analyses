
################################################################################

# GENE-LEVEL: Compare RCCA vs GRCCA results in the context of the SCZ literature

################################################################################

## Load data & functions
base_dir <- "~/Documents/PhD/projects/sgacc_wgcna_grcca/"
figures_dir <- paste0(base_dir, "outputs/figures/supplement/")
#tables_dir <- paste0(base_dir, "outputs/tables/")
analysis_objects_dir <- paste0(base_dir, "outputs/objects/")

source(paste0(base_dir, "scripts/load_WGCNA_res.R"))


## Load RCCA results
load(paste0(base_dir, "objects/RCCA_results.Rdata"))

## Also load GRCCA results for comparison
GRCCAres_dir <- paste0(base_dir, "outputs/updated_figures/objects/Fig3_GRCCAres/")
GRCCAenrich_dir <- paste0(base_dir, "outputs/updated_figures/objects/Fig4_characterizeGRCCA/")
for (dir in c(GRCCAres_dir, GRCCAenrich_dir)) {
    objects <- list.files(dir)
    objects <- objects[!str_detect(objects, "null|archive")]
    for (obj in objects){
        print(paste0("loading ", obj, "...."))
        load(paste0(dir, obj))
    }
}



# Combine X & Y results in GRCCA & RCCA -----------------------------------


## Y
df_y_res_all <- df_y_res_rcca %>% 
    mutate(analysis = "rcca", .before = 1) %>% 
    bind_rows(
        df_y_res %>% 
            mutate(analysis = "grcca", .before = 1)
        
    )

## X
df_x_res_all <- df_x_res_rcca %>% 
    mutate(analysis = "rcca", .before = 1) %>% 
    bind_rows(
        df_x_res %>% 
            mutate(analysis = "grcca", .before = 1)
        
    )


# Test hypergeometric overlap of GRCCA and RCCA hits across thresholds with each benchmark list ------------------


### Loop through and run hypergeometric tests across structure correlation thresholds ###

thresholds <- seq(0.01, 0.25, by = 0.01)
df_rcca_grcca_hypergeometric_benchmark <- expand_grid(
    analysis = c("rcca", "grcca"),
    benchmark_list = names(benchmarking_lists),
    threshold = thresholds
) %>% 
    mutate(
        hypergeometric_res = pmap(
            .l = list(..1 = analysis, ..2 = threshold, ..3 = benchmark_list),
            .f = ~ f_hypergeometric(gene_list1 = df_x_res_all %>% 
                                        dplyr::filter(analysis == ..1) %>% 
                                        top_n(n = ..2*length(ensembl_gene_ids), wt = abs(pearsons_r)) %>% 
                                        pull(ensembl_gene_id), 
                                    gene_list2 = benchmarking_lists[[..3]],
                                    gene_universe = ensembl_gene_ids
            )
        ),
        p_value = map(
            .x = hypergeometric_res,
            .f = ~ .x$p_value
        )
    ) %>% 
    unnest(cols = c(p_value)) %>%
    mutate(p_adj = p.adjust(p_value, method = "fdr"))



# Save objects for figures ------------------------------------------------

save(df_x_res_all, df_y_res_all, df_rcca_grcca_hypergeometric_benchmark,
     file = paste0(analysis_objects_dir, "compare_GRCCA_RCCA.Rdata")
)

